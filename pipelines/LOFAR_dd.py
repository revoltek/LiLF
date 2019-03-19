#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-dd.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_dd','parset_dir')
maxniter = parset.getint('LOFAR_dd','maxniter')
userReg = parset.get('model','userReg')

####################################################
MSs_self = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs_self.getListObj()[0].getPhaseCentre()
MSs_self.getListObj()[0].makeBeamReg('self/beam.reg', to_null=True) # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null
beamReg = 'self/beam.reg'

##########################
logger.info('Cleaning...')
lib_util.check_rm('ddcal')
os.makedirs('ddcal/masks')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/skymodels')

def clean(p, MSs, size, res='normal', apply_beam=False):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.getListObj()[0].getResolution()/2. # weighting lower the resolutions a bit, therefore a /2 should be enough

    # TODO: test uneven size
    size = np.max(size)*1.05 # add 5%
    imsize = int(size/(pixscale/3600.))
    if imsize % 2 != 0: imsize += 1 # make even
    if imsize < 512: imsize = 512 # prevent supersmall images

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    if res == 'normal':
        weight = 'briggs 0'
        maxuv_l = 1e30
    elif res == 'high':
        weight = 'uniform'
        maxuv_l = 1e30
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 4000

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/ddcal-'+str(p)
    lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=256, \
            auto_threshold=20, join_channels='', fit_spectral_pol=2, channels_out=8)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 3)

    # clean 2
    # TODO: add -parallel-deconvolution when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/ddcalM-'+str(p)
    if apply_beam:
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            auto_threshold=0.1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)
    else:
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
            baseline_averaging=5, auto_threshold=0.1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')


############################################################
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs_self.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg, beamReg = beamReg)
mosaic_image.selectCC()
rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    lib_util.check_rm('img')
    os.makedirs('img')
    os.makedirs('ddcal/masks/regions-c%02i' % c)
    mask_voro = 'ddcal/masks/facets%02i.fits' % c

    ### group into patches corresponding to the mask islands
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group(mosaic_image.maskname, root='Isl')

    ### select bright sources
    lsm.select('I >= 2.0 Jy', aggregate='sum')
    lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
    directions = lsm.getPatchPositions()
    patchNames = lsm.getPatchNames()
    patchFluxes = lsm.getColValues('I', aggregate='sum')

    logger.info("Created %i bright sources" % len(directions))
    logger.debug("Island info:")
    for i, patchname in enumerate(patchNames):
        logger.debug("%s: Flux=%f (coord: %s)" % ( patchname, patchFluxes[i], str(directions[patchname]) ) )
    
    tot_flux = np.sum(patchFluxes)
    logger.info("Total flux of bright sources %i Jy" % tot_flux)
    
    # write file
    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
    skymodel_cl_plot = 'ddcal/masks/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # convert to blob
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_cl_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)
    
    ### select the rest of the sources to be subtracted
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group(mosaic_image.maskname, root='Isl')
    lsm.select('I < 2.0 Jy', aggregate='sum')
    lsm.ungroup()
    rest_field = lsm.getColValues('I')
    rest_field = np.sum(rest_field)
    logger.info("Total flux in rest field %i Jy" % rest_field)

    # write file
    skymodel_rest = 'ddcal/masks/skymodel%02i_rest.txt' % c
    lsm.write(skymodel_rest, format='makesourcedb', clobber=True)
    skymodel_rest_plot = 'ddcal/masks/skymodel%02i_rest.png' % c
    lsm.plot(fileName=skymodel_rest_plot, labelBy='patch')
       
    # convert to blob
    skymodel_rest_skydb = skymodel_rest.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_rest_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_rest, skymodel_rest_skydb), log='makesourcedb_rest.log', commandType='general')
    s.run(check=True)
    
    ### create regions (using cluster directions)
    logger.info("Create regions.")
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/masks/voronoi%02i.png' % c)
    lsm.group('facet', facet=mask_voro, root='Isl_patch')
    sizes = lib_dd.sizes_from_mask_voro(mask_voro)
    directions = lib_dd.directions_from_mask_voro(mask_voro)

    # write file
    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
    skymodel_voro_plot = 'ddcal/masks/skymodel%02i_voro.png' % c
    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')

    # convert to blob
    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_voro_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
    s.run(check=True)

    del lsm
    ################################################################

    #Predict - ms:MODEL_DATA
    logger.info('Add rest_field to MODEL_DATA...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb,log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    # Empty dataset from faint sources (TODO: better corrupt with DDE solutions when available before subtract)
    logger.info('Set SUBTRACTED_DATA = DATA - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')

    # Smoothing - ms:SUBTRACTED_DATA -> ms:SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -f 1.0 -r -i SUBTRACTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')    

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solDD.parset msin=$pathMS ddecal.h5parm=$pathMS/cal-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')
    logger.info('Calibrating2...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solDDg.parset msin=$pathMS ddecal.h5parm=$pathMS/calg-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solDDg-c'+str(c)+'.log', commandType='DPPP')
    logger.info('Calibrating3...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solDDg.parset msin=$pathMS ddecal.h5parm=$pathMS/calgsmooth-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb+" ddecal.smoothnessconstraint=2e6", \
            log='$nameMS_solDDg-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    lib_util.run_losoto(s, 'c'+str(c), [MS+'/cal-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot.parset'])
    os.system('mv plots-c'+str(c)+'* ddcal/plots')

    ##############################################################
    # low S/N DIE corrections
    # TODO: add amp and FR sol + correction here after ft() a DDE-corrupted model

    ###########################################################
    # Empty the dataset
    logger.info('Set SUBTRACTED_DATA = DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql2-c'+str(c)+'.log', commandType='general')

    logger.info('Subtraction...')
    ##MSs.run('DPPP '+parset_dir+'/DPPP-sub.parset msin=$pathMS sub.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 sub.sourcedb='+skymodel_voro_skydb, \
    ##               log='$nameMS_sub-c'+str(c)+'.log', commandType='DPPP')

    for i, p in enumerate(patchNames):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+p+': predict...')
        #pre.applycal.h5parm='+ms+'/cal-c'+str(c)+'.h5 pre.applycal.direction='+p, \
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p,log='$nameMS_pre1-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
    
        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+p+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+']', \
                log='$nameMS_corrupt1-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
        
        logger.info('Patch '+p+': subtract...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql3-c'+str(c)+'-p'+str(p)+'.log', commandType='general')

    ##  for patch, phasecentre in directions.iteritems():
    ##      # add back single path - ms:SUBTRACTED_DATA -> ms:CORRECTED_DATA
    ##      logger.info('Patch '+patch+': add back...')
    ##      MSs.run('DPPP '+parset_dir+'/DPPP-add.parset msin=$pathMS add.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 add.sourcedb='+skymodel_voro_skydb+' add.directions=[['+patch+']]', \
    ##          log='$nameMS_add-c'+str(c)+'-p'+str(patch)+'.log', commandType='DPPP')

    # Add back 
    for i, p in enumerate(patchNames):
        #TODO: see if we can phase shift and average before predict-corrupt=correct
        # predict - ms:MODEL_DATA
        logger.info('Patch '+p+': predict...')
        #pre.applycal.h5parm='+ms+'/cal-c'+str(c)+'.h5 pre.applycal.direction='+p, \
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p, \
                   log='$nameMS_pre2-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+p+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+']', \
                 log='$nameMS_corrupt2-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': add...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql4-c'+str(c)+'-p'+str(p)+'.log', commandType='general')

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+p+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+']', \
               log='$nameMS_cor-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(directions[p][0])+'deg,'+str(directions[p][1])+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
        
        logger.info('Patch '+p+': imaging...')
        clean(p, lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=sizes[p], apply_beam = c==maxniter )

        # if one wants to make a low-res pathc
        if p == 'Isl_patch_663': 
            logger.info('Patch '+p+': imaging high-res...')
            clean(p+'high', lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=sizes[p], res='high')
            logger.info('Patch '+p+': predict high-res...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=img/ddcalM-'+p+'high-MFS-image.fits pre.sources='+p,log='$nameMS_pre1-c'+str(c)+'-p'+p+'.log', commandType='DPPP')
            logger.info('Patch '+p+': subtract high-res...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql4-c'+str(c)+'-p'+p+'.log', commandType='general')
            logger.info('Patch '+p+': imaging low-res...')
            clean(p+'low', lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=sizes[p], res='low', apply_beam = c==maxniter )

    ##############################################################
    # Mosaiching
    images = []
    for patchName in patchNames:
        image = lib_img.Image('img/ddcalM-%s-MFS-image.fits' % patchName, userReg = userReg)
        image.selectCC()
        # restrict skymodel to facet
        lsm = lsmtool.load(image.skymodel_cut)
        lsm.group('facet', facet=mask_voro )
        lsm.select('PatchName == %i' % int(patchName[10:]) )
        lsm.write(image.skymodel_cut, format='makesourcedb', clobber=True)
        images.append(image)

    logger.info('Mosaic: image...')
    image_files = ' '.join([image.imagename for image in images])
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    # prepare new skymodel
    lsm = lsmtool.load(images[0].skymodel_cut)
    for image in images[1:]:
        lsm2 = lsmtool.load(image.skymodel_cut)
        lsm.concatenate(lsm2)
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)

    os.system('cp img/*M*MFS-image.fits img/mos-MFS-image.fits img/mos-MFS-residual.fits ddcal/images/c%02i' % c )
    mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise
