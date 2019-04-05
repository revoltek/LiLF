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
#MSs_self.getListObj()[0].makeBeamReg('self/beam.reg', to_null=True) # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null
#beamReg = 'self/beam.reg'
fwhm = MSs_self.getListObj()[0].getFWHM(freq='min')

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
    pixscale = MSs.getListObj()[0].getResolution() 
    # weighting lower the resolutions a bit, therefore a /2 should be enough
    if res == 'normal':
        pixscale /= 2.
    elif res == 'high':
        pixscale /= 4.
    elif res == 'low':
        pixscale /= 1. # no change

    imsize = [0,0]
    imsize[0] = int(size[0]*1.05/(pixscale/3600.)) # add 5%
    imsize[1] = int(size[1]*1.05/(pixscale/3600.)) # add 5%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 64: imsize[0] == 64
    if imsize[1] < 64: imsize[1] == 64

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    if res == 'normal':
        weight = 'briggs 0'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.7'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500

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
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imsize, scale=str(pixscale)+'arcsec', pol='IQUV', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            auto_threshold=0.1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)
    else:
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            baseline_averaging=5, auto_threshold=0.1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')


############################################################
# TEST: use SUBTRACTED_DATA (no pre-correction) or CORRECTED_DATA (DIE iono correction)?
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs_self.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg)
mosaic_image.selectCC()
# TEST:
#mosaic_image = lib_img.Image('ddcal/images/c00/mos-MFS-image.fits', userReg = userReg)
rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)
    directions = []

    lib_util.check_rm('img')
    os.makedirs('img')
    if not os.path.exists('ddcal/masks/regions-c%02i' % c): os.makedirs('ddcal/masks/regions-c%02i' % c)
    if not os.path.exists('ddcal/images/c%02i' % c): os.makedirs('ddcal/images/c%02i' % c)
    mask_voro = 'ddcal/masks/facets%02i.fits' % c

    ### TTESTTESTTEST: DIE image
    if c == 0:
        clean('init', MSs, size=(fwhm,fwhm), res='normal')
    ############################

    ### group into patches corresponding to the mask islands
    # TODO: aggregate nearby sources. Expand mask?
    mask_cl = mosaic_image.imagename.replace('MFS-image.fits', 'mask-cl.fits')
    # this mask is with no user region, done isolate only bight compact sources
    lib_img.make_mask.make_mask(image_name=mosaic_image.imagename, mask_name=mask_cl, threshisl=5)
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group(mask_cl, root='Isl')

    ### select bright sources
    lsm.select('I >= 2.0 Jy', aggregate='sum')
    lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
    for name, flux in zip(lsm.getPatchNames(), lsm.getColValues('I', aggregate='sum')):
        direction = lib_dd.Direction(name)
        position = [ lsm.getPatchPositions()[name][0].deg, lsm.getPatchPositions()[name][1].deg ]
        direction.set_position( position, cal=True )
        direction.set_flux(flux, cal=True)
        directions.append(direction)

    logger.info("Created %i bright sources" % len(directions))
    tot_flux = np.sum([d.flux_cal for d in directions])
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
    del lsm
    
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
    lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, \
            outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/masks/voronoi%02i.png' % c)
    lsm.group('facet', facet=mask_voro, root='Isl_patch')
    [ d.add_mask_voro(mask_voro) for d in directions ]

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

    logger.debug("Islands' info:")
    for i, d in enumerate(directions):
        logger.info("%s: Flux=%f (coord: %s - size: %s deg)" % ( d.name, d.flux_cal, str(d.position_cal), str(d.size) ) )

    ################################################################

    #Predict - ms:MODEL_DATA
    logger.info('Add rest_field to MODEL_DATA...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb,log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    # Empty dataset from faint sources (TODO: better corrupt with DDE solutions when available before subtract)
    logger.info('Set SUBTRACTED_DATA = DATA - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

    ### TESTTESTTEST: empty image with cals
    MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    clean('onlycals-c'+str(c), MSs, size=(fwhm,fwhm), res='normal')
    ########################################

    # Smoothing - ms:SUBTRACTED_DATA -> ms:SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -f 1.0 -r -i SUBTRACTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')    

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solDD.parset msin=$pathMS ddecal.h5parm=$pathMS/cal-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    lib_util.run_losoto(s, 'c'+str(c), [MS+'/cal-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot.parset'])
    os.system('mv plots-c'+str(c)+'* ddcal/plots')

    ##############################################################
    # low S/N DIE corrections
    if c>0:
        logger.info('DIE calibration...')
        # predict and corrupt each facet
        logger.info('Reset MODEL_DATA...')
        MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        for i, d in enumerate(directions):
            # predict - ms:MODEL_DATA
            logger.info('Patch '+d.name+': predict+corrupt...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.operation=add pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name+ \
                    'pre.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 pre.applycal.direction=['+d.name+'] pre.applycal.correction=tec000', \
                log='$nameMS_preDIE-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        # Smoothing - ms:DATA -> ms:SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -f 1.0 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')    
    
        # DIE Calibration - ms:SMOOTHED_DATA
        logger.info('Calibrating DIE...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solDDg.parset msin=$pathMS ddecal.h5parm=$pathMS/calG-c'+str(c)+'.h5', \
                log='$nameMS_solDDg-c'+str(c)+'.log', commandType='DPPP')
    
        # Plot solutions
        lib_util.run_losoto(s, 'G-c'+str(c), [MS+'/calG-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset'])
        os.system('mv plots-G-c'+str(c)+'* ddcal/plots')

    ###########################################################
    ## Empty the dataset
    logger.info('Set SUBTRACTED_DATA = DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

    logger.info('Subtraction...')
    ## TODO: use this command once DP3 is fixed
    ##MSs.run('DPPP '+parset_dir+'/DPPP-sub.parset msin=$pathMS sub.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 sub.sourcedb='+skymodel_voro_skydb, \
    ##               log='$nameMS_sub-c'+str(c)+'.log', commandType='DPPP')

    for i, d in enumerate(directions):
        
        # TODO: use this command once tested
        #MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=SUBTRACTED_DATA pre.operation=sub pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name+ \
        #            'pre.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 pre.applycal.direction=['+d.name+']', \
        #            log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d.name+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name,log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
    
        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+d.name+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS cor.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor.direction=['+d.name+']', \
                log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
        logger.info('Patch '+d.name+': subtract...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

    ### TESTTESTTEST: empty image
    MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    clean('empty-c'+str(c), MSs, size=(fwhm,fwhm), res='normal')
    ##############################

    # Add back 
    logger.info('Facet imaging...')
    ##  for patch, phasecentre in directions.iteritems():
    ##      # add back single path - ms:SUBTRACTED_DATA -> ms:CORRECTED_DATA
    ##      logger.info('Patch '+patch+': add back...')
    ##      MSs.run('DPPP '+parset_dir+'/DPPP-add.parset msin=$pathMS add.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 add.sourcedb='+skymodel_voro_skydb+' add.directions=[['+patch+']]', \
    ##          log='$nameMS_add-c'+str(c)+'-p'+str(patch)+'.log', commandType='DPPP')
    for i, d in enumerate(directions):
        #TODO: see if we can phase shift and average before predict-corrupt=correct
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d.name+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                   log='$nameMS_pre2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+d.name+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS cor.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor.direction=['+d.name+']', \
                 log='$nameMS_corrupt2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': add...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+d.name+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor.direction=['+d.name+']', \
               log='$nameMS_cor-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(d.position_facet[0])+'deg,'+str(d.position_facet[1])+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
        logger.info('Patch '+d.name+': imaging...')
        clean(d.name, lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size, apply_beam = c==maxniter )

        # TEST: if one wants to make a low-res patch
        if c>1:
            logger.info('Patch '+d.name+': imaging high-res...')
            clean(d.name+'-high', lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size, res='high')
            logger.info('Patch '+d.name+': predict high-res...')
            # predict - ms:MODEL_DATA
            s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % \
                    ('img/ddcalM-'+d.name+'-high-sources.txt', 'img/ddcalM-'+d.name+'-high-sources.skydb'), log='makesourcedb_'+d.name+'.log', commandType='general' )
            s.run(check=True)
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=img/ddcalM-'+d.name+'-high-sources.skydb pre.sources='+d.name, \
                    log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
            logger.info('Patch '+d.name+': subtract high-res...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')
            logger.info('Patch '+d.name+': imaging low-res...')
            clean(d.name+'-low', lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size, res='low', apply_beam = c==maxniter )

    ##############################################################
    # Mosaiching

    # reorder in increasing isl_num order
    isl_nums = [d.isl_num for d in directions]
    directions = [d for _, d in sorted(zip(isl_nums,directions))]

    for d in directions:
        d.image = lib_img.Image('img/ddcalM-%s-MFS-image.fits' % d.name, userReg = userReg)
        d.image_res = lib_img.Image('img/ddcalM-%s-MFS-residual.fits' % d.name, userReg = userReg)
        d.image_low = lib_img.Image('img/ddcalM-%s-low-MFS-image.fits' % d.name, userReg = userReg)
        d.image_high = lib_img.Image('img/ddcalM-%s-high-MFS-image.fits' % d.name, userReg = userReg)

        # restrict skymodel to facet
        d.image.selectCC()
        lsm = lsmtool.load(d.image.skymodel_cut)
        lsm.group('facet', facet=mask_voro, root='Isl_patch' )
        lsm.select('Patch = Isl_patch_%i' % d.isl_num )
        lsm.write(d.image.skymodel_cut, format='makesourcedb', clobber=True)

    logger.info('Mosaic: image...')
    image_files = ' '.join([d.image.imagename for d in directions])
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    logger.info('Mosaic: residual image...')
    image_files = ' '.join([d.image_res.imagename for d in directions])
    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    if c>1:
        logger.info('Mosaic: low-res image...')
        image_files = ' '.join([d.image_low.imagename for d in directions])
        mosaic_residual = 'img/mos-low-MFS-image.fits'
        s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual, log='mosaic-img-low-c'+str(c)+'.log', commandType='python')
        s.run(check=True)
    
        logger.info('Mosaic: high-res image...')
        image_files = ' '.join([d.image_high.imagename for d in directions])
        mosaic_residual = 'img/mos-high-MFS-image.fits'
        s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual, log='mosaic-img-high-c'+str(c)+'.log', commandType='python')
        s.run(check=True)

    # prepare new skymodel
    lsm = lsmtool.load(directions[0].image.skymodel_cut)
    for image in [d.image for d in directions[1:]]:
        lsm2 = lsmtool.load(image.skymodel_cut)
        lsm.concatenate(lsm2)
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)
    del lsm

    os.system('cp img/*M*MFS-image.fits img/mos*.fits ddcal/images/c%02i' % c )
    mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)
    mosaic_image.makeMask(threshisl=3, atrous_do=True)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > rms_noise_pre: break
    rms_noise_pre = rms_noise
