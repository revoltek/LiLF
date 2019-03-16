#!/usr/bin/env python
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

parset = lib_util.getParset()
parset_dir = parset.get('uGMRT_self','parset_dir')
sourcedb = parset.get('model','sourcedb') # relative to tgts dir "./tgts/xxx/model.txt
userReg = parset.get('model','userReg') # relative to tgts dir "./tgts/xxx/region.ref"

tgts = glob.glob('tgts/*')
MSs = lib_ms.AllMSs( glob.glob('mss/*.MS'), s )

############################################################################
# Clear
logger.info('Cleaning...')
lib_util.check_rm('img')
os.makedirs('img')
lib_util.check_rm('self')
os.makedirs('self/images')
os.makedirs('self/solutions')

# make beam
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg')
beamReg = 'self/beam.reg'

#################################################################
# Get online model
if sourcedb is None:
    if not os.path.exists('tgts.skydb'):
        fwhm = MSs.getListObj()[0].getFWHM()
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "http://172.104.228.177/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f"' % (radeg, decdeg, fwhm/2.))
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<0.1')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')
        apparent = False

    sourcedb = 'tgts.skydb'

##################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/'+sourcedb_basename)
    logger.debug('Copy: '+sourcedb+' -> '+MS)
    os.system('cp -r '+sourcedb+' '+MS)

logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA,SUBTRACTED_DATA,CORRECTED_DATA_DIE -i DATA', log="$nameMS_addcol.log", commandType="python")

logger.info('Add model to MODEL_DATA...')
MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=$pathMS/'+sourcedb_basename, log='$nameMS_pre.log', commandType='DPPP')

# Smooth DATA -> SMOOTHED_DATA
#logger.info('BL-based smoothing...')
#MSs.run('BLsmooth.py -r -f 0.2 -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType='python')

#####################################################################################################
# Self-cal cycle
for c in range(3):

    logger.info('Start selfcal cycle: '+str(c))

    # solve - concat*.MS:SMOOTHED_DATA
    logger.info('Solving G...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/gs.h5 \
             sol.solint=1 sol.nchan=1 sol.mode=complexgain sol.smoothnessconstraint=1e6', \
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')

    # LoSoTo plot
    # TODO: add some smoothing on Ga? Like normalization!
    lib_util.run_losoto(s, 'gs'+str(c), [MS+'/gs.h5' for MS in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-flag.parset'])
    os.system('mv plots-gs'+str(c)+'* self/solutions/')
    os.system('mv cal-gs'+str(c)+'*.h5 self/solutions/')

    # correct phases - MS:DATA -> MS:CORRECTED_DATA
    logger.info('Correcting Gp...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.correction=phase000', \
                log='$nameMS_corGp-c'+str(c)+'.log', commandType='DPPP')
    if c >= 2:
        # correct amplitudes - MS:DATA -> MS:CORRECTED_DATA
        logger.info('Correcting Ga...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.correction=amplitude000', \
                log='$nameMS_corGa-c'+str(c)+'.log', commandType='DPPP')

    # set image size at 1.5 * FWHM
    imgsizepix = 1.5*MSs.getListObj()[0].getFWHM()*3600/2.
    if c>=2: imgsizepix *= 2 # last cycle make a very large image to catch source in the sidelobes

    # clean mask clean
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='2arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', mgain=0.9, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=20, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl=4, atrous_do=False)
    
    # baseline averaging possible as we cut longest baselines (also it is in time, where smearing is less problematic)
    # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imgsizepix, scale='2arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', mgain=0.9, \
            #multiscale='', multiscale_scales='0,5,10,20,40', \
            baseline_averaging=5, auto_threshold=1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')

    if c != 2:

        im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
        im.makeMask(threshisl=5, atrous_do=False)
        im.selectCC()

        # predict - ms: MODEL_DATA
        logger.info('Predict model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in range(3) ]
os.system('mv img/wideM-2-sources.txt self/images')

sys.exit()

# final, large self-cal image
image_field = lib_img.Image('self/images/wideM-2-MFS-image.fits', userReg=userReg)
image_field.makeMask(threshisl=5, atrous_do=True)
image_field.selectCC()

# Move DIE-corrected data into CORRECTED_DATA_DIE
logger.info('Set CORRECTED_DATA_DIE = CORRECTED_DATA...')
MSs.run('taql "update $pathMS set CORRECTED_DATA_DIE = CORRECTED_DATA"', log='$nameMS_taql2.log', commandType='general')

logger.info('== Starting DDE cal ==')
logger.info('Cleaning...')
lib_util.check_rm('ddcal')
os.makedirs('ddcal/masks')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/skymodels')

# TODO: make sources outside primary beam withougt a "patch" and reduce the size of the facet mask

#####################################################################################################
# DDE-cal cycle
for c in range(3):

    logger.info('Start DDE-cal cycle: '+str(c))

    # Prepare facets
    os.makedirs('ddcal/masks/regions-c%02i' % c)
    mask_voro = 'ddcal/masks/facets%02i.fits' % c

    ### group into patches corresponding to the mask islands
    lsm = lsmtool.load(image_field.skymodel_cut)
    lsm.group(image_field.maskname, root='Isl')

    ### select bright sources
    # TODO: aggregate nearby sources
    lsm.select('I >= 0.2 Jy', aggregate='sum')
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
    skymodel_cl_plot = 'ddcal/skymodels/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # convert to blob
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_cl_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)

    ### select the rest of the sources to be subtracted
    lsm = lsmtool.load(image_field.skymodel_cut)
    lsm.group(image_field.maskname, root='Isl')
    lsm.select('I < 0.2 Jy', aggregate='sum')
    lsm.ungroup()
    rest_field = lsm.getColValues('I')
    rest_field = np.sum(rest_field)
    logger.info("Total flux in rest field %i Jy" % rest_field)

    # write file
    skymodel_rest = 'ddcal/skymodels/skymodel%02i_rest.txt' % c
    lsm.write(skymodel_rest, format='makesourcedb', clobber=True)
    skymodel_rest_plot = 'ddcal/skymodels/skymodel%02i_rest.png' % c
    lsm.plot(fileName=skymodel_rest_plot, labelBy='patch')

    # convert to blob
    skymodel_rest_skydb = skymodel_rest.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_rest_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_rest, skymodel_rest_skydb), log='makesourcedb_rest.log', commandType='general')
    s.run(check=True)

    ### create regions (using cluster directions)
    logger.info("Create regions.")
    lsm = lsmtool.load(image_field.skymodel_cut)
    # use the cycle=1 image that is as large as the beam
    directions_in, directions_out = lib_dd.split_directions(directions,'self/images/wideM-1-MFS-image.fits')
    lib_dd.make_voronoi_reg(directions_in, 'self/images/wideM-1-MFS-image.fits', \
        outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/skymodels/voronoi%02i.png' % c)
    # TODO: check if group ignore sources outside mask_voro
    lsm.group('facet', facet=mask_voro, root='Isl_patch')
    sizes = lib_dd.sizes_from_mask_voro(mask_voro)
    directions = lib_dd.directions_from_mask_voro(mask_voro)

    # add sizes and directions for sources outside mask
    for direction in directions_out:
        sizes[direction] = [0.1,0.1]
        directions[direction] = directions_out[direction]

    # write file
    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
    skymodel_voro_plot = 'ddcal/skymodels/skymodel%02i_voro.png' % c
    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')

    # convert to blob
    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_voro_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
    s.run(check=True)

    del lsm

    ###############################################################
    # Calibration

    #Predict - ms:MODEL_DATA
    logger.info('Add rest_field to MODEL_DATA...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb,log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    # Empty dataset from faint sources (TODO: better corrupt with DDE solutions when available before subtract)
    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA_DIE - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA_DIE - MODEL_DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')

    # Smoothing - ms:SUBTRACTED_DATA -> ms:SMOOTHED_DATA
    # TODO: check if it makes sense
    #logger.info('BL-based smoothing...')
    #MSs.run('BLsmooth.py -f 1.0 -r -i SUBTRACTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solGdd.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA sol.h5parm=$pathMS/cal-dd-c'+str(c)+'.h5 sol.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    lib_util.run_losoto(s, 'dd-c'+str(c), [MS+'/cal-dd-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-amp-dd.parset', parset_dir+'/losoto-plot-ph-dd.parset'])
    os.system('mv plots-dd-c'+str(c)+'* ddcal/plots')

    ###########################################################
    # Empty the dataset
    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA_DIE...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA_DIE"', log='$nameMS_taql2-c'+str(c)+'.log', commandType='general')

    logger.info('Subtraction...')

    for i, p in enumerate(patchNames):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+p+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p,log='$nameMS_pre1-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        # TODO: corrupt also for amplitudes?
        logger.info('Patch '+p+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-dd-c'+str(c)+'.h5 cor.direction=['+p+'] cor.correction=phase000 cor.invert=False \
                 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA', \
                log='$nameMS_corrupt1-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': subtract...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql3-c'+str(c)+'-p'+str(p)+'.log', commandType='general')

    ###########################################################
    # Facet imaging
    for i, p in enumerate(patchNames):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+p+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p, \
                   log='$nameMS_pre2-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        # TODO: corrupt also for amplitudes?
        logger.info('Patch '+p+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-dd-c'+str(c)+'.h5 cor.direction=['+p+'] cor.correction=phase000 cor.invert=False \
                 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA', \
                 log='$nameMS_corrupt2-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': add...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql4-c'+str(c)+'-p'+str(p)+'.log', commandType='general')

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+p+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-dd-c'+str(c)+'.h5 cor.direction=['+p+'] cor.correction=phase000 \
               msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED__DATA', \
               log='$nameMS_cor-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(directions[p][0])+'deg,'+str(directions[p][1])+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': imaging...')

        # set pixscale and imsize
        MSs_shift = lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s )
        pixscale = MSs_shift.getListObj()[0].getResolution()/2. # weighting lower the resolutions a bit, therefore a /2 should be enough
    
        # TODO: test uneven size
        size = np.max(sizes[p])*1.05 # add 5%
        imsize = int(size/(pixscale/3600.))
        if imsize < 512: imsize = 512
        if imsize % 2 == 1: imsize += 1 # make even
    
        logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))
    
        # clean 1
        logger.info('Cleaning ('+str(p)+')...')
        imagename = 'img/ddcal-'+str(p)
        lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs_shift.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
                weight='briggs 0.', niter=10000, no_update_model_required='', minuv_l=30, mgain=0.85, \
                baseline_averaging=5, parallel_deconvolution=256, \
                auto_threshold=20, join_channels='', fit_spectral_pol=2, channels_out=8)
    
        # make mask
        im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
        im.makeMask(threshisl = 3)
    
        # clean 2
        # TODO: add -parallel-deconvolution when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
        logger.info('Cleaning w/ mask ('+str(p)+')...')
        imagename = 'img/ddcalM-'+str(p)
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs_shift.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
                weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
                baseline_averaging=5, auto_threshold=0.1, fits_mask=im.maskname, \
                join_channels='', fit_spectral_pol=2, channels_out=8)
    
        os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')

    ##############################################################
    # Mosaiching
    images = []
    for patchName in patchNames:
        image = lib_img.Image('img/ddcalM-%s-MFS-image.fits' % patchName, userReg = userReg)
        image.selectCC()
        # restrict skymodel to facet
        lsm = lsmtool.load(image.skymodel_cut)
        lsm.select('%s == %i' % ( mask_voro, int(patchName[10:]) ))
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
    image_field = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)
    # TODO: how to keep bright sources in the sidelobes into the model?

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise

logger.info("Done.")
