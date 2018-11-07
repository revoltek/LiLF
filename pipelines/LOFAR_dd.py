#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
lib_log.set_logger('pipeline-dd.logger')
logger = lib_log.logger
s = lib_util.Scheduler(dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('dd','parset_dir')
maxniter = parset.getint('dd','maxniter')
userReg = parset.get('model','userReg')

####################################################
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg', to_null=True) # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null
beamReg = 'self/beam.reg'

##########################
logger.info('Cleaning...')
#lib_util.check_rm('ddcal')
#os.makedirs('ddcal/masks')
#os.makedirs('ddcal/plots')
#os.makedirs('ddcal/images')
#os.makedirs('ddcal/skymodels')

def clean(p, MSs, size=2.):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.getListObj()[0].getResolution()/3.
    imsize = int(size/(pixscale/3600.))

    if imsize < 512:
        imsize = 512

    if imsize % 2 == 1: imsize += 1 # make even

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/ddcal-'+str(p)
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -auto-threshold 20 \
            -join-channels -fit-spectral-pol 2 -channels-out 10 '+MSs.getStrWsclean(), \
            log='wsclean-'+str(p)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 3)

    # clean 2
    #-multiscale -multiscale-scale-bias 0.5 \
    #-auto-mask 3 -rms-background-window 40 -rms-background-method rms-with-min \
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/ddcalM-'+str(p)
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0.0 -niter 1000000 -no-update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -auto-threshold 0.1 -fits-mask '+im.maskname+' \
            -join-channels -fit-spectral-pol 2 -channels-out 10 -save-source-list '+MSs.getStrWsclean(), \
            log='wscleanM-'+str(p)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-'+str(p)+'.log | grep "background noise"')

    return imagename


#############################################################
#logger.info('Copy data...')
#if not os.path.exists('mss-dd'):
#    os.makedirs('mss-dd')
#    MSs.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
#                log='$nameMS_avg.log', commandType='DPPP')
#MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
#       
#logger.info('Add columns...')
#MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA', log='$nameMS_addcol.log', commandType='python')
#
###############################################################
#logger.info('BL-based smoothing...')
#MSs.run('BLsmooth.py -f 1.0 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')
#
## setup initial model
#mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg)
#mosaic_image.selectCC()
#rms_noise_pre = np.inf

for c in xrange(maxniter):
#    logger.info('Starting cycle: %i' % c)
#
#    lib_util.check_rm('img')
#    os.makedirs('img')
#    os.makedirs('ddcal/images/c%02i/regions' % c)
#    mask = 'ddcal/masks/facets%02i.fits' % c
#
#    ### group into patches of similar flux
#    lsm = lsmtool.load(mosaic_image.skymodel_cut)
#    lsm.group('tessellate', targetFlux='20Jy', root='Dir', applyBeam=False, method = 'wmean', pad_index=True)
#    lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
#    directions = lsm.getPatchPositions()
#    patches = lsm.getPatchNames()
#    logger.info("Created %i directions." % len(patches))
#
#    # write file
#    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
#    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
#    skymodel_cl_plot = 'ddcal/skymodels/skymodel%02i_cluster.png' % c
#    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')
#
#    # convert to blob
#    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
#    lib_util.check_rm(skymodel_cl_skydb)
#    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
#    s.run(check=True)
#
#    ### create regions (using cluster directions)
#    logger.info("Create regions.")
#    lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, outdir_reg='ddcal/images/c%02i/regions/' % c, out_mask=mask, png='ddcal/skymodels/voronoi%02i.png' % c)
#    lsm.group('facet', facet=mask, root='Dir')
#    lsm.setPatchPositions(method='mid') # recalculate the patch centre as mid point for imaging
#    directions = lsm.getPatchPositions()
#    sizes = lsm.getPatchSizes(units='degree')
#
#    # write file
#    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
#    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
#    skymodel_voro_plot = 'ddcal/skymodels/skymodel%02i_voro.png' % c
#    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')
#
#    # convert to blob
#    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
#    lib_util.check_rm(skymodel_voro_skydb)
#    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
#    s.run(check=True)
#
#    del lsm
#
#    ################################################################
#    # Calibration
#    logger.info('Calibrating...')
#    MSs.run('DPPP '+parset_dir+'/DPPP-solDD.parset msin=$pathMS ddecal.h5parm=$pathMS/cal-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
#            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')
#
#    # Plot solutions
#    lib_util.run_losoto(s, 'c'+str(c), [MS+'/cal-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot.parset'])
#    os.system('mv plots-c'+str(c)+'* ddcal/plots')
#
#   ###########################################################
#   # Empty the dataset
#    logger.info('Set SUBTRACTED_DATA = DATA...')
#    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')
#
#    logger.info('Subtraction...')
#    MSs.run('DPPP '+parset_dir+'/DPPP-sub.parset msin=$pathMS sub.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 sub.sourcedb='+skymodel_voro_skydb, \
#                   log='$nameMS_sub-c'+str(c)+'.log', commandType='DPPP')
#
#    ## TODO: test
#    #logger.info('Empty imaging')
#    #s.add('wsclean -reorder -name img/testSUB -size 5000 5000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
#    #        -scale 10arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.9 \
#    #        -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 20 -minuv-l 30 -data-column SUBTRACTED_DATA '+MSs.getStrWsclean(), \
#    #        log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors='max')
#    #s.run(check=True)
#
#    for i, p in enumerate(patches):
#
#        # add back single path - ms:SUBTRACTED_DATA -> ms:CORRECTED_DATA
#        logger.info('Patch '+p+': add back...')
#        MSs.run('DPPP '+parset_dir+'/DPPP-add.parset msin=$pathMS add.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 add.sourcedb='+skymodel_voro_skydb+' add.directions=[['+p+']]', \
#                   log='$nameMS_add-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#
#        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
#        logger.info('Patch '+p+': correct...')
#        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+']', \
#               log='$nameMS_cor-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#
#        logger.info('Patch '+p+': phase shift and avg...')
#        lib_util.check_rm('mss-dir')
#        os.makedirs('mss-dir')
#        phasecentre = directions[p]
#        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
#                shift.phasecenter=['+str(phasecentre[0].degree)+'deg,'+str(phasecentre[1].degree)+'deg\]', \
#                log='$nameMS_shift-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#        
#        logger.info('Patch '+p+': imaging...')
#        clean(p, lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=sizes[i])

    ##############################################################
    # Mosaiching

    images = []
    for image, region in zip( sorted(glob.glob('img/ddcalM-Dir*MFS-image.fits')), sorted(glob.glob('ddcal/images/c%02i/regions/Dir*' % c)) ):
        images.append( lib_img.Image(image, facetReg = region, userReg = userReg) )

    logger.info('Mosaic: image...')
    image_files = ' '.join([image.imagename for image in images])
    print image_files
    sys.exit()
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    logger.info('Mosaic: residuals...')
    image_files = ' '.join([image.imagename.replace('image', 'residual') for image in images])
    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    # prepare new skymodel
    skymodels = []
    for image in images:
        image.select_cc()
        skymodels.append(image.skymodel_cut)
    lsm = lsmtool.load(skymodels[0])
    for skymodel in skymodels[1:]:
        lsm2 = lsmtool.load(skymodel)
        lsm.concatenate(lsm2)
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)

    os.system('cp img/*M*MFS-image.fits img/mos-MFS-image.fits img/mos-MFS-residual.fits ddcal/images/c%02i' % c )
    mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise
