#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

parset_dir = '/home/fdg/scripts/autocal/parset_dd'
maxniter = 10 # max iteration if not converged

##########################################################################################

import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *
import pyrap.tables as pt
from make_mask import make_mask
import lsmtool

logger = set_logger('pipeline-dd.logger')
check_rm('logs')
s = Scheduler(dry=False)
mss = sorted(glob.glob('mss/TC*[0-9].MS'))
phasecentre = get_phase_centre(mss[0])
check_rm('ddcal')
os.makedirs('ddcal/regions')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/skymodels')
os.makedirs('logs/mss')

# set user masks
if 'tooth' in os.getcwd():
    user_mask = '/home/fdg/scripts/autocal/regions/tooth.reg'
else:
    user_mask = None

def clean(c, mss, size=2.):
    """
    c = cycle/name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = scale_from_ms(mss[0])
    imsize = int(size/(pixscale/3600.)*1.3)

    if imsize < 512:
        imsize = 512

    trim = int(imsize*0.9)

    if imsize % 2 == 1: imsize += 1 # make even
    if trim % 2 == 1: trim += 1 # make even

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    # clean 1
    logger.info('Cleaning ('+str(c)+')...')
    imagename = 'img/ddcal-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -trim '+str(trim)+' '+str(trim)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0. -niter 100000 -no-update-model-required -mgain 0.9 -pol I \
            -joinchannels -fit-spectral-pol 2 -channelsout 10 \
            -auto-threshold 20 -minuv-l 30 '+' '.join(mss), \
            log='wsclean-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # make mask
    im = Image(imagename+'-MFS-image.fits', user_mask=user_mask)
    im.make_mask(threshisl = 3)

    # clean 2
    #-multiscale -multiscale-scale-bias 0.5 \
    #-auto-mask 3 -rms-background-window 40 -rms-background-method rms-with-min \
    logger.info('Cleaning w/ mask ('+str(c)+')...')
    imagename = 'img/ddcalM-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -trim '+str(trim)+' '+str(trim)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0. -niter 1000000 -no-update-model-required -mgain 0.8 -pol I \
            -joinchannels -fit-spectral-pol 2 -channelsout 10 \
            -auto-threshold 0.1 -save-source-list -minuv-l 30 -fitsmask '+im.maskname+' '+' '.join(mss), \
            log='wscleanM-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-c'+str(c)+'.log | grep "background noise"')

    return imagename


############################################################
logger.info('Copy data...')
for ms in mss:
    msout = ms.replace('.MS','-cp.MS')
    if os.path.exists(msout): continue
    s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log=msout.split('/')[-1]+'_cp.log', cmd_type='NDPPP')
s.run(check=True)
mss = sorted(glob.glob('mss/TC*-cp.MS'))
       
logger.info('Add columns...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
s.run(check=True)

###############################################################
logger.info('BL-based smoothing...')
for ms in mss:
    s.add('BLsmooth.py -f 1.0 -r -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True, max_threads=6)

mosaic_image = Image(sorted(glob.glob('self/images/wide*-[0-9]-MFS-image.fits'))[-1], user_mask = user_mask)
mosaic_image.select_cc()
rms_noise_pre = np.inf

for c in xrange(maxniter):
    logger.info('Starting cycle: %i' % c)

    check_rm('img')
    os.makedirs('img')

    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group('tessellate', targetFlux='10Jy', root='Dir', applyBeam=False, method = 'wmean', pad_index=True)
    directions_clusters = lsm.getPatchPositions()
    patches = lsm.getPatchNames()
    logger.info("Created %i directions." % len(patches))

    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
    skymodel_cl_plot = 'ddcal/skymodels/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # voronoi tessellation of skymodel for imaging
    lsm.group('voronoi', root='Dir', applyBeam=False, method='mid')
    directions_shifts = lsm.getPatchPositions()
    sizes = lsm.getPatchSizes(units='degree')

    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
    skymodel_voro_plot = 'ddcal/skymodels/skymodel%02i_voro.png' % c
    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')
    del lsm

    # create regions (using cluster directions)
    make_voronoi_reg(directions_clusters, mosaic_image.imagename, outdir='ddcal/regions/', beam_reg='', png='ddcal/skymodels/voronoi%02i.png' % c)

    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    check_rm(skymodel_cl_skydb)
    s.add('run_env.sh makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', cmd_type='general' )
    s.run(check=True)

    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
    check_rm(skymodel_voro_skydb)
    s.add('run_env.sh makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', cmd_type='general')
    s.run(check=True)

    ################################################################
    # Calibration
    logger.info('Calibrating...')
    for ms in mss:
        check_rm(ms+'/cal-c'+str(c)+'.h5')
        s.add('run_env.sh NDPPP '+parset_dir+'/NDPPP-solDD.parset msin='+ms+' ddecal.h5parm='+ms+'/cal-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
                log=ms+'_solDD-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # Plot solutions
    # TODO: concat h5parm into a single file
    logger.info('Running losoto...')
    for i, ms in enumerate(mss):
        s.add('losoto -v '+ms+'/cal-c'+str(c)+'.h5 '+parset_dir+'/losoto-plot.parset', log=ms+'_losoto-c'+str(c)+'.log', cmd_type='python', processors='max')
        s.run(check=True)
        os.system('mv plots ddcal/plots/plots-c'+str(c)+'-t'+str(i))

    ############################################################
    # Empty the dataset
    logger.info('Set SUBTRACTED_DATA = DATA...')
    for ms in mss:
        s.add('taql "update '+ms+' set SUBTRACTED_DATA = DATA"', log=ms+'_taql1-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    logger.info('Subtraction...')
    for i, p in enumerate(patches):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+p+': predict...')
        #pre.applycal.h5parm='+ms+'/cal-c'+str(c)+'.h5 pre.applycal.direction='+p, \
        for ms in mss:
            s.add('run_env.sh NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p, \
                   log=ms+'_pre1-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+p+': corrupt...')
        for ms in mss:
            s.add('run_env.sh NDPPP '+parset_dir+'/NDPPP-corrupt.parset msin='+ms+' cor1.parmdb='+ms+'/cal-c'+str(c)+'.h5 cor1.direction=['+p+'] cor2.parmdb='+ms+'/cal-c'+str(c)+'.h5 cor2.direction=['+p+']', \
                 log=ms+'_corrupt1-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        logger.info('Patch '+p+': subtract...')
        for ms in mss:
            s.add('taql "update '+ms+' set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log=ms+'_taql2-c'+str(c)+'-p'+str(p)+'.log', cmd_type='general')
        s.run(check=True)

    ##############################################################
    # Imaging
    logger.info('Imaging...')

    ## TODO: test
    #logger.info('Empty imaging')
    #s.add('wsclean -reorder -name img/empty-c'+str(c)+' -datacolumn SUBTRACTED_DATA -size 3000 3000 \
    #        -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
    #        -scale 10arcsec -weight briggs 0.0 -niter 0 -no-update-model-required -mgain 1 -pol I '+' '.join(mss), \
    #        log='wscleanEmpty-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    #s.run(check=True)

    for i, p in enumerate(patches):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+p+': predict...')
        #pre.applycal.h5parm='+ms+'/cal-c'+str(c)+'.h5 pre.applycal.direction='+p, \
        for ms in mss:
            s.add('run_env.sh NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p, \
                   log=ms+'_pre2-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+p+': corrupt...')
        for ms in mss:
            s.add('run_env.sh NDPPP '+parset_dir+'/NDPPP-corrupt.parset msin='+ms+' cor1.parmdb='+ms+'/cal-c'+str(c)+'.h5 cor1.direction=['+p+'] cor2.parmdb='+ms+'/cal-c'+str(c)+'.h5 cor2.direction=['+p+']', \
                 log=ms+'_corrupt2-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        logger.info('Patch '+p+': add...')
        for ms in mss:
            s.add('taql "update '+ms+' set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log=ms+'_taql2-c'+str(c)+'-p'+str(p)+'.log', cmd_type='general')
        s.run(check=True)

        ### TEST
        #logger.info('Patch '+p+': phase shift and avg...')
        #check_rm('mss_dd')
        #os.makedirs('mss_dd')
        #for ms in mss:
        #    msout = 'mss_dd/'+os.path.basename(ms)
        #    phasecentre = directions_shifts[p]
        #    s.add('NDPPP '+parset_dir+'/NDPPP-shiftavg.parset msin='+ms+' msout='+msout+' shift.phasecenter=['+str(phasecentre[0].degree)+'deg,'+str(phasecentre[1].degree)+'deg\]', \
        #        log=ms+'_shift-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        #s.run(check=True)
        #logger.info('Patch '+p+': corrupted image...')
        #clean('uncor-'+p, glob.glob('mss_dd/*MS'), size=sizes[i])
        ### end TEST

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+p+': correct...')
        for ms in mss:
            s.add('run_env.sh NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor1.parmdb='+ms+'/cal-c'+str(c)+'.h5 cor1.direction=['+p+'] cor2.parmdb='+ms+'/cal-c'+str(c)+'.h5 cor2.direction=['+p+']', \
               log=ms+'_cor-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        logger.info('Patch '+p+': phase shift and avg...')
        check_rm('mss_dd')
        os.makedirs('mss_dd')
        for ms in mss:
            msout = 'mss_dd/'+os.path.basename(ms)
            phasecentre = directions_shifts[p]
            s.add('NDPPP '+parset_dir+'/NDPPP-shiftavg.parset msin='+ms+' msout='+msout+' shift.phasecenter=['+str(phasecentre[0].degree)+'deg,'+str(phasecentre[1].degree)+'deg\]', \
                log=ms+'_shift-c'+str(c)+'-p'+str(p)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        logger.info('Patch '+p+': imaging...')
        clean(p, glob.glob('mss_dd/*MS'), size=sizes[i])

    ##############################################################
    # Mosaiching

    os.makedirs('ddcal/images/c'+str(c))
    directions = []
    for image, region in zip( sorted(glob.glob('img/ddcalM-Dir*MFS-image.fits')), sorted(glob.glob('ddcal/regions/Dir*')) ):
        directions.append( Image(image, region_facet = region, user_mask = user_mask) )

    logger.info('Mosaic: image...')
    images = ' '.join([image.imagename for image in directions])
    masks = ' '.join([image.region_facet for image in directions])
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+images+' --masks '+masks+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', cmd_type='python')
    s.run(check=True)

    logger.info('Mosaic: residuals...')
    images = ' '.join([image.imagename.replace('image', 'residual') for image in directions])
    masks = ' '.join([image.region_facet for image in directions])
    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+images+' --masks '+masks+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', cmd_type='python')
    s.run(check=True)

    # prepare new skymodel
    skymodels = []
    for image in directions:
        image.select_cc()
        skymodels.append(image.skymodel_cut)
    lsm = lsmtool.load(skymodels[0])
    for skymodel in skymodels[1:]:
        lsm2 = lsmtool.load(skymodel)
        lsm.concatenate(lsm2)
    lsm.write('ddcal/images/c'+str(c)+'/mos-sources-cut.txt', format='makesourcedb', clobber=True)

    os.system('cp img/*M*MFS-image.fits img/mos-MFS-image.fits img/mos-MFS-residual.fits ddcal/images/c'+str(c))
    mosaic_image = Image('ddcal/images/c'+str(c)+'/mos-MFS-image.fits', user_mask = user_mask)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = get_noise_img(mosaic_residual)
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise
