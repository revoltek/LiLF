#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for single facet self calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool


#TODO:
size=[5.5,3.0]

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-dd.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)#, maxThreads = 4)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_facet_self','parset_dir')
maxniter = parset.getint('LOFAR_facet_self','maxniter')
userReg = parset.get('model','userReg')

####################################################
MSs_self = lib_ms.AllMSs( glob.glob('mss-facet/TC*[0-9].MS'), s )

############################
logger.info('Cleaning...')

def clean(p, MSs, size, res='normal', apply_beam=False):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.getListObj()[0].getResolution() 
    if res == 'normal':
        pixscale = float('%.1f'%(pixscale/2.5))
    elif res == 'high':
        pixscale = float('%.1f'%(pixscale/3.5))
    elif res == 'low':
        pass # no change

    imsize = [0,0]
    imsize[0] = int(size[0]*1.05/(pixscale/3600.)) # add 5%
    imsize[1] = int(size[1]*1.05/(pixscale/3600.)) # add 5%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 64: imsize[0] == 64
    if imsize[1] < 64: imsize[1] == 64

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    if res == 'normal':
        weight = 'briggs -0.1'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.6'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/facet-'+str(p)
    lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=512, auto_threshold=5, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 3)

    # clean 2
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/facetM-'+str(p)
    if apply_beam:

        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40,80', \
            parallel_deconvolution=512, local_rms='', auto_threshold=1., auto_mask=2., fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    else:

        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), do_predict=True, name=imagename, size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40,80', 
            baseline_averaging=5, parallel_deconvolution=512, local_rms='', auto_threshold=1., auto_mask=2., fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    os.system('cat logs/wscleanB-'+str(p)+'.log | grep "background noise"')

MSs = lib_ms.AllMSs( glob.glob('mss-facet/*MS'), s )

# initial imaging to get the model in the MODEL_DATA
logger.info('Initial imaging...')
clean('init', MSs, s ), size=size )

# Smoothing - ms:DATA -> ms:SMOOTHED_DATA
logger.info('BL-based smoothing...')
MSs.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')

rms_noise_pre = np.inf
for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Gain calibration...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS \
            sol.h5parm=$pathMS/cal-g-c'+str(c)+'.h5', \
            log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    lib_util.run_losoto(s, 'g-c'+str(c), [ms+'/cal-g-c'+str(c)+'.h5' for ms in MSs.getListStr()], \
                [parset_dir+'/losoto-amp.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset'])

    # correct G - ms:DATA -> ms:CORRECTED_DATA
    logger.info('Patch  correct...')
    MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS \
                cor.parmdb=cal-g-c'+str(c)+'.h5 cor.correction=phase000', \
                log='$nameMS_correct-c'+str(c)+'.log', commandType='DPPP') 
    if c>0:
        MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn = CORRECTED_DATA \
                cor.parmdb=cal-g-c'+str(c)+'.h5 cor.correction=amplitude000', \
                log='$nameMS_correct-c'+str(c)+'.log', commandType='DPPP') 

    logger.info('Imaging...')
    clean('%02i' % c, MSs, size=size, apply_beam = c==maxniter )

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = mosaic_image.getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > rms_noise_pre: break
    rms_noise_pre = rms_noise
