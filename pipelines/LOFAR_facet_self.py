#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for single facet self calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

#TODO, to extract from regionfile:
lastcycle = 0 # get from somewhere
target_reg = 'target.reg'

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-facet_self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)#, maxThreads = 4)
w = lib_util.Walker('pipeline-facet-self.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_facet_self','parset_dir')
maxniter = parset.getint('LOFAR_facet_self','maxniter')
userReg = parset.get('model','userReg')
mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % lastcycle)

############################
with w.if_todo('cleaning'):

    logger.info('Cleaning...')
    lib_util.check_rm('plot*')
    lib_util.check_rm('cal*h5')
    lib_util.check_rm('img')
    lib_util.check_rm('facet')
    os.makedirs('img')
    os.makedirs('facet')
    lib_util.check_rm('mss-facet')
    if not os.path.exists('mss-facet'): os.system('cp -r mss-dd mss-facet')
    

### DONE

MSs = lib_ms.AllMSs( glob.glob('mss-facet/*MS'), s )

def clean(p, MSs, size, res='normal', apply_beam=False):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.resolution
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
        maxuv_l = None

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/facet-'+str(p)
    lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            baseline_averaging='', parallel_deconvolution=512, auto_threshold=5, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshpix = 5)

    # clean 2
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/facetM-'+str(p)
    if apply_beam:

        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40', \
            parallel_deconvolution=512, local_rms='', auto_threshold=0.5, auto_mask=1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    else:

        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), do_predict=True, name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40', \
            baseline_averaging='', parallel_deconvolution=512, local_rms='', auto_threshold=0.5, auto_mask=1., fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

        lib_util.run_wsclean(s, 'wscleanBlow-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename+'-low', size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, taper_gaussian='30asec', mgain=0.85, \
            multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40', \
            baseline_averaging='', parallel_deconvolution=512, local_rms='', auto_threshold=0.5, auto_mask=1., fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    os.system('cat logs/wscleanB-'+str(p)+'.log | grep "background noise"')

# Load facet mask and set target region to 0
mask_voro = 'ddcal/masks/facets%02i.fits' % lastcycle
os.system('cp %s facet/facets.fits' % mask_voro)
lib_img.blank_image_reg('facet/facets.fits', target_reg, blankval=0)

# mosaic the skymodel, Isl_patch_000 will be the target of interest
lsm = lsmtool.load(mosaic_image.skymodel_cut)
lsm.group('facet', facet='facet/facets.fits', root='Isl_patch')
lsm.setPatchPositions(method='mid') # center of the facets
directions = set(lsm.getColValues('patch'))
coord = [c.deg for c in lsm.getPatchPositions('Isl_patch_0')['Isl_patch_0']]
logger.info("Facet centre: "+str(coord))

# calculate region size (TODO: maybe better using regionfile?)
ramin = np.min(lsm.getColValues('RA')[lsm.getColValues('Patch')=='Isl_patch_0'])
ramax = np.max(lsm.getColValues('RA')[lsm.getColValues('Patch')=='Isl_patch_0'])
decmin = np.min(lsm.getColValues('Dec')[lsm.getColValues('Patch')=='Isl_patch_0'])
decmax = np.max(lsm.getColValues('Dec')[lsm.getColValues('Patch')=='Isl_patch_0'])
#size = [abs(ramax-ramin),abs(decmax-decmin)]
size = [1,1]

# write skymodel
lsm.write('facet/skymodel_init.txt', format='makesourcedb', clobber=True)
lib_util.check_rm('facet/skymodel_init.skydb')
s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % ('facet/skymodel_init.txt', 'facet/skymodel_init.skydb'), \
	log='makesourcedb.log', commandType='general')
s.run(check=True)

del lsm

with w.if_todo('prepare_col'):
    logger.info('Subtraction...')
    # Copy DATA -> SUBTRACTED_DATA
    logger.info('Add columns...')
    MSs.run('addcol2ms.py -m $pathMS -c SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')
    logger.info('Set SUBTRACTED_DATA = DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql.log', commandType='general')


### DONE

# subtract all sources outside the region of interest
with w.if_todo('subtract'):
    for d in directions:
        if d == 'Isl_patch_0':
            # this is the target of interest
            continue
    
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=facet/skymodel_init.skydb pre.sources='+d, \
                log='$nameMS_pre-'+d+'.log', commandType='DPPP')
    
        # corrupt G - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+d+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt1.parset msin=$pathMS \
                cor.parmdb=ddcal/solutions/cal-g-c'+str(lastcycle)+'.h5 cor.correction=phase000 cor.direction=['+d+']', \
                log='$nameMS_corrupt1-'+d+'.log', commandType='DPPP')
        #MSs.run('DPPP '+parset_dir+'/DPPP-corrupt1.parset msin=$pathMS \
        #        cor.parmdb=ddcal/solutions/cal-g-c'+str(lastcycle)+'.h5 cor.correction=amplitude000 cor.direction=['+d+']', \
        #        log='$nameMS_corrupt2-'+d+'.log', commandType='DPPP')
    
        logger.info('Patch '+d+': subtract...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql-'+d+'.log', commandType='general')
 

### DONE

# Phase shift in the target location
with w.if_todo('phaseshift'):

    logger.info('Phase shift and avg...')
    lib_util.check_rm('mss-facet/*MS-small')
    MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-facet/$nameMS.MS-small msin.datacolumn=SUBTRACTED_DATA \
            shift.phasecenter=['+str(coord[0])+'deg,'+str(coord[1])+'deg\]', \
            log='$nameMS_avgshift.log', commandType='DPPP')

### DONE

MSs = lib_ms.AllMSs( glob.glob('mss-facet/*MS-small'), s )

# initial imaging to get the model in the MODEL_DATA
with w.if_todo('image_init'):
    logger.info('Initial imaging...')
    clean('init', MSs, size=size )


### DONE

# Smoothing - ms:DATA -> ms:SMOOTHED_DATA
with w.if_todo('smooth'):
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')


### DONE

rms_noise_pre = np.inf
for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    with w.if_todo('solve-c%02i' % c):
        # Calibration - ms:SMOOTHED_DATA
        logger.info('Gain calibration...')
        solint = max(8-c,1)
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.solint='+str(solint)+'\
                sol.h5parm=$pathMS/cal-g-c'+str(c)+'.h5', \
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
    
        # Plot solutions
        lib_util.run_losoto(s, 'g-c'+str(c), [ms+'/cal-g-c'+str(c)+'.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-amp.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset'])
        os.system('mv plots-g-c%i facet' % c)
    

    ### DONE

    with w.if_todo('cor-c%02i' % c):
        # correct G - ms:DATA -> ms:CORRECTED_DATA
        logger.info('Ph correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS \
                    cor.parmdb=cal-g-c'+str(c)+'.h5 cor.correction=phase000', \
                    log='$nameMS_correctPH-c'+str(c)+'.log', commandType='DPPP') 
        if c>1:
            logger.info('Amp correct...')
            MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                    cor.parmdb=cal-g-c'+str(c)+'.h5 cor.correction=amplitude000', \
                    log='$nameMS_correctAMP-c'+str(c)+'.log', commandType='DPPP') 
    

    ### DONE

    with w.if_todo('image-c%02i' % c):
        logger.info('Imaging...')
        clean('c%02i' % c, MSs, size=size, apply_beam = c==maxniter )


    ### DONE

    # get noise, if larger than 95% of prev cycle: break
    facet_image = lib_img.Image('img/facetM-c%02i-MFS-image.fits' % c)
    rms_noise = facet_image.getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > rms_noise_pre and c>=5: break
    rms_noise_pre = rms_noise
