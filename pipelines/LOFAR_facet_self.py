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

def clean(p, MSs, size, res='normal', apply_beam=False, empty=False):
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

    imsize = [int(size[0]*1.5/(pixscale/3600.)), int(size[1]*1.5/(pixscale/3600.))] # add 50%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 256: imsize[0] = 256
    if imsize[1] < 256: imsize[1] = 256

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
    else:
        logger.error('Wrong "res": %s.' % str(res))
        sys.exit()

    if empty:

        logger.info('Cleaning empty ('+str(p)+')...')
        imagename = 'img/empty-'+str(p)
        lib_util.run_wsclean(s, 'wscleanE-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, data_column='SUBTRACTED_DATA', \
                             size=imsize, scale=str(pixscale)+'arcsec', \
                             weight=weight, niter=0, no_update_model_required='', minuv_l=30, mgain=0, \
                             baseline_averaging='')
    else:
        # clean 1
        logger.info('Cleaning ('+str(p)+')...')
        imagename = 'img/facet-'+str(p)
        lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
                weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
                baseline_averaging='', parallel_deconvolution=512, auto_threshold=5, \
                join_channels='', fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3)

        # make mask
        im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
        im.makeMask(threshisl = 5, rmsbox=(70,5))

        # clean 2
        logger.info('Cleaning w/ mask ('+str(p)+')...')
        imagename = 'img/facetM-'+str(p)
        if apply_beam:
            lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
                weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
                use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
                multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40,80', \
                parallel_deconvolution=512, local_rms='', auto_threshold=0.5, auto_mask=1, fits_mask=im.maskname, \
                join_channels='', fit_spectral_pol=3, channels_out=ch_out) #, deconvolution_channels=3)
        else:
            lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), do_predict=True, name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
                weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
                multiscale='', multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80', fits_mask=im.maskname,\
                baseline_averaging='', parallel_deconvolution=512, local_rms='', auto_threshold=0.75, auto_mask=1.5,  \
                join_channels='', fit_spectral_pol=3, channels_out=ch_out) #, deconvolution_channels=3)
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
    MSs.run('BLsmooth.py -c 8 -n 6 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python', maxThreads=1)
### DONE


# get initial noise and set iterators for timeint solutions
image = lib_img.Image('img/facetM-init-MFS-image.fits', userReg=userReg)
rms_noise_pre = image.getNoise();
rms_noise_init = rms_noise_pre
mm_ratio_pre = image.getMaxMinRatio();
mm_ratio_init = mm_ratio_pre
doamp = False
# usually there are 3600/30=120 or 3600/15=240 timesteps, try to use multiple numbers
iter_ph_solint = lib_util.Sol_iterator([4, 1])
iter_amp_solint = lib_util.Sol_iterator([120, 60, 30, 10])
iter_amp2_solint = lib_util.Sol_iterator([120, 60, 30])
logger.info('RMS noise (init): %f' % (rms_noise_pre))
logger.info('MM ratio (init): %f' % (mm_ratio_pre))
rms_noise_pre = np.inf

for c in range(20):
    logger.info('Starting cycle: %i' % c)

    with w.if_todo('smooth-c%02i' % c):
        logger.info('BL-based smoothing on MODEL_DATA...')
        MSs.run('BLsmooth.py -c 8 -n 6 -r -i MODEL_DATA -o MODEL_DATA $pathMS', log='$nameMS_smoothM-c%02i.log' % c, commandType='python', maxThreads=1)

    h5ph = 'facet/cal-ph-c%02i.h5' % c
    solint_ph = next(iter_ph_solint)
    if doamp:
        h5amp1 = 'facet/cal-amp1-c%02i.h5' % c
        solint_amp = next(iter_amp_solint)
        h5amp2 = 'facet/cal-amp2-c%02i.h5' % c
        solint_amp2 = next(iter_amp2_solint)

    logger.info('Phase calibration...')
    antconstr = '[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,' \
                'CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,' \
                'CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,' \
                'CS302LBA,CS401LBA,CS501LBA]]'
    with w.if_todo('cal-ph-c%02i' % c):
        MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
            sol.mode=scalarcomplexgain sol.solint=' + str(solint_ph) + ' sol.nchan=1 sol.smoothnessconstraint=5e6 \
            sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]', \
                    log='$nameMS_solGph-%s.log' % c, commandType='DPPP')
        lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs.getListStr()],
                            [parset_dir + '/losoto-plot1.parset'],
                            plots_dir='facet/plots-%s' % c)
        os.system('mv cal-ph.h5 %s' % h5ph)

    with w.if_todo('cor-ph-c%02i' % c):
        # correct ph - ms:DATA -> ms:CORRECTED_DATA
        logger.info('Correct ph...')
        MSs.run('DPPP ' + parset_dir + '/DPPP-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                     cor.parmdb=' + h5ph + ' cor.correction=phase000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DPPP')
    if doamp:
        with w.if_todo('cal-amp1-c%02i' % c):
            logger.info('Gain amp calibration 1 (solint: %i)...' % solint_amp)
            # Calibration - ms:CORRECTED_DATA
            # possible to put nchan=6 if less channels are needed in the h5parm (e.g. for IDG)
            MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp1.h5 \
                sol.mode=diagonal sol.solint=' + str(solint_amp) + ' sol.nchan=1 sol.uvmmin=100 sol.smoothnessconstraint=4e6 sol.minvisratio=0.5\
                sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]]', \
                        log='$nameMS_solGamp1-c%02i.log' % c, commandType='DPPP')

            losoto_parsets = [parset_dir + '/losoto-clip.parset', parset_dir + '/losoto-norm.parset',
                                  parset_dir + '/losoto-plot2.parset']
            lib_util.run_losoto(s, 'amp1', [ms + '/cal-amp1.h5' for ms in MSs.getListStr()], losoto_parsets,
                                plots_dir='facet/plots-%s' % c)
            os.system('mv cal-amp1.h5 %s' % h5amp1)

        with w.if_todo('cor-amp1-c%02i' % c):
            logger.info('Correct amp 1...')
            # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs.run('DPPP ' + parset_dir + '/DPPP-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp1 + ' cor.correction=amplitude000',
                        log='$nameMS_correct-c%02i.log' % c, commandType='DPPP')

        with w.if_todo('cal-amp2-c%02i' % c):
            logger.info('Gain amp calibration 2 (solint: %i)...' % solint_amp2)
            # Calibration - ms:SMOOTHED_DATA
            MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp2.h5 \
                sol.mode=diagonal sol.solint=' + str(
                solint_amp2) + ' sol.nchan=6 sol.uvmmin=100 sol.smoothnessconstraint=10e6 sol.minvisratio=0.5', \
                        log='$nameMS_solGamp2-c%02i.log' % s, commandType='DPPP')

            losoto_parsets = [parset_dir + '/losoto-clip2.parset', parset_dir + '/losoto-norm.parset',
                              parset_dir + '/losoto-plot3.parset']
            lib_util.run_losoto(s, 'amp2', [ms + '/cal-amp2.h5' for ms in MSs.getListStr()], losoto_parsets,
                                plots_dir='facet/plots-%s' % c)
            os.system('mv cal-amp2.h5 %s' % h5amp2)

        with w.if_todo('cor-amp2-c%02i' % c):
            logger.info('Correct amp 2...')
            # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs.run('DPPP ' + parset_dir + '/DPPP-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp2 + ' cor.correction=amplitude000',
                        log='$nameMS_correct-c%02i.log' % c, commandType='DPPP')
    ### DONE

    with w.if_todo('image-c%02i' % c):
        logger.info('Imaging...')
        clean('c%02i' % c, MSs, size=size, apply_beam = c==maxniter )


    ### DONE

    # get noise, if larger than 95% of prev cycle: break
    facet_image = lib_img.Image('img/facetM-c%02i-MFS-image.fits' % c)
    # get noise, if larger than prev cycle: break
    rms_noise = facet_image.getNoise()
    mm_ratio = facet_image.getMaxMinRatio()
    logger.info('RMS noise (c:%02i): %f' % (c, rms_noise))
    logger.info('MM ratio (c:%02i): %f' % (c, mm_ratio))
    if rms_noise > 0.99 * rms_noise_pre and mm_ratio < 1.01 * mm_ratio_pre and c >4:
        if (mm_ratio < 10 and c >= 2) or \
                (mm_ratio < 20 and c >= 3) or \
                (mm_ratio < 30 and c >= 4) or \
                (c >= 5): break

    if c >= 4 and mm_ratio >= 30:
        logger.info('Start amplitude calibration in next cycle...')
        doamp = True

    rms_noise_pre = rms_noise
    mm_ratio_pre = mm_ratio