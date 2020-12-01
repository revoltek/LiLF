#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for extraction of target region after dd-serial

import sys, os, glob, re
import numpy as np
import pyregion
import lsmtool


#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-extract.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-extract.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_extract','parset_dir')
maxniter = parset.getint('LOFAR_extract','maxniter')
target_reg_file = parset.getstr('LOFAR_extract','extractRegion') # default 'target.reg'
userReg = parset.get('model','userReg')

############################
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('extract')
    os.makedirs('extract')
    lib_util.check_rm('extract/init')
    os.makedirs('extract/init')
    os.system('cp ddcal/c01/images/wideDD-c01.app.restored.fits extract/init/')
    os.system('cp ddcal/c01/solutions/interp.h5 extract/init/')
    lib_util.check_rm('mss-extract')
    if not os.path.exists('mss-extract'): os.system('cp -r mss-avg mss-extract')

# for now, region must be a single ds9 circle
target_reg = pyregion.open(target_reg_file)
if len(target_reg) > 1 or target_reg[0].name != 'circle':
    raise ImportError('Only single circle regions supported.')
target_reg = target_reg[0]
logger.info("Extract target region: {:.2f}deg, {:.2f}deg; radius: {:.3f}deg".format(*target_reg.coord_list))

MSs = lib_ms.AllMSs( glob.glob('mss-extract/*MS'), s )
ch_out = MSs.getChout(4e6)  # chout from dd-serial
fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
# read image, h5parm, make mask
wideDD_image = lib_img.Image('extract/init/wideDD-c01.app.restored.fits')
dde_h5parm = 'extract/init/interp.h5'
# make mask for subtraction
mask_ddcal = wideDD_image.imagename.replace('.fits', '_mask-ddcal.fits')  # this is used to find calibrators
wideDD_image.makeMask(threshpix=5, atrous_do=True, maskname=mask_ddcal, write_srl=True, write_ds9=True)


def clean(p, MSs, res='normal', size=[1, 1], empty=False, imagereg=None):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.resolution

    if res == 'normal':
        pixscale = float('%.1f' % (pixscale / 2.5))
    elif res == 'high':
        pixscale = float('%.1f' % (pixscale / 3.5))
    elif res == 'low':
        pass  # no change

    imsize = [int(size[0] * 1.5 / (pixscale / 3600.)), int(size[1] * 1.5 / (pixscale / 3600.))]  # add 50%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 256: imsize[0] = 256
    if imsize[1] < 256: imsize[1] = 256

    logger.debug('Image size: ' + str(imsize) + ' - Pixel scale: ' + str(pixscale))

    if res == 'normal':
        weight = 'briggs -0.3'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.6'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500
    else:
        logger.error('Wrong "res": %s.' % str(res))
        sys.exit()

    if empty:

        logger.info('Cleaning empty (' + str(p) + ')...')
        imagename = 'img/empty-' + str(p)
        lib_util.run_wsclean(s, 'wscleanE-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
                             data_column='SUBTRACTED_DATA',
                             size=imsize, scale=str(pixscale) + 'arcsec',
                             weight=weight, niter=0, no_update_model_required='', minuv_l=30, mgain=0,
                             baseline_averaging='')
    else:

        # clean 1
        logger.info('Cleaning (' + str(p) + ')...')
        imagename = 'img/extract-' + str(p)
        lib_util.run_wsclean(s, 'wscleanA-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
                             size=imsize, scale=str(pixscale) + 'arcsec',
                             weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l,
                             mgain=0.85,
                             baseline_averaging='', parallel_deconvolution=512, auto_threshold=5,
                             join_channels='', fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3)

        # make mask
        im = lib_img.Image(imagename + '-MFS-image.fits', userReg=userReg)
        try:
            im.makeMask(threshpix=10, rmsbox=(70, 5))
        except:
            logger.warning('Fail to create mask for %s.' % imagename + '-MFS-image.fits')
            return

        if imagereg is not None:
            lib_img.blank_image_reg(im.maskname, imagereg, inverse=True, blankval=0., )

        # clean 2
        # TODO: add deconvolution_channels when bug fixed
        logger.info('Cleaning w/ mask (' + str(p) + ')...')
        imagenameM = 'img/extractM-' + str(p)
        lib_util.run_wsclean(s, 'wscleanB-' + str(p) + '.log', MSs.getStrWsclean(), name=imagenameM, do_predict=True,
                             size=imsize, save_source_list='', scale=str(pixscale) + 'arcsec', reuse_psf=imagename,
                             reuse_dirty=imagename,
                             weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l,
                             mgain=0.85,
                             multiscale='', multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80',
                             baseline_averaging='', local_rms='', auto_threshold=0.75, auto_mask=1.5,
                             fits_mask=im.maskname,
                             join_channels='', fit_spectral_pol=3, channels_out=ch_out)  # , deconvolution_channels=3)

        os.system('cat logs/wscleanB-' + str(p) + '.log | grep "background noise"')


### DONE

with w.if_todo('predict_rest'):
    # DDF predict+corrupt in MODEL_DATA of everything BUT the calibrator
    indico = wideDD_image.root + '.DicoModel'
    outdico = indico + '-' + target_reg_file.split('.')[0] # use prefix of target reg
    inmask = sorted(glob.glob(wideDD_image.root + '.mask*.fits'))[-1]
    outmask = outdico + '.mask'
    lib_img.blank_image_reg(inmask, target_reg_file, outfile=outmask, inverse=False, blankval=0.)
    s.add('MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s' % (outmask, indico, outdico),
          log='MaskDicoModel.log', commandType='python', processors='max')
    s.run(check=True)

    ddf_parms = {
        'Data_MS': MSs.getStrDDF(),
        'Data_ColName': 'CORRECTED_DATA',
        'Data_Sort': 1,
        'Output_Mode': 'Predict',
        'Predict_InitDicoModel': outdico,
        'Predict_ColName': 'MODEL_DATA',
        'Deconv_Mode': 'HMP',
        'Image_NPix': 8775,
        'CF_wmax': 50000,
        'CF_Nw': 100,
        'Beam_CenterNorm': 1,
        'Beam_Smooth': 0,
        'Beam_Model': 'LOFAR',
        'Beam_LOFARBeamMode': 'A',
        'Beam_NBand': 1,
        'Beam_DtBeamMin': 5,
        'Output_Also': 'onNeds',
        'Image_Cell': 3.,
        'Freq_NDegridBand': ch_out,
        'Freq_NBand': ch_out,
        'Facets_DiamMax': 1.5,
        'Facets_DiamMin': 0.1,
        'Weight_ColName': 'WEIGHT_SPECTRUM',
        'Comp_BDAMode': 1,
        'DDESolutions_DDModeGrid': 'AP',
        'DDESolutions_DDModeDeGrid': 'AP',
        'RIME_ForwardMode': 'BDA-degrid',
        'DDESolutions_DDSols': dde_h5parm+':sol000/phase000+amplitude000'
    }
    logger.info('Predict corrupted rest-of-the-sky...')
    lib_util.run_DDF(s, 'ddfacet-pre.log', **ddf_parms, Cache_Reset=1)
### DONE

with w.if_todo('subtract_rest'):
    # Remove corrupted data from CORRECTED_DATA
    logger.info('Add columns...')
    MSs.run('addcol2ms.py -m $pathMS -c SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')
    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
            log='$nameMS_subtract.log', commandType='general')
    ## TTESTTESTTEST: empty image
    clean('but-target.fits', MSs, size=(fwhm,fwhm), res='normal', empty=True)


# Phase shift in the target location
with w.if_todo('phaseshift'):
    logger.info('Phase shift and avg...')
    lib_util.check_rm('mss-extract/*MS-small')
    MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-extract/$nameMS.MS-extract msin.datacolumn=SUBTRACTED_DATA \
            shift.phasecenter=['+str(target_reg.coord_list[0])+'deg,'+str(target_reg.coord_list[1])+'deg\]', \
            log='$nameMS_shiftavg.log', commandType='DPPP')
### DONE

MSs = lib_ms.AllMSs( glob.glob('mss-extract/*MS-extract'), s )

# initial imaging to get the model in the MODEL_DATA
with w.if_todo('image_init'):
    logger.info('Initial imaging...')
    clean('init', MSs, size=2*target_reg.coord_list[2])
### DONE

# Smoothing - ms:DATA -> ms:SMOOTHED_DATA
with w.if_todo('smooth'):
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -c 8 -n 6 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python') #, maxThreads=1)
### DONE

# get initial noise and set iterators for timeint solutions
image = lib_img.Image('img/extractM-init-MFS-image.fits', userReg=userReg)
rms_noise_pre = image.getNoise();
rms_noise_init = rms_noise_pre
mm_ratio_pre = image.getMaxMinRatio();
mm_ratio_init = mm_ratio_pre
doamp = False
# usually there are 3600/30=120 or 3600/15=240 timesteps, try to use multiple numbers
# integration time twice as long as DD-serial
iter_ph_solint = lib_util.Sol_iterator([8, 8, 4, 2, 1])
iter_amp_solint = lib_util.Sol_iterator([60, 30, 15, 10])
iter_amp2_solint = lib_util.Sol_iterator([120, 60, 30])
logger.info('RMS noise (init): %f' % (rms_noise_pre))
logger.info('MM ratio (init): %f' % (mm_ratio_pre))
rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    with w.if_todo('smooth-c%02i' % c):
        logger.info('BL-based smoothing on MODEL_DATA...')
        MSs.run('BLsmooth.py -c 8 -n 6 -r -i MODEL_DATA -o MODEL_DATA $pathMS', log='$nameMS_smoothM-c%02i.log' % c, commandType='python') #, maxThreads=1

    h5ph = 'extract/cal-ph-c%02i.h5' % c
    solint_ph = next(iter_ph_solint)
    if doamp:
        h5amp1 = 'extract/cal-amp1-c%02i.h5' % c
        solint_amp = next(iter_amp_solint)
        h5amp2 = 'extract/cal-amp2-c%02i.h5' % c
        solint_amp2 = next(iter_amp2_solint)

    logger.info('Phase calibration...')
    with w.if_todo('cal-ph-c%02i' % c):
        MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
            sol.mode=scalarcomplexgain sol.solint=' + str(solint_ph) + ' sol.nchan=1 sol.smoothnessconstraint=5e6 \
            sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]', \
                    log='$nameMS_solGph-%s.log' % c, commandType='DPPP')
        lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs.getListStr()],
                            [parset_dir + '/losoto-plot1.parset'],
                            plots_dir='extract/plots-%s' % c)
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
                                plots_dir='extract/plots-%s' % c)
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
                                plots_dir='extract/plots-%s' % c)
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
        clean('c%02i' % c, MSs, size=2*target_reg.coord_list[2], apply_beam = c==maxniter ) # size 2 times radius
    ### DONE

    # get noise, if larger than 95% of prev cycle: break
    extract_image = lib_img.Image('img/extractM-c%02i-MFS-image.fits' % c)
    # get noise, if larger than prev cycle: break
    rms_noise = extract_image.getNoise()
    mm_ratio = extract_image.getMaxMinRatio()
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