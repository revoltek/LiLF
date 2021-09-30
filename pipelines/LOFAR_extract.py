#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for extraction of target region after LOFAR_dd-serial.
# Provide a extractRegion in lilf.config. This pipeline will subtract
# sources outside of the region and perform subsequent self-calibration.
# A userReg may be specified as clean mask.
# phSolMode can be used to solve either using phases or phaseandtec.

import sys, os, glob, re
import numpy as np
import pyregion
from astropy.io import fits
from losoto.h5parm import h5parm

def get_ddf_parms_from_header(img):
    """
    Parse the HISTORY header of a DDFacet image and return a dict containing the options used to create the image.
    Will replace '-' by '_'.
    Parameters
    ----------
    img: str, filename of image

    Returns
    -------
    params_dict: dict,
    """
    params_dict = dict()
    hdr = fits.open(img)[0].header['HISTORY']
    for line in hdr:
        if line.count('=') == 1 and line.count('-') > 0:
            _key, _value = line.replace(' ', '').split('=')
            _key = _key.replace('-', '_')
            params_dict[_key] = _value
        else:
            continue
    return params_dict


#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-extract.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-extract.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_extract','parset_dir')
maxniter = parset.getint('LOFAR_extract','maxniter')
target_reg_file = parset.get('LOFAR_extract','extractRegion')  # default 'target.reg'
phSolMode = parset.get('LOFAR_extract','phSolMode')  # default: tecandphase
if phSolMode not in ['tecandphase', 'phase']:
    logger.error('phSolMode {} not supported. Choose tecandphase, phase.')
    sys.exit()
userReg = parset.get('model','userReg')

############################
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('extract')
    lib_util.check_rm('img')
    os.makedirs('img')
    os.makedirs('extract/init')
    os.system('cp ddcal/c01/images/wideDD-c01.app.restored.fits extract/init/') # copy ddcal image
    os.system('cp ddcal/c01/images/wideDD-c01.DicoModel extract/init/') # copy dico model
    os.system('cp ddcal/c01/solutions/interp.h5 extract/init/') # copy fnal dde sol
    lib_util.check_rm('mss-extract')
    if not os.path.exists('mss-extract'):
        logger.info('Copy MS...')
        os.system('cp -r mss-avg mss-extract')

# region must be a list of ds9 circles and polygons (other shapes can be implemented in lib_util.Rgion_helper()
target_reg = lib_util.Region_helper(target_reg_file)
center = target_reg.get_center() # center of the extract region
MSs = lib_ms.AllMSs( glob.glob('mss-extract/*MS'), s )

ch_out = MSs.getChout(4e6)  # chout from dd-serial
fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
phase_center = MSs.getListObj()[0].getPhaseCentre()
# read image, h5parm, make mask
wideDD_image = lib_img.Image('extract/init/wideDD-c01.app.restored.fits')
dde_h5parm = 'extract/init/interp.h5'
# make mask for subtraction
mask_ddcal = wideDD_image.imagename.replace('.fits', '_mask-ddcal.fits')  # this is used to find calibrators
wideDD_image.makeMask(threshpix=5, atrous_do=True, maskname=mask_ddcal, write_srl=True, write_ds9=True)


def clean(p, MSs, res='normal', size=[1, 1], empty=False, imagereg=None, apply_beam=False):
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
        arg_dict = dict()
        if apply_beam:
            arg_dict['use_idg'] = ''
            arg_dict['grid_with_beam'] = ''
        else:
            arg_dict['baseline_averaging'] = ''
        # clean 1
        logger.info('Cleaning (' + str(p) + ')...')
        imagename = 'img/extract-' + str(p)

        lib_util.run_wsclean(s, 'wscleanA-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
                             size=imsize, scale=str(pixscale) + 'arcsec',
                             weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l,
                             mgain=0.85, parallel_deconvolution=512, auto_threshold=5, join_channels='',
                             fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3, **arg_dict)

        # make mask
        im = lib_img.Image(imagename + '-MFS-image.fits', userReg=userReg)
        try:
            im.makeMask(threshpix=10, rmsbox=(70, 5))
        except:
            logger.warning('Fail to create mask for %s.' % imagename + '-MFS-image.fits')
            return

        if imagereg is not None:
            lib_img.blank_image_reg(im.maskname, imagereg, inverse=True, blankval=0.)

        # clean 2
        # TODO: add deconvolution_channels when bug fixed
        logger.info('Cleaning w/ mask (' + str(p) + ')...')
        imagenameM = 'img/extractM-' + str(p)
        if apply_beam: arg_dict['reuse_primary_beam'] = imagename
        lib_util.run_wsclean(s, 'wscleanB-' + str(p) + '.log', MSs.getStrWsclean(), name=imagenameM, do_predict=True,
                             size=imsize, scale=str(pixscale) + 'arcsec',
                             weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l,
                             reuse_psf=imagename, reuse_dirty=imagename,  mgain=0.85, multiscale='',
                             parallel_deconvolution=512, auto_threshold=0.7, auto_mask=1.5,
                             fits_mask=im.maskname, join_channels='', fit_spectral_pol=3, channels_out=ch_out, **arg_dict)  # , deconvolution_channels=3)
        os.system('cat logs/wscleanB-' + str(p) + '.log | grep "background noise"')


with w.if_todo('predict_rest'):
    # DDF predict+corrupt in MODEL_DATA of everything BUT the calibrator
    indico = wideDD_image.root + '.DicoModel'
    outdico = indico + '-' + target_reg_file.split('.')[0] # use prefix of target reg
    inmask = sorted(glob.glob(wideDD_image.root + '*_mask-ddcal.fits'))[-1]
    outmask = outdico + '.mask'
    lib_img.blank_image_reg(inmask, target_reg_file, outfile=outmask, inverse=False, blankval=0.)
    s.add('MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s' % (outmask, indico, outdico),
          log='MaskDicoModel.log', commandType='python', processors='max')
    s.run(check=True)

    # get DDF parameters used to create the image/model
    ddf_parms = get_ddf_parms_from_header(wideDD_image.imagename)
    # change for PREDICT
    ddf_parms['Data_MS'] = MSs.getStrDDF()
    ddf_parms['Data_ColName'] = 'CORRECTED_DATA'
    ddf_parms['Output_Mode'] = 'Predict'
    ddf_parms['Predict_InitDicoModel'] = outdico
    ddf_parms['Beam_Smooth'] = 1
    ddf_parms['Cache_Reset'] = 1
    ddf_parms['DDESolutions_DDSols'] = dde_h5parm + ':sol000/phase000+amplitude000'
    if 'Misc_ParsetVersion' in ddf_parms.keys(): del ddf_parms['Misc_ParsetVersion']

    logger.info('Predict corrupted rest-of-the-sky...')
    lib_util.run_DDF(s, 'ddfacet-pre.log', **ddf_parms)
    ### DONE

with w.if_todo('subtract_rest'):
    # Remove corrupted data from CORRECTED_DATA
    logger.info('Add columns...')
    MSs.run('addcol2ms.py -m $pathMS -c SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')
    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
            log='$nameMS_subtract.log', commandType='general')
    ### DONE

## TTESTTESTTEST: empty image
if not os.path.exists('img/empty-but-target-image.fits'):
    clean('but-target', MSs, size=(fwhm,fwhm), res='normal', empty=True)
    ### DONE

# Phase shift in the target location
with w.if_todo('phaseshift'):
    logger.info('Phase shift and avg...')
    MSs.run('DP3 '+parset_dir+'/DP3-shiftavg.parset msin=$pathMS msout=mss-extract/$nameMS.MS-extract msin.datacolumn=SUBTRACTED_DATA \
            shift.phasecenter=['+str(center[0])+'deg,'+str(center[1])+'deg\] avg.freqstep=8 avg.timestep=4', \
            log='$nameMS_shiftavg.log', commandType='DP3')
    ### DONE

MSs = lib_ms.AllMSs( glob.glob('mss-extract/*MS-extract'), s )

with w.if_todo('beamcorr'):
    logger.info('Correcting beam...') # TODO is this correct?
    # Convince DP3 that DATA is corrected for the beam in the phase centre
    MSs.run('DP3 ' + parset_dir + '/DP3-beam.parset msin=$pathMS setbeam.direction=\[' + str(phase_center[0]) + 'deg,'
            + str(phase_center[1]) + 'deg\] corrbeam.direction=\[' + str(center[0]) + 'deg,' + str(
        str(center[1])) + 'deg\]', log='$nameMS_beam-.log', commandType='DP3')
    ### DONE

# apply init - closest DDE sol
# TODO: this assumes phase000 and optionally, amplitude000
with w.if_todo('apply_init'):
    h5init = h5parm(dde_h5parm)
    solset_dde = h5init.getSolset('sol000')
    # get closest dir to target reg center
    dirs = np.array([solset_dde.getSou()[k] for k in solset_dde.getSoltab('phase000').dir])
    dir_dist = lib_util.distanceOnSphere(dirs[:,0], dirs[:,1],*np.deg2rad(center), rad=True)
    closest = solset_dde.getSoltab('phase000').dir[np.argmin(dir_dist)]
    logger.info('Init apply: correct closest DDE solutions ({})'.format(closest))
    logger.info('Correct init ph...')
    MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA '
                                   'msout.datacolumn=CORRECTED_DATA cor.parmdb=' + dde_h5parm + ' cor.correction=phase000 cor.direction='+closest,
            log='$nameMS_init-correct.log', commandType='DP3')
    if 'amplitude000' in solset_dde.getSoltabNames():
        logger.info('Correct init amp...')
        MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                     cor.parmdb=' + dde_h5parm + ' cor.correction=amplitude000 cor.direction=' + closest,
                log='$nameMS_init-correct.log', commandType='DP3')
    h5init.close()
    ### DONE

# initial imaging to get the model in the MODEL_DATA (could also be done using the Dico DDFacet model
with w.if_todo('image_init'):
    logger.info('Initial imaging...')
    clean('init', MSs, size=(1.1*target_reg.get_width(),1.1*target_reg.get_height()))
    ### DONE

# Smoothing - ms:DATA -> ms:SMOOTHED_DATA
with w.if_todo('smooth'):
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -c 1 -n 8 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python', maxThreads=1)
    ### DONE

# get initial noise and set iterators for timeint solutions
image = lib_img.Image('img/extractM-init-MFS-image.fits', userReg=userReg)
rms_noise_pre, mm_ratio_pre = image.getNoise(), image.getMaxMinRatio()
rms_noise_init, mm_ratio_init = rms_noise_pre, mm_ratio_pre
doamp = False

# Per default we have 32s exposure
iter_ph_solint = lib_util.Sol_iterator([8, 8, 4, 2, 1])
iter_amp_solint = lib_util.Sol_iterator([60, 60, 30, 15])
iter_amp2_solint = lib_util.Sol_iterator([120, 60, 60, 30])
logger.info('RMS noise (init): %f' % (rms_noise_pre))
logger.info('MM ratio (init): %f' % (mm_ratio_pre))
rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    h5ph = 'extract/cal-ph-c%02i.h5' % c
    solint_ph = next(iter_ph_solint)
    if doamp:
        h5amp1 = 'extract/cal-amp1-c%02i.h5' % c
        solint_amp = next(iter_amp_solint)
        h5amp2 = 'extract/cal-amp2-c%02i.h5' % c
        solint_amp2 = next(iter_amp2_solint)

    if phSolMode == 'phase':
        logger.info('Phase calibration...')
        with w.if_todo('cal-ph-c%02i' % c):
            MSs.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
                     sol.mode=scalarphase sol.solint=' + str(solint_ph) + ' sol.smoothnessconstraint=5e6 sol.smoothnessreffrequency=54e6 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]',
                     log='$nameMS_solGph-c%02i.log' % c, commandType='DP3')
            lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs.getListStr()],
                                [parset_dir + '/losoto-plot1.parset'],
                                plots_dir='extract/plots-%s' % c)
            os.system('mv cal-ph.h5 %s' % h5ph)

        with w.if_todo('cor-ph-c%02i' % c):
            # correct ph - ms:DATA -> ms:CORRECTED_DATA
            logger.info('Correct ph...')
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=phase000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DP3')
    else:
        logger.info('Tecandphase calibration...')
        with w.if_todo('cal-tecandph-c%02i' % c):
            MSs.run('DP3 ' + parset_dir + '/DP3-soltecandphase.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
                     sol.mode=tecandphase sol.solint=' + str(solint_ph) + ' sol.nchan=1 sol.smoothnessconstraint=5e6 sol.smoothnessreffrequency=54e6 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]',
                     log='$nameMS_soltecandphase-c%02i.log' % c, commandType='DP3')
            lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs.getListStr()],
                                [parset_dir + '/losoto-plottecandphase.parset'],
                                plots_dir='extract/plots-%s' % c)
            os.system('mv cal-ph.h5 %s' % h5ph)
        logger.info('tecandphase calibration...')

        with w.if_todo('cor-tecandph-c%02i' % c):
            # correct ph - ms:DATA -> ms:CORRECTED_DATA
            logger.info('Correct tec...')
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=tec000',
                        log='$nameMS_correct-c%02i.log' % c, commandType='DP3')
            logger.info('Correct ph...')
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=phase000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DP3')

    if doamp:
        with w.if_todo('cal-amp1-c%02i' % c):
            logger.info('Gain amp calibration 1 (solint: %i)...' % solint_amp)
            # Calibration - ms:CORRECTED_DATA
            # possible to put nchan=6 if less channels are needed in the h5parm (e.g. for IDG)
            MSs.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp1.h5 \
                sol.mode=diagonal sol.solint=' + str(solint_amp) + ' sol.nchan=1 sol.smoothnessconstraint=4e6 sol.minvisratio=0.5 \
                sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]]', \
                        log='$nameMS_solGamp1-c%02i.log' % c, commandType='DP3')

            losoto_parsets = [parset_dir + '/losoto-clip.parset', parset_dir + '/losoto-norm.parset',
                                  parset_dir + '/losoto-plot2.parset']
            lib_util.run_losoto(s, 'amp1', [ms + '/cal-amp1.h5' for ms in MSs.getListStr()], losoto_parsets,
                                plots_dir='extract/plots-%s' % c)
            os.system('mv cal-amp1.h5 %s' % h5amp1)

        with w.if_todo('cor-amp1-c%02i' % c):
            logger.info('Correct amp 1...')
            # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp1 + ' cor.correction=amplitude000',
                        log='$nameMS_correct-c%02i.log' % c, commandType='DP3')

        with w.if_todo('cal-amp2-c%02i' % c):
            logger.info('Gain amp calibration 2 (solint: %i)...' % solint_amp2)
            # Calibration - ms:SMOOTHED_DATA
            MSs.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp2.h5 \
                sol.mode=diagonal sol.solint=' + str(
                solint_amp2) + ' sol.nchan=6  sol.smoothnessconstraint=10e6 sol.minvisratio=0.5', \
                        log='$nameMS_solGamp2-c%02i.log' % c, commandType='DP3')

            losoto_parsets = [parset_dir + '/losoto-clip2.parset', parset_dir + '/losoto-norm.parset',
                              parset_dir + '/losoto-plot3.parset']
            lib_util.run_losoto(s, 'amp2', [ms + '/cal-amp2.h5' for ms in MSs.getListStr()], losoto_parsets,
                                plots_dir='extract/plots-%s' % c)
            os.system('mv cal-amp2.h5 %s' % h5amp2)

        with w.if_todo('cor-amp2-c%02i' % c):
            logger.info('Correct amp 2...')
            # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp2 + ' cor.correction=amplitude000',
                        log='$nameMS_correct-c%02i.log' % c, commandType='DP3')
    ### DONE

    with w.if_todo('image-c%02i' % c):
        logger.info('Imaging...')
        # TODO: Apply beam for last iteration
        clean('c%02i' % c, MSs, size=(1.1*target_reg.get_width(),1.1*target_reg.get_height())) # size 2 times radius  , apply_beam = c==maxniter
    ### DONE

    # get noise, if larger than 95% of prev cycle: break
    extract_image = lib_img.Image('img/extractM-c%02i-MFS-image.fits' % c)
    # get noise, if larger than prev cycle: break
    rms_noise = extract_image.getNoise()
    mm_ratio = extract_image.getMaxMinRatio()
    logger.info('RMS noise (c:%02i): %f' % (c, rms_noise))
    logger.info('MM ratio (c:%02i): %f' % (c, mm_ratio))
    if rms_noise < rms_noise_pre:
        best_iter = c
    else: best_iter = c - 1

    if rms_noise > 0.99 * rms_noise_pre and mm_ratio < 1.01 * mm_ratio_pre and c >4:
        if (mm_ratio < 10 and c >= 2) or \
                (mm_ratio < 20 and c >= 3) or \
                (mm_ratio < 30 and c >= 4) or \
                (c >= 5):
            break
    # save best iteration

    if c >= 4 and mm_ratio >= 30:
        logger.info('Start amplitude calibration in next cycle...')
        doamp = True

    rms_noise_pre = rms_noise
    mm_ratio_pre = mm_ratio
    ### DONE

# Finally:
with w.if_todo('apply_final'):
    if best_iter != c: # If last iteration was NOT the best iteration, apply best iteration.
        logger.info('Best iteration: second to last cycle ({})'.format(best_iter))
        h5ph = 'extract/cal-ph-c%02i.h5' % best_iter
        h5amp1 = 'extract/cal-amp1-c%02i.h5' % best_iter
        h5amp2 = 'extract/cal-amp2-c%02i.h5' % best_iter
        # correct ph - ms:DATA -> ms:CORRECTED_DATA
        logger.info('Correct ph...')
        MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                     cor.parmdb=' + h5ph + ' cor.correction=phase000',
                log='$nameMS_correct-final.log', commandType='DP3')

        if phSolMode == 'tecandphase':
            logger.info('Correct tec...')
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=tec000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DP3')

        if os.path.exists(h5amp1):
            logger.info('Correct amp 1...')
            # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp1 + ' cor.correction=amplitude000',
                    log='$nameMS_correct-final.log', commandType='DP3')

        if os.path.exists(h5amp2):
            logger.info('Correct amp 2...')
            # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp2 + ' cor.correction=amplitude000',
                    log='$nameMS_correct-final.log', commandType='DP3')
    else:
        logger.info('Best ieration: last cycle ({})'.format(best_iter))

with w.if_todo('imaging_final'):
    logger.info('Final imaging w/ beam correction...')
    clean('final', MSs, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=True)  # size 2 times radius  , apply_beam = c==maxniter
    logger.info('Done.')
