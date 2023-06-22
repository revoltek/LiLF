#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Peel a bright source outside of the PB from the data
# Two possible ways to run this pipeline:
# 1. You provide either a ds9 circle/polygon peelReg (default 'peel.reg') and either a sourcedb or a fits_model
#    in the [model] section. The pipeline will solve against this initial model and then continue w/ self-cal.
# The MS files need to be in "./mss/"

import sys, os, glob, re
from shutil import copy2, copytree, move
import numpy as np
import astropy.io.fits as fits
import casacore.tables as pt
import lsmtool as lsm

from LiLF import lib_ms, lib_img, lib_util, lib_log
#############################################################################
logger_obj = lib_log.Logger('pipeline-peel.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-peel.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_peel','parset_dir')
cal_dir = parset.get('LOFAR_peel','cal_dir')
data_dir = parset.get('LOFAR_peel','data_dir')
peelReg = parset.get('LOFAR_peel','peelReg')  # default peel.reg
predictReg = parset.get('LOFAR_peel','predictReg')  # default None
sourcedb = parset.get('model','sourcedb')
fits_model = parset.get('model','fits_model')
bl2flag = parset.get('flag','stations')

# Some script-wide definitions
uvlambdamin = 30 # why not 100? What is used in ddf-pipeline?
t_int = 8
final_chan_per_sb = 2
pixscale = '1.5arcsec'
imgsize = 1600
#############################################################################

def solve_and_apply(MSs_object, suffix, sol_factor_t=1, sol_factor_f=1, column_in='DATA'):
    """
    Run a scalarphase-solve and a fulljones-solve on MSs_object and apply the solutions to the data each time.
    Parameters
    ----------
    MSs_object
    suffix: string, suffix for the logs, i.e. 'c00' will give us xxx-c00.log
    sol_factor_t
    sol_factor_f
    column_in: string, optional. default = DATA
    """
    with w.if_todo(f'solve_iono_{suffix}'):
        logger.info('Solving scalarphase...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn={column_in} sol.h5parm=$pathMS/iono.h5 '
                      f'sol.mode=scalarcomplexgain sol.nchan=1 sol.smoothnessconstraint={sol_factor_f * 1e6} sol.solint={sol_factor_t:d} '
                      f'sol.uvlambdamin={uvlambdamin} sol.usemodelcolumn=True',
                      log=f'$nameMS_sol_iono-{suffix}.log', commandType="DP3")

        lib_util.run_losoto(s, f'iono-{suffix}', [ms + '/iono.h5' for ms in MSs_object.getListStr()], \
                            [   parset_dir + '/losoto-flag.parset',
                                parset_dir + '/losoto-plot-scalaramp.parset',
                                parset_dir + '/losoto-plot-scalarph.parset'])

        move(f'cal-iono-{suffix}.h5', 'peel/solutions/')
        move(f'plots-iono-{suffix}', 'peel/plots/')

        # Correct all CORRECTED_DATA -> CORRECTED_DATA2
    with w.if_todo(f'cor_iono_{suffix}'):
        logger.info('Scalarphase correction...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={column_in} msout.datacolumn=CORRECTED_DATA2 cor.updateweights=False '
                      f'cor.parmdb=peel/solutions/cal-iono-{suffix}.h5 cor.correction=phase000', \
                      log=f'$nameMS_cor_iono-{suffix}.log', commandType='DP3')
        ### DONE

    with w.if_todo(f'solve_fulljones_{suffix}'):
        logger.info('Solving full-Jones...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA2 '
                      f'sol.h5parm=$pathMS/fulljones.h5 sol.mode=fulljones sol.nchan=1 sol.solint={sol_factor_t * 128 // t_int:d} '
                      f'sol.smoothnessconstraint={sol_factor_f * 2.0e6} sol.uvlambdamin={uvlambdamin} '
                      f'sol.usemodelcolumn=True', log=f'$nameMS_sol_fulljones-{suffix}.log',
                      commandType="DP3")

        lib_util.run_losoto(s, f'fulljones-{suffix}', [ms + '/fulljones.h5' for ms in MSs_object.getListStr()], \
                            [  parset_dir + '/losoto-plot-fulljones.parset'])
        move(f'cal-fulljones-{suffix}.h5', 'peel/solutions/')
        move(f'plots-fulljones-{suffix}', 'peel/plots/')

        # Correct gain amp and ph CORRECTED_DATA2 -> CORRECTED_DATA2
    with w.if_todo(f'cor_fulljones_{suffix}'):
        logger.info('Full-Jones correction...')
        MSs_object.run(
            f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.correction=fulljones msin.datacolumn=CORRECTED_DATA2 msout.datacolumn=CORRECTED_DATA2 '
            f'cor.parmdb=peel/solutions/cal-fulljones-{suffix}.h5 cor.soltab=\[amplitude000,phase000\]',
            log=f'$nameMS_cor_gain-{suffix}.log', commandType='DP3')

def corrupt_subtract_testimage(MSs_object, suffix, sol_suffix=None, column_in='DATA'):
    """
    Use solutions of a suffix to corrupt and subtract MODEL_DATA, also create an empty test image.
    (NOTE: for field_8, predict 2 directions, solve 2 dir (also for 3c270), then use MODEL_DATA and m87 sols.
    Solved scalar phase and diagonal, use only diagonal amps for subtraction.)
    Parameters
    ----------
    MSs_object
    suffix
    """
    if sol_suffix is None:
        sol_suffix = suffix
    with w.if_todo(f'subtract-{suffix}'):
        logger.info(f'Scalarphase corruption... ({suffix})')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA '
                      f'msout.datacolumn=CORRUPTED_MODEL_DATA cor.updateweights=False '
                      f'cor.parmdb=peel/solutions/cal-iono-{sol_suffix}.h5 cor.correction=phase000 cor.invert=False',
                      log=f'$nameMS_corrup_iono-{sol_suffix}.log', commandType='DP3')

        logger.info(f'Full-Jones corruption... ({suffix})')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRUPTED_MODEL_DATA '
                      f'msout.datacolumn=CORRUPTED_MODEL_DATA cor.correction=fulljones '
                      f'cor.parmdb=peel/solutions/cal-fulljones-{sol_suffix}.h5 cor.soltab=\[amplitude000,phase000\] '
                      f'cor.invert=False', log=f'$nameMS_corrup_gain-{sol_suffix}.log', commandType='DP3')

        logger.info(f'SET SUBTRACTED_DATA = {column_in} - CORRUPTED_MODEL_DATA ({suffix})')
        MSs_object.addcol('SUBTRACTED_DATA', column_in)
        MSs_object.run(f'taql "UPDATE $pathMS SET SUBTRACTED_DATA = {column_in} - CORRUPTED_MODEL_DATA"',
                      log='$nameMS_taql_subtract.log',
                      commandType='general')

    with w.if_todo(f'correct-subtracted-{suffix}'):
        logger.info(f'Scalarphase correction... ({suffix})')
        MSs_object.run(
            f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA '
            f'msout.datacolumn=SUBTRACTED_DATA cor.updateweights=False '
            f'cor.parmdb=peel/solutions/cal-iono-{sol_suffix}.h5 cor.correction=phase000', \
            log=f'$nameMS_cor_iono_subtracted-{sol_suffix}.log', commandType='DP3')

    with w.if_todo(f'test-image-empty-{suffix}'):
        logger.info('Test empty... (SUBTRACTED_DATA)')
        lib_util.run_wsclean(s, f'wsclean-peel.log', MSs_object.getStrWsclean(), weight='briggs -0.5',
                             data_column='SUBTRACTED_DATA', channels_out=3,
                             name='img/test-empty-' + suffix, scale='2.0arcsec', size=2000, niter=0, nmiter=0,
                             no_update_model_required='', minuv_l=uvlambdamin)

def predict_fits_model(MSs_object, model_basename, stepname='init_model', predict_reg=None, apply_beam=True):
    """
    Predict a fits_model into one or more MS files. Optionally, blank fits model with 'predict_reg'
    Parameters
    ----------
    MSs_object
    model_basename
    stepname
    predict_reg
    """
    with w.if_todo(stepname):
        # addcol MODEL_DATA if not there
        # with pt.table(MSs_object.getListStr()[0]) as t:
        #     if not 'MODEL_DATA' in t.colnames(): # assume true for all
        #         MSs_object.addcol('MODEL_DATA', 'DATA', False)
        # TODO should copy here not to change original model
        if predict_reg:
            for model_img in glob.glob(model_basename + '*-model.fits'):
                lib_img.blank_image_reg(model_img, predict_reg, inverse=True, blankval=0.)
        assert os.path.isfile(model_basename + '-0000-model.fits')
        logger.info(f'Using fits model {model_basename}.')
        n = len(glob.glob(model_basename + '-[0-9]*-model.fits'))
        logger.info('Predict (wsclean: %s - chan: %i)...' % ('model', n))
        _str = ' -grid-with-beam -use-idg -use-differential-lofar-beam ' if apply_beam else '' # -use-idg '
        # -use-wgridder
        s.add(f'wsclean -predict  -name {model_basename}  -j {s.max_processors} -channels-out {n} {_str} '
              f'{MSs_object.getStrWsclean()}', log=f'wscleanPRE-{stepname}.log', commandType='wsclean', processors='max')
        s.run(check=True)

def predict_sourcedb_model(MSs_object, stepname='init_model', apply_beam=True):
    """
    Predict a sky model into one or more MS files.
    Parameters
    ----------
    MSs_object: object, All_MSs object
    stepname : str, name of step in walker
    apply_beam : bool, default = True. Apply beam corruption during predict?
    """
    with w.if_todo(stepname):
        lib_util.check_rm('peel.skydb')
        os.system('makesourcedb outtype="blob" format="<" in=peel.skymodel out=peel.skydb')
        for MS in MSs_object.getListStr():
            lib_util.check_rm(MS + '/peel.skydb')
            logger.debug('Copy: peel.skydb -> ' + MS)
            copy2('peel.skydb', MS)
        MSs_object.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/peel.skydb pre.usebeammodel={apply_beam}',
                log='$nameMS_predict.log', commandType='DP3')

def do_testimage(MSs_files, datacolumn='CORRECTED_DATA2'):
    # debug image
    imagename = f'img/test-corrected'

    wsclean_params = {
        'scale': f'1.5arcsec',
        'size': 1200,
        'weight': 'briggs -0.6',
        'join_channels': '',
        'fit_spectral_pol': 5,
        'channels_out': 12, # len(MSs.getFreqs()) // 24,
        'minuv_l': uvlambdamin,
        'name': imagename,
        'no_update_model_required': '',
        'mgain': 0.85,
        'multiscale': '',
        'auto_mask': 3.0,
        'auto_threshold': 1.0,
        'baseline_averaging': 10,
        'data_column' : datacolumn
    }
    with w.if_todo('imaging_corrected'):
        logger.info(f'Cleaning corrected...')
        basemask = 'img/mask-corrected'
        if not os.path.exists(basemask):
            logger.info('Create masks...')
            # dummy clean to get image -> mask
            lib_util.run_wsclean(s, 'wsclean-mask-corrected.log', MSs_files.getStrWsclean(), niter=0, channel_range='0 1', no_reorder='',
                                 interval='0 10', name=imagename, scale=wsclean_params['scale'], size=wsclean_params['size'], nmiter=0)
            # create basemask
            copy2(f'{parset_dir}/masks/VirAhba.reg', f'{peelReg.filename}')
            copy2(f'{imagename}-image.fits', f'{basemask}')
            lib_img.blank_image_reg(basemask, peelReg.filename, inverse=True, blankval=0.)
            lib_img.blank_image_reg(basemask, peelReg.filename, inverse=False, blankval=1.)
            logger.info('Cleaning corrected...')
            lib_util.run_wsclean(s, f'wsclean-corrected.log', MSs_files.getStrWsclean(), niter=1000000,
                                 multiscale_scales='0,20,40,80,160,320', fits_mask=basemask, **wsclean_params)
            os.system(f'cat logs/wsclean-corrected.log | grep "background noise"')

#############################################################################
# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')
    # here images, models, solutions for each group will be saved
    lib_util.check_rm('peel')
    os.makedirs('peel/solutions')
    os.makedirs('peel/plots')
    os.makedirs('peel/masks')
### DONE

MSs = lib_ms.AllMSs( glob.glob(data_dir + '/*.MS'), s, check_flags=False)
if not MSs.isHBA:
    logger.error('Only HBA measurement sets supported for now.')
    sys.exit(1)
try:
    MSs.print_HAcov()
except:
    logger.error('Problem with HAcov, continue anyway.')
logger.info(f'Resolution {MSs.resolution}\'\', using uvlambdamin={uvlambdamin}, averaging to {final_chan_per_sb}chan/SB in freq and to an '
            f'integration time of t_int={t_int}')

peelMask = 'peel/masks/peelmask.fits'
phasecentre = MSs.getListObj()[0].getPhaseCentre()
# region must be a list of ds9 circles and polygons (other shapes can be implemented in lib_util.Rgion_helper()
peelReg = lib_util.Region_helper(peelReg)
# get Pointing centre from either fits model or sourcedb
if fits_model != '':
    assert os.path.isfile(fits_model + '-0000-model.fits')
    model_hdr = fits.open(glob.glob(fits_model + '-[0-9]*-model.fits')[0])[0].header
    model_centre = np.array([model_hdr['CRVAL1'], model_hdr['CRVAL2']])
elif sourcedb != '':
    skymodel = lsm.load(sourcedb)
    skymodel.group('single', method='wmean')
    patch_pos = skymodel.getPatchPositions()['Patch']
    skymodel.write('peel.skymodel', clobber=True)
    model_centre = np.array([patch_pos[0].value, patch_pos[1].value])
else:
    raise ValueError('Please provide either fits_model or sourcedb in the lilf.config [model] section.')

pointing_distance = lib_util.distanceOnSphere(*model_centre, *phasecentre)
logger.info(f"Distance between model and MS phase centre: {pointing_distance:.5f}deg")
##################################################
with w.if_todo('apply'):
    # Find solutions to apply
    cal_dirs = glob.glob(cal_dir + 'id*_-_3[c|C]196') + glob.glob(cal_dir + 'id*_-_3[c|C]295')
    if len(cal_dirs) == 0:
        logger.error(f'No calibrators found in cal dir: {cal_dir}')
        sys.exit(1)

    cal_times = []  # mean times of the cal
    for cal in cal_dirs:
        cal_ms = lib_ms.MS(glob.glob(f'{cal}/*.MS')[0]) # assume all MS are equal
        assert cal_ms.isCalibrator()
        cal_times.append(np.mean(cal_ms.getTimeRange()))
    obs_time = np.mean(MSs.getListObj()[0].getTimeRange())
    delta_t = np.abs(obs_time - np.array(cal_times))
    cal_id = np.argmin(delta_t)
    cal_dir =  cal_dirs[cal_id]
    if delta_t[cal_id] < 5*3600:
        logger.info(f'Found calibrator dir {cal_dirs[cal_id]} which is {delta_t[cal_id]/3600:.2f}h from mean observation time.')
    else:
        logger.error(f'Found calibrator dir {cal_dirs[cal_id]} which is {delta_t[cal_id]/3600:.2f}h from mean observation time!')

    logger.info('Calibrator directory: %s' % cal_dir)
    h5_pa = cal_dir + '/cal-pa.h5'
    h5_amp = cal_dir + '/cal-amp.h5'
    h5_iono = cal_dir + '/cal-iono.h5'
    if not os.path.exists(h5_pa) or not os.path.exists(h5_amp) or not os.path.exists(h5_iono):
        logger.error("Missing solutions in %s" % cal_dir)
        sys.exit()

    # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
    logger.info('Apply solutions (pa)...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA '
            f'msout.storagemanager=dysco cor.parmdb={h5_pa} cor.correction=polalign',
            log='$nameMS_cor1.log', commandType='DP3')

    # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
    logger.info('Apply solutions (amp/ph)...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
            f'msout.datacolumn=CORRECTED_DATA msout.storagemanager=dysco '
            f'cor.steps=[amp,clock] cor.amp.parmdb={h5_amp} cor.amp.correction=amplitudeSmooth '
            f'cor.amp.updateweights=True cor.clock.parmdb={h5_iono} cor.clock.correction=clockMed000',
            log='$nameMS_cor2.log', commandType='DP3')

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
    logger.info('Beam correction...')
    MSs.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
            f'msout.datacolumn=CORRECTED_DATA  msout.storagemanager=dysco '
            f'corrbeam.updateweights=True', log='$nameMS_beam.log',
            commandType='DP3')
### DONE

###################################################

with w.if_todo('clip_ateam'):
    logger.info('Clip A-Team: predict...')
    clip_model = os.path.dirname(__file__) + '/../models/A-Team_lowres.skydb'
    MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
            f'pre.sourcedb={clip_model} pre.sources=[TauA,CasA,CygA]',
            log='$nameMS_pre_clipAteam.log', commandType='DP3')

    logger.info('Clip A-Team: flagging...')
    MSs.run('Ateamclipper.py $pathMS', log='$nameMS_ateamclipper.log', commandType='python')

    MSs.run('plot_Ateamclipper.py logs/Ateamclipper.txt peel/plots/Ateamclipper.png', log='$nameMS_ateamclipper.log', commandType='python')

# Initial flagging
with w.if_todo('flag'):
    logger.info('Flagging...')
    MSs.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS ant.baseline=\"{bl2flag}\" msin.datacolumn=CORRECTED_DATA '
            f'aoflagger.strategy={parset_dir}/HBAdefaultwideband.lua uvmin.uvlambdamin={uvlambdamin}',
            log='$nameMS_flag.log', commandType='DP3')
    logger.info('Remove bad timestamps...')
    MSs.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    # need to make this work for more than one MS!
    # logger.info('Plot weights...')
    # MSs.run(f'reweight.py $pathMS -v -p -a "CS001HBA0"',
    #         log='$nameMS_weights.log', commandType='python')
    # os.system('move *.png self/plots')
### DONE
# Also cut frequencies above 168 MHz such that the number of channels is multiple of 4
with pt.table(MSs.getListObj()[0].pathMS + "/OBSERVATION") as tab:
    field = tab[0]["LOFAR_TARGET"][0]

################################################################################################
################################################################################################
# Start peeling - Two different cases depending on MS - M87 separation
sol_factor_f = 1 if pointing_distance < 3 else 2
if 1/3600 < pointing_distance: # CASE 1 -> model and MS not aligned, peel
    if not os.path.exists('mss-shift'):
        timeint_init = MSs.getListObj()[0].getTimeInt()
        nchan_init = len(MSs.getFreqs())
        nchan = np.sum(np.array(MSs.getFreqs()) < 168.3e6)  # only use 120-168 MHz
        logger.info(f'{nchan_init} channels, {nchan} of which are above 168MHz')
        # nchan = (48*freqstep) * (nchan // (48*freqstep)) # multiple of 48 after average
        lib_util.check_rm('mss-shift')
        os.makedirs('mss-shift')
        logger.info(f'Phase shifting @{timeint_init}s, using {nchan} channels (from {nchan_init})')

        for MS in MSs.getListObj():
            nchan_thisms = int(np.sum(np.array(MS.getFreqs()) < 168.3e6))  # only use 120-168 MHz
            if nchan_thisms == 0:
                logger.warning(f"Skipping {MS.nameMS} (above 168.3MHz)")
                continue
            commandCurrent = MS.concretiseString(
                f'DP3 {parset_dir}/DP3-shift.parset msin=$pathMS msout=mss-shift/$nameMS.MS msin.datacolumn=CORRECTED_DATA '
                f'msin.nchan={nchan_thisms} shift.phasecenter=[{model_centre[0]}deg,{model_centre[1]}deg]')
            logCurrent = MS.concretiseString('$nameMS_initshift.log')
            s.add(cmd=commandCurrent, log=logCurrent, commandType='DP3', )
        s.run(check=True)
    MSs_shift = lib_ms.AllMSs( glob.glob('mss-shift/*.MS'), s, check_flags=False )

    with w.if_todo('apply_beam'):
        logger.info('Correcting beam: DATA -> DATA...')
        MSs_shift.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.invert=True', log='$nameMS_beam.log',
                      commandType='DP3')

    name_msavg = field + '-combined.MS'
    if not os.path.exists('mss-shiftavg'):
        timeint = MSs_shift.getListObj()[0].getTimeInt()
        t_avg_factor = round(16/timeint)
        nchan_shift = len(MSs_shift.getFreqs())
        nchan_shiftavg = len(MSs_shift.getListObj()) / 2 # 0.5 chan/SB
        f_avg_factor = round(nchan_shift/nchan_shiftavg)
        lib_util.check_rm('mss-shiftavg')
        os.makedirs('mss-shiftavg')
        # Average all MSs to one single MS (so that we can average to less than 1 chan/SB)
        logger.info(f'Averaging {timeint:.2f}s -> {t_avg_factor*timeint:.2f}s; {nchan_shift}chan -> {nchan_shiftavg}chan')
        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin=[{",".join(MSs_shift.getListStr())}] msin.datacolumn=DATA '
              f'msout=mss-shiftavg/{name_msavg} avg.timestep={t_avg_factor} avg.freqstep={f_avg_factor} numthreads={s.max_processors}',
              log=f'create-shiftavg.log', commandType='DP3')
        s.run(check=True)
    MSs_shiftavg = lib_ms.AllMSs(glob.glob(f'mss-shiftavg/{name_msavg}'), s, check_flags=False )

    # Add model to MODEL_DATA
    if fits_model != '':
        predict_fits_model(MSs_shiftavg, fits_model)
    else:
        predict_sourcedb_model(MSs_shiftavg)
    #####################################################################################################
    # Get mask -> required to blank model for DI calibration later on
    if not os.path.exists(peelMask):
        logger.info('Create mask...')
        # dummy clean to get image -> mask
        lib_util.run_wsclean(s, 'wsclean-mask.log', MSs_shiftavg.getStrWsclean(), niter=0, channel_range='0 1',
                             no_reorder='', interval='0 10', name='img/mask', scale=pixscale, size=imgsize, nmiter=0)
        # create peelMask
        copy2(f'img/mask-dirty.fits', f'{peelMask}')
        lib_img.blank_image_reg(peelMask, peelReg.filename, inverse=True, blankval=0.)
        lib_img.blank_image_reg(peelMask, peelReg.filename, inverse=False, blankval=1.)
    ##################################################################################
    # solve+apply scalarphase and fulljones
    solve_and_apply(MSs_shiftavg, field, sol_factor_f=sol_factor_f)
    # do_testimage(MSs_shiftavg)
    # corrupt and subtract with above solutions
    # corrupt_subtract_testimage(MSs_shiftavg, field) # testing...
    ##################################################################################
    # Predict model to MSs_shift
    if fits_model != '':
        predict_fits_model(MSs_shift, fits_model, stepname='predict_fnal', predict_reg=predictReg)
    else:
        predict_sourcedb_model(MSs_shift, stepname='predict_fnal')
    # Subtract the model
    corrupt_subtract_testimage(MSs_shift, field) # --> SUBTRACT_DATA
    MSs_subtracted = MSs_shift
# #### CASE 2: Model and MS are distant. Demix
# elif pointing_distance > 5:
#     logger.info(f'Large separation ({pointing_distance:2f}deg). Demixing using {demix_skymodel} (assuming intrinsic fluxes).')
#     if not os.path.exists(peelMask):
#         logger.info('Create mask...')
#         # dummy clean to get image -> mask
#         lib_util.run_wsclean(s, 'wsclean-mask.log', MSs.getStrWsclean(), niter=0, channel_range='0 1',
#                              no_reorder='', interval='0 10', name='img/mask', scale=pixscale, size=imgsize, nmiter=0)
#         # create peelMask
#         copy2(f'img/mask-dirty.fits', f'{peelMask}')
#         lib_img.blank_image_reg(peelMask, peelReg.filename, inverse=True, blankval=0.)
#         lib_img.blank_image_reg(peelMask, peelReg.filename, inverse=False, blankval=1.)
#
#     sm_combined = lsmtool.load(demix_skymodel)
#     demix_patches = sm_combined.getPatchNames()
#     if not os.path.exists('demix_combined.skydb'):
#         # get a rough field model (no beam applied)
#         os.system(f'wget -O demix_field.skymodel \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv4.cgi?coord={phasecentre[0]},'
#                   f'{phasecentre[1]}&radius={MSs.getListObj()[0].getFWHM("min") / 2.}&unit=deg&deconv=y\' ')  # This assumes first MS is lowest in freq
#         sm_field = lsmtool.load('demix_field.skymodel', MSs.getListStr()[0])
#         sm_field.remove('I<0.1', applyBeam=True) # remove faint sources
#         sm_field.remove(f'{peelMask} == True')  # This should remove the peel source from the field model.
#         sm_field.group('single', root='field')
#         sm_field.setColValues('LogarithmicSI', ['true']*len(sm_field))
#         sm_combined.concatenate(sm_field)
#         sm_combined.write('demix_combined.skymodel', clobber=True)
#         lib_util.check_rm('demix_combined.skydb')
#         os.system('makesourcedb outtype="blob" format="<" in=demix_combined.skymodel out=demix_combined.skydb')
#         for MS in MSs.getListStr():
#             lib_util.check_rm(MS + '/demix_combined.skydb')
#             logger.debug('Copy: demix_combined.skydb -> ' + MS)
#             copy2('demix_combined.skydb', MS)
#
#     logger.debug(f'Demix sources: {demix_patches}')
#     with w.if_todo('demix'): # Demix CORRECTED_DATA -> SUBTRACTED_DATA
#         MSs.run(f'DP3 {parset_dir}/DP3-demix.parset msin=$pathMS demixer.skymodel=$pathMS/demix_combined.skydb '
#                 f'demixer.instrumentmodel=$pathMS/instrument_demix '
#                 f'demixer.subtractsources=\[{",".join(demix_patches)}\]',
#                 log='$nameMS_demix.log', commandType='DP3')
#     MSs_subtracted = MSs
#### CASE 2: Model and MS are aligned. Solve
else:
    if fits_model != '':
        predict_fits_model(MSs, fits_model, apply_beam=False)
    else:
        predict_sourcedb_model(MSs)
    solve_and_apply(MSs, field, column_in='CORRECTED_DATA')
    # delete
    # TODO if col exists
    # MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CORRECTED_DATA2"',
    #         log='$nameMS_taql_delcol1.log', commandType='general')
    ### debug:
    # do_testimage(MSs)
    corrupt_subtract_testimage(MSs, field, column_in='CORRECTED_DATA')  # --> SUBTRACTED_DATA
    MSs_subtracted = MSs
    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CORRECTED_DATA, CORRUPTED_MODEL_DATA"',
            log='$nameMS_taql_delcol2.log', commandType='general')

###################################################################################################################
###################################################################################################################

if not os.path.exists('mss-peel'):
    timeint_init = MSs_subtracted.getListObj()[0].getTimeInt()
    avgtimeint = int(round(t_int / timeint_init))
    nchan_init = len(MSs_subtracted.getFreqs())
    nchan = np.sum(np.array(MSs_subtracted.getFreqs()) < 168.3e6)  # only use 120-168 MHz
    chan_per_sb = int(nchan_init / len(MSs_subtracted.getListStr()))
    freqstep =  int(chan_per_sb / final_chan_per_sb)
    lib_util.check_rm('mss-peel')
    os.makedirs('mss-peel')
    logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (
        timeint_init, timeint_init * avgtimeint, nchan_init, nchan // freqstep))
    logger.info('Concatenating MSs in groups of 10, fill holes.')
    for sb_start in np.arange(104, 354, 10):
        base_msname = MSs_subtracted.getListStr()[0].split('SB')[0]
        list_mss = [f'{base_msname}SB{sb_start+i}.MS' for i in range(10)]
        freq_group = int((121189880 + (sb_start-104) * 195312.5) * 1e-6)
        out_msname = f'mss-peel/{base_msname.split("/")[-1]}{freq_group}MHz.ms'
        s.add(f'DP3 {parset_dir}/DP3-shiftavg.parset msin=[{",".join(list_mss)}] msout={out_msname} msin.datacolumn=SUBTRACTED_DATA '
                    f'msin.nchan={chan_per_sb*10} avg.timestep={avgtimeint} avg.freqstep={freqstep} shift.phasecenter=[{phasecentre[0]}deg,{phasecentre[1]}deg] ',
                    log=f'{freq_group}MHz-initavg.log', commandType='DP3')
    s.run(check=True)
    os.system("rename.ul .MS .ms mss-peel/*")
MSs_peel = lib_ms.AllMSs(glob.glob('mss-peel/*.ms'), s, check_flags=False)

with w.if_todo('peel-corrbeam'):
    # check if beam needs to be corrected in phase center - assume all MSs have same Beam Keyword.
    do_corrbeam = True
    with pt.table(MSs_peel.getListStr()[-1]) as t:
        phasecentre_rad = np.squeeze(t.FIELD.getcol("PHASE_DIR"))
        kws = t.getcolkeywords('DATA')
        try:
            if kws['LOFAR_APPLIED_BEAM_MODE'] == 'Full':
                print(kws['LOFAR_APPLIED_BEAM_DIR']['m0'], kws['LOFAR_APPLIED_BEAM_DIR']['m1'], phasecentre_rad)
                dist = ((kws['LOFAR_APPLIED_BEAM_DIR']['m0']['value'] - phasecentre_rad[0])**2 +
                        (kws['LOFAR_APPLIED_BEAM_DIR']['m1']['value'] - phasecentre_rad[1])**2)**0.5
                print(dist)
                if dist < 1/3600:
                    do_corrbeam = False
        except KeyError:
            logger.warning('Beam Keywords not found.')

    if do_corrbeam:
        logger.info('Correcting beam: DATA -> DATA...')
        MSs_peel.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS', log='$nameMS_beam.log', commandType='DP3')

# If distant pointing -> derive and apply new phase-solutions against field model
# if pointing_distance > 1.0:
#     with w.if_todo('get_field_model'):
#         os.system(f'wget -O field.skymodel \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv4.cgi?coord={phasecentre[0]},'
#                   f'{phasecentre[1]}&radius={MSs_peel.getListObj()[0].getFWHM("min") / 2.}&unit=deg&deconv=y\' ')  # This assumes first MS is lowest in freq
#         lsm = lsmtool.load('field.skymodel', MSs_peel.getListStr()[0])
#         lsm.remove('I<0.1', applyBeam=True)
#         lsm.remove(f'{peelMask} == True')  # This should remove the peel source from the field model.
#         lsm.group('single', root='field')
#         lsm.write('field.skymodel', clobber=True)
#         lib_util.check_rm('field.skydb')
#         os.system('makesourcedb outtype="blob" format="<" in=field.skymodel out=field.skydb')
#         for MS in MSs_peel.getListStr():
#             lib_util.check_rm(MS + '/field.skydb')
#             logger.debug('Copy: field.skydb -> ' + MS)
#             copy2('field.skydb', MS)
#
#     with w.if_todo('sol_di'): # Solve DATA vs. MODEL_DATA
#         logger.info('Solving direction-independent (scalarphase)...')
#         MSs_peel.run(
#             f'DP3 {parset_dir}/DP3-sol-sourcedb.parset msin=$pathMS msin.datacolumn=DATA '
#             f'sol.h5parm=$pathMS/di.h5 sol.mode=scalarphase sol.uvlambdamin={uvlambdamin} sol.sourcedb=$pathMS/field.skydb',
#             log='$nameMS_sol_di.log', commandType="DP3")
#
#         lib_util.run_losoto(s, f'di', [ms + '/di.h5' for ms in MSs_peel.getListStr()], \
#                             [ # parset_dir+'/losoto-flag.parset',
#                               # parset_dir + '/losoto-plot-scalaramp.parset',
#                               parset_dir + '/losoto-plot-scalarph.parset'])
#         move(f'cal-di.h5', 'peel/solutions/')
#         move(f'plots-di', 'peel/plots/')
#
#     # Correct DATA -> CORRECTED_DATA
#     with w.if_todo('cor_di'):
#         logger.info('Direction-independent correction...')
#         MSs_peel.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.updateweights=False '
#                       f'cor.parmdb=peel/solutions/cal-di.h5 cor.correction=phase000', \
#                       log=f'$nameMS_cor_di.log', commandType='DP3')
#
with w.if_todo(f'test-image-wide'):
    logger.info('Cleaning wide...') # IMAGING - either on DATA or if present, CORRECTED_DATA
    lib_util.run_wsclean(s, f'wsclean-wide.log', MSs_peel.getStrWsclean(), weight='briggs -0.0',
                         name='img/peel-wide',  parallel_deconvolution=1024, scale='4.0arcsec', size=5000, niter=500000,
                         nmiter=8, minuv_l=uvlambdamin, channels_out=3, join_channels='',  mgain=0.9, auto_threshold=1.0, auto_mask=4.0,
                         no_update_model_required='', do_predict=False, local_rms='')
    os.system(f'cat logs/wsclean-wide.log | grep "background noise"')

logger.info("Done.")
