#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# pipeline for a bright a-team source in LOFAR HBA
# they need to be in "./mss/"

import sys, os, glob, re
from shutil import copy2, copytree, move
import numpy as np
import casacore.tables as pt
import lsmtool as lsm

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-virgo.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-virgo.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_virgo','parset_dir')
cal_dir = parset.get('LOFAR_virgo','cal_dir')
data_dir = parset.get('LOFAR_virgo','data_dir')
init_model = parset.get('model','fits_model')
userReg = parset.get('model','userReg')
bl2flag = parset.get('flag','stations')

#############################################################################
m87_model = '/beegfs/p1uy068/virgo/models/m87/m87'
m87_model_hires = '/beegfs/p1uy068/virgo/models/m87/m87'
field_model = '/beegfs/p1uy068/virgo/models/m87_field/'
# TODO
# Station/antennaconstrains? For FR just all CS or maybe only superterp or so? For phases any constraint?
# updateweights
# amp normalization -> fulljones???

# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')
    # here images, models, solutions for each group will be saved
    lib_util.check_rm('self')
    os.makedirs('self/solutions')
    os.makedirs('self/plots')
    os.makedirs('self/masks')
### DONE

MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*.MS'), s, check_flags=False )
if not MSs.isHBA:
    logger.error('Only HBA measurement sets supported.')
    sys.exit(1)
try:
    MSs.print_HAcov()
except:
    logger.error('Problem with HAcov, continue anyway.')

is_IS = False
if MSs.resolution < 2.0:
    is_IS = True
if is_IS:
    uvlambdamin = 100
    uvlambdamin_sol = 500
    freqstep = 2  # might wanna decrease later!
    t_int = 8 # might wanna decrease later!
    final_chan_per_sb = 9999
else:
    uvlambdamin = 100
    uvlambdamin_sol = 100
    t_int = 8
    final_chan_per_sb = 2

logger.info(f'Resolution {MSs.resolution}\'\', using uvlambdamin={uvlambdamin}, averaging to {final_chan_per_sb}chan/SB in freq and to an '
            f'integration time of t_int={t_int}')

with pt.table(MSs.getListObj()[0].pathMS + '/OBSERVATION', ack=False) as t:
    if 'M87' in t.getcell('LOFAR_TARGET',0):
        target = 'VirA'
    else:
        raise ValueError('Target is not M87. Only M87 is implemented.')

basemask = 'self/masks/basemask.fits'
basemaskC = 'self/masks/basemaskC.fits'
baseregion = 'self/masks/baseregion.reg'
baseregionC = 'self/masks/baseregionC.reg'
phasecentre = MSs.getListObj()[0].getPhaseCentre()
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

    for MS in MSs.getListObj():
        obs_time = np.mean(MS.getTimeRange())
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

        # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign, amp, ph, beam + reweight)
        if 'LBA' in MS.getAntennaSet():
            s.add(MS.concretiseString(f'DP3 {parset_dir}/DP3-corcal.parset msin=$pathMS cor.pa.parmdb={h5_pa} '
                                      f'cor.amp.parmdb={h5_amp} cor.ph.parmdb={h5_iono} cor.ph.correction=phaseOrig000'),
                  log=MS.concretiseString('$nameMS_corcal.log'),
                  commandType='DP3')
        elif 'HBA' in MS.getAntennaSet():
            s.add(MS.concretiseString(f'DP3 {parset_dir}/DP3-corcal.parset msin=$pathMS cor.pa.parmdb={h5_pa} '
                                      f'cor.amp.parmdb={h5_amp} cor.ph.parmdb={h5_iono} cor.ph.correction=clockMed000'),
                  log=MS.concretiseString('$nameMS_corcal.log'),
                  commandType='DP3')

    logger.info('Running - apply solutions (pa, amp, ph/clock, beam)')
    s.run(check=True)
    ### DONE

#################################################################################################
# # Clip A-Team
# with w.if_todo('clip_ateam'):
#     logger.info('Clip A-Team: predict...')
#     clip_model = os.path.dirname(__file__) + '/../models/A-Team_lowres.skydb'
#     MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
#             f'pre.sourcedb={clip_model} pre.sources=[TauA,CasA,CygA]',
#             log='$nameMS_pre_clipAteam.log', commandType='DP3')
#
#     logger.info('Clip A-Team: flagging...')
#     MSs.run('Ateamclipper.py $pathMS', log='$nameMS_ateamclipper.log', commandType='python')
#     MSs.run('plot_Ateamclipper.py logs/Ateamclipper.txt peel/plots/Ateamclipper.png', log='$nameMS_ateamclipper.log', commandType='python')

# Initial flagging
with w.if_todo('flag'):
    logger.info('Flagging...')
    MSs.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS  ant.baseline=\"{bl2flag}\" msin.datacolumn=CORRECTED_DATA '
            f'aoflagger.strategy={parset_dir}/HBAdefaultwideband.lua uvmin.uvlambdamin={uvlambdamin}',
            log='$pathMS_flag.log', commandType='DP3')
    logger.info('Remove bad timestamps...')
    MSs.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    # TODO need to make this work for more than one MS!
    # logger.info('Plot weights...')
    # MSs.run(f'reweight.py $pathMS -v -p -a "CS001HBA0"',
    #         log='$nameMS_weights.log', commandType='python')
    # os.system('move *.png self/plots')
### DONE

if not os.path.exists('mss-avg'):
    timeint_init = MSs.getListObj()[0].getTimeInt()
    nchan_init = len(MSs.getFreqs())
    nchan = np.sum(np.array(MSs.getFreqs()) < 168.3e6)  # only use 120-168 MHz
    logger.info(f'{nchan_init} channels, {nchan} of which are above 168MHz')
    # nchan = (48*freqstep) * (nchan // (48*freqstep)) # multiple of 48 after average
    lib_util.check_rm('mss-avg')
    os.makedirs('mss-avg')
    logger.info(f'Averaging @{timeint_init}s, using {nchan} channels (from {nchan_init})')

    for MS in MSs.getListObj():
        nchan_thisms = int(np.sum(np.array(MS.getFreqs()) < 168.3e6))  # only use 120-168 MHz
        if nchan_thisms == 0:
            logger.warning(f"Skipping {MS.nameMS} (above 168.3MHz)")
            continue
        commandCurrent = MS.concretiseString(
            f'DP3 {parset_dir}/DP3-avg.parset msin=$pathMS msout=mss-avg/$nameMS.MS msin.datacolumn=CORRECTED_DATA '
            f'msin.nchan={nchan_thisms} avg.timestep=2 avg.freqstep=4')
        logCurrent = MS.concretiseString('$nameMS_initavg.log')
        s.add(cmd=commandCurrent, log=logCurrent, commandType='DP3', )
    s.run(check=True)
MSs = lib_ms.AllMSs( glob.glob('mss-avg/*.MS'), s, check_flags=False )


# Add model to MODEL_DATA
with w.if_todo('init_model'):
    model = m87_model_hires if is_IS else m87_model
    n = len(glob.glob(f'{model}-[0-9]*-model.fits'))
    logger.info('Predict (wsclean: %s - chan: %i)...' % (model, n))
    s.add(f'wsclean -predict -no-reorder -name {model} -j {s.max_processors} -use-wgridder -channels-out {n} {MSs.getStrWsclean()}',
          log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    s.run(check=True)
    # else:
    #     copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
        # sourcedb_basename = init_model.split('/')[-1]
        # for MS in MSs.getListStr():
        #     lib_util.check_rm(MS + '/' + sourcedb_basename)
        #     logger.debug('Copy: ' + init_model + ' -> ' + MS)
        #     copy2(init_model, MS)
        #
        # note: do not add MODEL_DATA or the beam is transported from DATA, while we want it without beam applied
        # logger.info('Predict (DP3: %s - %s))...' % (sourcedb_basename, target))
        # MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.usebeammodel=false pre.sources={target} '
        #         f'pre.sourcedb=$pathMS/{sourcedb_basename}', log='$nameMS_pre.log', commandType='DP3')
with w.if_todo('addcol-residual'):
    MSs.addcol('RESIDUAL_DATA', 'DATA')
#####################################################################################################
# Self-cal cycle
# field_subtracted = False
for c in range(100):
    with w.if_todo(f'solve_iono_c{c:02}'):
        logger.info('Solving scalarphase...')
        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA '
                f'sol.h5parm=$pathMS/iono.h5 sol.mode=scalarcomplexgain sol.nchan=1 sol.solint=1 sol.smoothnessconstraint=1e6 '
                f'sol.uvlambdamin={uvlambdamin_sol}', log=f'$nameMS_sol_iono-c{c:02}.log', commandType="DP3")

        lib_util.run_losoto(s, f'iono-c{c:02}', [ms + '/iono.h5' for ms in MSs.getListStr()], \
                            [#parset_dir + '/losoto-flag.parset', # disabled to not flag too much in many iterations...
                             parset_dir + '/losoto-plot-scalaramp.parset',
                             parset_dir + '/losoto-plot-scalarph.parset'])
        move(f'cal-iono-c{c:02}.h5', 'self/solutions/')
        move(f'plots-iono-c{c:02}', 'self/plots/')

    with w.if_todo(f'cor_iono_c{c:02}'):
        logger.info('Scalarphase correction...')
        MSs.run(
            f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA cor.updateweights=False '
            f'cor.parmdb=self/solutions/cal-iono-c{c:02}.h5 cor.correction=phase000', \
            log=f'$nameMS_cor_iono.log', commandType='DP3')
    ### DONE

    with w.if_todo('solve_gain_c%02i' % c):
        logger.info('Solving full-Jones...')
        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA ' 
                f'sol.h5parm=$pathMS/fulljones.h5 sol.mode=fulljones sol.nchan=10 sol.solint={8*30//t_int} '# sol.nchan=8 4*64//t_int #sol.smoothnessconstraint=1.5e6 '
                f'sol.uvlambdamin={uvlambdamin_sol}', log=f'$nameMS_sol_fulljones-c{c}.log',
                commandType="DP3")

        lib_util.run_losoto(s, f'fulljones-c{c:02}', [ms + '/fulljones.h5' for ms in MSs.getListStr()], \
                            [parset_dir + '/losoto-fulljones.parset'])

        move(f'cal-fulljones-c{c:02}.h5', 'self/solutions/')
        move(f'plots-fulljones-c{c:02}', 'self/plots/')

    # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
    with w.if_todo('cor_gain_c%02i' % c):
        logger.info('Full-Jones correction...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                f'cor.correction=fulljones cor.parmdb=self/solutions/cal-fulljones-c{c:02}.h5 '
                f'cor.soltab=\[amplitude000,phase000\]', log=f'$nameMS_cor_gain-c{c:02}.log', commandType='DP3')
    # clip on residuals
    with w.if_todo('clip_residuals_c%02i' % c):
        logger.info('Not clipping on RESIDUAL_DATA')
        # MSs.run('taql "UPDATE $pathMS SET RESIDUAL_DATA = CORRECTED_DATA-MODEL_DATA"', log=f'$nameMS_residual-c{c:02}.log')
        # for MS in MSs.getListObj():
        #     with pt.table(MS.pathMS, readonly=False) as t:
        #         residuals = t.getcol('RESIDUAL_DATA')
        #         flags = t.getcol('FLAG')
        #         ant1 = t.getcol('ANTENNA1')
        #         ant2 = t.getcol('ANTENNA2')
        #         sigma = np.nanstd(residuals[~flags])
        #         newflags = (np.abs(residuals) > 5 * sigma) | flags
        #         logger.info(f'({MS.nameMS}) Using sigma {sigma:2e}. Flagged data: before {np.sum(flags)/flags.size:.3%}; after {np.sum(newflags)/flags.size:.3%}')
        #         t.putcol('FLAG', newflags)

    ###################################################################################################################
    # clean CORRECTED_DATA
    imagename = f'img/img-c{c:02}'

    wsclean_params = {
        'scale': '0.2arcsec' if is_IS else '1.0arcsec',
        'size': 6000 if is_IS else 2400,
        'weight': 'briggs -0.6' if is_IS  else'briggs -1.0',
        'join_channels': '',
        # 'fit_spectral_pol': 12,
        'channels_out':  12, # len(MSs.getFreqs()) // 24,
        # 'deconvolution_channels': len(MSs.getFreqs()) // 96 if is_IS else len(MSs.getFreqs()) // 48,
        'minuv_l': uvlambdamin,
        'name': imagename,
        'no_update_model_required': '',
        'do_predict': True,
        'use_wgridder': '',
        'wgridder_accuracy': 1e-04 if is_IS else 1e-05,
        # 'use_idg': '',
        # 'grid_with_beam':'',
        # 'idg_mode': 'cpu',
        # 'no_small_inversion': '',
        'mgain': 0.8,
        'multiscale': '',
        # 'save_source_list': '',
        # 'auto_mask': 3.0,
        'auto_threshold': 0.5,
        'baseline_averaging': 6,
    }
    with w.if_todo('imaging_c%02i' % c):
        logger.info(f'Cleaning (cycle: {c}; size: {wsclean_params["size"]}pix scale: {wsclean_params["scale"]})...')
        if not os.path.exists(basemask) or not os.path.exists(basemaskC):
            logger.info('Create masks...')
            # dummy clean to get image -> mask
            lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), niter=0, channel_range='0 1', no_reorder='',
                                 interval='0 10', name=imagename, scale=wsclean_params['scale'], size=wsclean_params['size'], nmiter=0)
            # create basemask
            copy2(f'{parset_dir}/masks/VirAhba.reg', f'{baseregion}')
            copy2(f'{imagename}-image.fits', f'{basemask}')
            copy2(f'{parset_dir}/masks/VirAChba.reg', f'{baseregionC}')
            copy2(f'{imagename}-image.fits', f'{basemaskC}')
            lib_img.blank_image_reg(basemask, baseregion, inverse=True, blankval=0.)
            lib_img.blank_image_reg(basemask, baseregion, inverse=False, blankval=1.)
            lib_img.blank_image_reg(basemaskC, baseregionC, inverse=True, blankval=0.)
            lib_img.blank_image_reg(basemaskC, baseregionC, inverse=False, blankval=1.)
            if userReg != '':
                lib_img.blank_image_reg(basemask, userReg, inverse=True, blankval=1.)

        if not is_IS:
            logger.info('Cleaning...')
            lib_util.run_wsclean(s, f'wsclean-c{c}.log', MSs.getStrWsclean(), niter=1000000,
                                 multiscale_scales='0,10,17,30,50,80,140,200,340', fits_mask=basemask, **wsclean_params)
            os.system(f'cat logs/wsclean-c{c}.log | grep "background noise"')
        else: # International LOFAR Telescope
            logger.info('Cleaning...')
            wsclean_params['mgain'] = 0.5
            lib_util.run_wsclean(s, f'wsclean-c{c}.log', MSs.getStrWsclean(), niter=1000000,
                                 auto_mask= 3.0, fits_mask=basemask, multiscale_scales='0,10,20,40,80,160,320', mem=95,
                                 taper_gaussian='0.8asec', maxuv_l=250000,
                                  **wsclean_params) # taper_gaussian='2.0asec', # 250000 should be a 0.8asec or so
            os.system(f'cat logs/wsclean-c{c}.log | grep "background noise"')


    # widefield_model = False
    # if c >= 5 and not field_subtracted:
    #     with w.if_todo('imaging_wide_c%02i' % c):
    #         logger.info('SET SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA')
    #         MSs.addcol('SUBTRACTED_DATA', 'CORRECTED_DATA')
    #         MSs.run('taql "UPDATE $pathMS SET SUBTRACTED_DATA = CORRECTED_DATA-MODEL_DATA"', log='$nameMS_taql_subtract.log', commandType='general')
    #
    #         logger.info('Cleaning Virgo A subtracted wide-field image...')
    #         lib_util.run_wsclean(s, f'wsclean-wide-c{c}.log', MSs.getStrWsclean(), weight='briggs -0.5', data_column='SUBTRACTED_DATA',
    #                              name=imagename+'-wide', parallel_deconvolution=1024, scale='2.0arcsec', size=9000, niter=500000,
    #                              join_channels='', channels_out=6, nmiter=15, fit_spectral_pol=3, minuv_l=uvlambdamin, multiscale='', multiscale_max_scales=5,
    #                              mgain=0.85, auto_threshold=1.0, auto_mask=4.0, baseline_averaging='', no_update_model_required='', do_predict=True, local_rms='',
    #                              fits_mask=parset_dir+'/masks/FieldMask.fits')
    #         widefield_model = True
    #         os.system(f'cat logs/wsclean-wide-c{c}.log | grep "background noise"')

            # logger.info('Cleaning Virgo A subtracted low-res wide-field image...')
            # TODO eventuall 5arcsec lowrs
            # lib_util.run_wsclean(s, f'wsclean-wide-lr-c{c}.log', MSs.getStrWsclean(), weight='briggs -0.5', taper_gaussian='30arcsec', data_column='SUBTRACTED_DATA',
            #                      name=imagename+'-wide-lr', parallel_deconvolution=2048, scale='6.0arcsec', size=3000, niter=500000,
            #                      join_channels='', channels_out=6, fit_spectral_pol=3, minuv_l=uvlambdamin, multiscale='', multiscale_max_scales=5,
            #                      mgain=0.85, auto_threshold=1.0, auto_mask=3.0, baseline_averaging='', no_update_model_required='')
            # os.system(f'cat logs/wsclean-wide-c{c}.log | grep "background noise"')

        # with w.if_todo('subtract_wide_c%02i' % c):
            # if not widefield_model:
            #     logger.info('Predict rest-field...')
            #     n = len(glob.glob(field_model + 'm87-field-[0-9]*-model.fits'))
            #     logger.info('Predict (wsclean: %s - chan: %i)...' % ('model-field', n))
            #     s.add(f'wsclean -predict -name {field_model}m87-field -j {s.max_processors} -channels-out {n} {MSs.getStrWsclean()}',
            #           log='wscleanPRE-field.log', commandType='wsclean', processors='max')
            #     s.run(check=True)
            #
            # logger.info('TEST EMPTY')
            # logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA')
            # MSs.run('taql "UPDATE $pathMS SET SUBTRACTED_DATA = SUBTRACTED_DATA-MODEL_DATA"', log='$nameMS_taql_subtract_empty.log',
            #         commandType='general')
            # lib_util.run_wsclean(s, f'wsclean-empty-c{c}.log', MSs.getStrWsclean(), weight='briggs -0.5', data_column='SUBTRACTED_DATA',
            #                      name=imagename+'-empty', scale='2.0arcsec', size=7200, minuv_l=uvlambdamin, no_update_model_required='')

            # logger.info('Corrupt widefield MODEL_DATA...')
            # logger.info('Scalarphase corruption (MODEL_DATA)...')
            # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA '
            #         f'cor.updateweights=False cor.parmdb=self/solutions/cal-iono-c{c:02}.h5 cor.correction=phase000 '
            #         f'cor.invert=False', log=f'$nameMS_corrupt_iono-c{c:02}.log', commandType='DP3')
            # logger.info('Full-Jones corruption (MODEL_DATA)...')
            # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA '
            #         f'cor.correction=fulljones cor.parmdb=self/solutions/cal-fulljones-c{c:02}.h5 '
            #         f'cor.soltab=\[amplitude000,phase000\] cor.invert=False', log=f'$nameMS_corrupt_gain-c{c:02}.log', commandType='DP3')
            #
            # logger.info('Set FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA')
            # MSs.run('taql "UPDATE $pathMS SET FR_CORRECTED_DATA = FR_CORRECTED_DATA-MODEL_DATA"', log='$nameMS_taql_subtract.log',
            #         commandType='general')
            #
            # logger.info('BL-smooth...')
            # MSs.run('BLsmooth.py -r -c 4 -n 8 -f 5e-3 -i FR_CORRECTED_DATA -o FR_SMOOTHED_DATA $pathMS',
            #         log='$nameMS_smooth.log', commandType='python', maxThreads=8)
            #
            # logger.info('Get back Virgo A MODEL_DATA...')
            # s.add(f'wsclean -predict -name {imagename} -j {s.max_processors} -channels-out {wsclean_params["channels_out"]} {MSs.getStrWsclean()}',
            #       log='wscleanPRE-field.log', commandType='wsclean', processors='max')
            # s.run(check=True)
logger.info("Done.")

# oXVL51G   # clean CORRECTED_DATA
#     imagename = f'img/img-c{c:02}'
#
#     wsclean_params = {
#         'scale': f'{0.3/5}arcsec' if is_IS else f'1.5arcsec',
#         'size': 3000 if is_IS else 1400,
#         'weight': 'briggs -1.0' if is_IS  else'briggs -0.4',
#         'join_channels': '',
#         'fit_spectral_pol': 3,
#         'channels_out': len(MSs.getFreqs()) // 48 if is_IS else len(MSs.getFreqs()) // 24,
#         # 'deconvolution_channels': len(MSs.getFreqs()) // 96 if is_IS else len(MSs.getFreqs()) // 48,
#         'minuv_l': uvlambdamin,
#         'name': imagename,
#         'no_update_model_required': '',
#         'do_predict': True,
#         'mgain': 0.8,
#         'save_source_list': ''
#         ''
#     }
#     with w.if_todo('imaging_c%02i' % c):
#         logger.info(f'Cleaning (cycle: {c}; size: {wsclean_params["size"]}pix scale: {wsclean_params["scale"]})...')
#
#         if not os.path.exists(basemask) or not os.path.exists(basemaskC):
#             logger.info('Create masks...')
#             # dummy clean to get image -> mask
#             lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), niter=0, channel_range='0 1',
#                                  interval='0 10', name=imagename, scale=wsclean_params['scale'], size=wsclean_params['size'], nmiter=0)
#             # create basemask
#             copy2(f'{parset_dir}/masks/VirAhba.reg', f'{baseregion}')
#             copy2(f'{imagename}-image.fits', f'{basemask}')
#             copy2(f'{parset_dir}/masks/VirAChba.reg', f'{baseregionC}')
#             copy2(f'{imagename}-image.fits', f'{basemaskC}')
#             lib_img.blank_image_reg(basemask, baseregion, inverse=True, blankval=0.)
#             lib_img.blank_image_reg(basemask, baseregion, inverse=False, blankval=1.)
#             lib_img.blank_image_reg(basemaskC, baseregionC, inverse=True, blankval=0.)
#             lib_img.blank_image_reg(basemaskC, baseregionC, inverse=False, blankval=1.)
#             if userReg != '':
#                 lib_img.blank_image_reg(basemask, userReg, inverse=True, blankval=1.)
#
#         if not is_IS:
#             # iterate cocoon - halo
#             last_it, _break = False, False
#             am = dict() # for auto-masking
#             switch_gain = 0.7
#             threshold = 0.2 # initial threshold
#             threshold_final = 0.0003
#             i = 0
#             if os.path.exists(f'{imagename}-sources-merged.txt'):
#                 sm = lsm.load(f'{imagename}-sources-merged.txt')
#
#             while not _break:
#                 logger.info(f'Iteration threshold = {threshold}')
#                 wsclean_params['threshold'] = threshold
#                 # Clean delta scale on cocoon
#                 with w.if_todo(f'imaging_cocoon_c{c:02}-{i:02}'):
#                     logger.info('Cleaning cocoon...')
#                     lib_util.run_wsclean(s, f'wsclean1-c{c}.log', MSs.getStrWsclean(), niter=5000000, fits_mask=basemaskC,
#                                          gain=0.07, **wsclean_params)
#                     os.system(f'cat logs/wsclean1-c{c}.log | grep "background noise"')
#                     if i == 0:
#                         logger.debug('create initial skymodel')
#                         move(f'{imagename}-sources.txt', f'{imagename}-sources-merged.txt')
#                         sm = lsm.load(f'{imagename}-sources-merged.txt')
#                     else:
#                         logger.debug('append skymodel')
#                         sm.concatenate(f'{imagename}-sources.txt')
#                         os.remove(f'{imagename}-sources.txt')
#                         sm.write(f'{imagename}-sources-merged.txt', clobber=True)
#                 try:
#                     del wsclean_params['mgain'] # only use C-S during first call where many orders of CC are found
#                 except KeyError: pass
#                 wsclean_params['cont'] = True
#                 wsclean_params['reuse-psf'] = imagename
#                 # Clean multi-scale halo
#                 if threshold < 4 * threshold_final:
#                     am['auto_mask'] = 3.0
#                 with w.if_todo(f'imaging_extended_c{c:02}-{i:02}'):
#                     logger.info('Cleaning (multi-scale) halo...')
#                     lib_util.run_wsclean(s, f'wsclean2-c{c}.log', MSs.getStrWsclean(), niter=400000, multiscale='',
#                                          multiscale_shape='gaussian', #multiscale_convolution_padding=1.2,
#                                          multiscale_scales='0,20,40,80,160,320', fits_mask=basemask, **{**am, **wsclean_params})
#                     os.system(f'cat logs/wsclean2-c{c}.log | grep "background noise"')
#
#                     logger.debug('append skymodel')
#                     sm.concatenate(f'{imagename}-sources.txt')
#                     os.remove(f'{imagename}-sources.txt')
#                     sm.write(f'{imagename}-sources-merged.txt', clobber=True)
#
#                 i += 1
#                 threshold = (1.0 - switch_gain) * threshold
#                 if last_it:
#                     _break = True # normal break messes with "with"
#                 elif threshold < threshold_final: # perform last iteration using the minimal threshold
#                     threshold = threshold_final
#                     last_it = True