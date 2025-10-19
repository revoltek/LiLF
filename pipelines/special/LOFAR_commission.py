#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation .
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import os, glob, re
import casacore.tables as pt
import numpy as np

########################################################
from LiLF import lib_ms, lib_util, lib_log, lib_h5

logger_obj = lib_log.Logger('pipeline-cal')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry=False)
w = lib_util.Walker('pipeline-cal.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: ' + str(dict(parset['LOFAR_cal'])))
parset_dir = parset.get('LOFAR_cal', 'parset_dir')
data_dir = parset.get('LOFAR_cal', 'data_dir')
skymodel = parset.get('LOFAR_cal', 'skymodel')
imaging = parset.getboolean('LOFAR_cal', 'imaging')
fillmissingedges = parset.getboolean('LOFAR_cal', 'fillmissingedges')
less_aggressive_flag = parset.getboolean('LOFAR_cal', 'less_aggressive_flag') # change flagging to hande data that uses only alternating sb
develop = parset.getboolean('LOFAR_cal', 'develop') # for development, don't delete files
use_shm = parset.getboolean('LOFAR_cal', 'use_shm') # use shared memory for wsclean
bl2flag = parset.get('flag', 'stations')
use_GNSS = parset.getboolean('LOFAR_cal', 'use_GNSS') # Use GNSS for pre-TEC and FR
#############################################################

#############################################################
# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('plots-preiono plots-preiono-cs plots-fr plots-bp plots-bp-sub plots-fj plots-iono plots-iono-cs plots-pa plots-amp plots-weights plots-test*}')
    lib_util.check_rm('cal*.h5')
    lib_util.check_rm('ionex*') # spinifex iono data
### DONE

# unpack tar files if present
for tarfile in glob.glob(data_dir + '/*tar'):
    if not os.path.exists(tarfile.replace('.tar','')):
        s.add(f'tar xf {tarfile} --one-top-level={data_dir}', log='tar.log', commandType='general')
if len(s.action_list) > 0:
    logger.info('Untar files...')
    s.run(check=True, maxProcs=5)

MSs = lib_ms.AllMSs(glob.glob(data_dir + '/*MS'), s, check_flags=False)

### This part is done so that missing subbands get concatenated correctly.
for i, msg in enumerate(np.array_split(sorted(glob.glob(data_dir + '*MS')), 1)):
    if fillmissingedges:
        min_nu = pt.table(MSs.getListStr()[0], ack=False).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MIN']
        max_nu = pt.table(MSs.getListStr()[0], ack=False).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MAX']
    else:
        min_nu = min(MSs.getFreqs()) / 1e6
        max_nu = max(MSs.getFreqs()) / 1e6
    num_init = lib_util.lofar_nu2num(min_nu) + 1  # +1 because FREQ_MIN/MAX somewhat have the lowest edge of the SB freq
    num_fin = lib_util.lofar_nu2num(max_nu) + 1
    prefix = re.sub('SB[0-9]*.MS', '', msg[0])
    msg = []
    for j in range(num_init, num_fin + 1):
        msg.append(prefix + 'SB%03i.MS' % j)

#msg = sorted(glob.glob(data_dir+'/*MS')) # does not fill gaps

calname = MSs.getListObj()[0].getNameField(checkCalName=True)
nchan = MSs.mssListObj[0].getNchan()
tint = MSs.mssListObj[0].getTimeInt()
freqres = MSs.mssListObj[0].getChanband()
uvlmin = 0.0

if skymodel == '':  # default case
    if MSs.hasIS:
        logger.warning('Sub-arcsecond models for LBA only partially tested! '
                       '3C196 is usable but not optimal and on the wrong scale. 3C380 is usable and scales between 40 and 70 MHz.'
                       'For 3C48 and 3C147, using S&H point sources which may give reasonable results.')
        if calname.lower() in ['3c196', '3c380']:
            skymodel = os.path.dirname(__file__) + '/../models/calib-highres.skymodel'
        elif calname.lower() in ['3c48', '3c147']:
            skymodel = os.path.dirname(__file__) + '/../models/calib-simple.skymodel'
    else:
        skymodel = os.path.dirname(__file__) + '/../models/calib-simple.skymodel'

logger.info(f"Initial time res: {tint:.1f}, nchan: {nchan}")

with w.if_todo('concat_all'):
    freqstep = 1  # keep all channels
    timestep = round(4 / tint)  # brings down to 4s
    # concat all SBs
    # SB.MS:DATA -> concat_all.MS:DATA
    logger.info('Concatenating data all (avg time %i)...' % timestep)
    lib_util.check_rm('concat_all.MS')
    s.add('DP3 ' + parset_dir + '/DP3-concat.parset msin="[' + ','.join(msg) + ']" msout=concat_all.MS \
              avg.freqstep=' + str(freqstep) + ' avg.timestep=' + str(timestep),
          log='concat.log', commandType='DP3')
    s.run(check=True)
### DONE

MSs_concat_all = lib_ms.AllMSs(['concat_all.MS'], s, check_flags=False)
MSs_concat_all.print_HAcov()
# for averaged data data set - the only purpose is to speed up the PA and FR solves
# small_freqres = 781248 # four SB - 390624 # two SB - 195312 # one SB
small_freqstep = round(781248/freqres)
small_timestep = 16 # to 64 s


######################################################
# rescale data to expected theoretical bandpass
if min(MSs_concat_all.getFreqs()) < 100e6:
    with w.if_todo('scale_bp'):
        # check weights
        MSs_concat_all.run('reweight.py $pathMS -v -p -a %s' % (MSs_concat_all.getListObj()[0].getAntennas()[0]),
                    log='$nameMS_weights.log', commandType='python')
        os.system('mkdir plots-weights; mv *png plots-weights/prebptheo.png')

        logger.info("Scale data to expected bandpass...")
        # Solve concat_all.MS:DATA
        # dummy call to create template
        MSs_concat_all.run(
            f'DP3 {parset_dir}/DP3-sol.parset msin=concat_all.MS msin.datacolumn=DATA sol.h5parm=cal-bp-theo.h5 sol.solint={round(20 / tint)} sol.nchan=1 \
                sol.maxiter=0 sol.mode=diagonalamplitude sol.modeldatacolumns="[DATA]" sol.datause=dual',
                log='$nameMS_bptemplate.log', commandType='DP3')
        lib_h5.create_h5bandpass(MSs_concat_all.getListObj()[0], 'cal-bp-theo.h5')
        # Correct avg BP concat_all:DATA -> DATA
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=DATA '
                           f'cor.parmdb=cal-bp-theo.h5 cor.correction=amplitude000 cor.updateweights=True',
                           log="$nameMS_bpscale.log", commandType="DP3")

        # check weights
        MSs_concat_all.run('reweight.py $pathMS -v -p -a %s' % (MSs_concat_all.getListObj()[0].getAntennas()[0]),
                    log='$nameMS_weights.log', commandType='python')
        os.system('mv *png plots-weights/postbptheo.png')

# flag bad stations, flags will propagate
with w.if_todo('flag'):
    logger.info("Flagging...")
    # MSs_concat_all.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS ant.baseline=\"{bl2flag}\" \
    #         aoflagger.strategy={parset_dir}/LBAdefaultwideband.lua steps=[ant]',
    #         log='$nameMS_flag.log', commandType='DP3')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS ant.baseline=\"{bl2flag}\" \
            aoflagger.strategy={parset_dir}/LBAdefaultwideband.lua',
            log='$nameMS_flag.log', commandType='DP3')
#
#     # extend flags on bad time/freq slots
#     logger.info('Remove bad time/freq stamps...')
#     MSs_concat_all.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
### DONE

# Predict cal concat_all.MS:MODEL_DATA
with w.if_todo('predict_all'):
    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    MSs_concat_all.run(
        f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}",
        log="$nameMS_pre.log", commandType="DP3")

###################################################

# 1: PRE-calibration: remove the fast-wrapping scalarphase (clock+disp. delay + 3rd order).
# This is necessary to increase the solution intervals/channels for the PA rot+diag step, that otherwise becomes
# too noisy for international stations. Do this on FR-corrected data to reduce glitches from FR+PA interplay.
with w.if_todo('pre_iono'):
    # Correct beam concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(
        f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn=DATA corrbeam.updateweights=True',
        log='$nameMS_beam.log', commandType="DP3")

    # check weights
    MSs_concat_all.run(
        'reweight.py $pathMS -v -p -a %s' % (MSs_concat_all.getListObj()[0].getAntennas()[0]),
        log='$nameMS_weights.log', commandType='python')
    os.system('mv *png plots-weights/postbeam.png')

    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # Solve concat_all.MS:SMOOTHED_DATA (only solve on CS-CS BL)
    # not more smoothing since in rare case some CS have strong delay!
    # solint = 8 * 4s = 32s
    logger.info('Calibrating IONO...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        sol.h5parm=$pathMS/preiono.h5 sol.mode=scalarphase sol.datause=single sol.solint=8 sol.nchan=2 \
                        sol.smoothnessconstraint=5.0e6 sol.smoothness_kernel_truncation=False sol.solutions_per_direction=[8] \
                        sol.antenna_averaging_factors=[CS*:8,RS*:1] sol.antenna_smoothness_factors=[CS*:1,RS*:0.05] \
                        sol.uvlambdamin={uvlmin}', log='$nameMS-preIONO_sol.log',
                       commandType="DP3")

    lib_util.run_losoto(s, 'preiono', [ms + '/preiono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset'])


########################################################

# 2: find PA
with w.if_todo('cal_pa'):
    # Correct pre-iono concat_all: DATA -> CORRECTED_DATA
    logger.info('Iono correction (preliminary)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono.h5 \
                msin.datacolumn=DATA cor.correction=phase000', log='$nameMS_cor-preIONO.log', commandType="DP3")
    # Smooth data concat_all:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # predict the element-corrupted model
    logger.info(f'Add beam-corrupted model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA_BEAMCOR...')
    MSs_concat_all.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_BEAMCOR \
                       pre.sourcedb={skymodel} pre.sources={calname} pre.beam_interval=120 \
                       pre.usebeammodel=True, pre.beammode=element", log="$nameMS_pre.log", commandType="DP3")

    # # Solve concat_all.MS:SMOOTHED_DATA (only solve)
    if MSs_concat_all.isLBA:
        logger.info('Calibrating PA...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/pa.h5 sol.modeldatacolumns=[MODEL_DATA_BEAMCOR] \
                   sol.mode=rotation+diagonal sol.rotationdiagonalmode=diagonalphase sol.solver=directionsolve sol.datause=full  \
                   sol.solint={small_timestep} sol.uvlambdamin={uvlmin} sol.nchan={small_freqstep} sol.minvisratio=0.1', log='$nameMS_solPA.log', commandType="DP3")

        lib_util.run_losoto(s, 'pa', [ms+'/pa.h5' for ms in MSs_concat_all.getListStr()],
                            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset',  parset_dir+'/losoto-pa.parset'])
    else:
        # Solve concat_all.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating PA...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/pa.h5 sol.modeldatacolumns=[MODEL_DATA_BEAMCOR] \
                   sol.mode=diagonal sol.datause=full  \
                   sol.solint={small_timestep} sol.uvlambdamin={uvlmin} sol.nchan={small_freqstep} sol.minvisratio=0.1', log='$nameMS_solPA.log', commandType="DP3")

        lib_util.run_losoto(s, 'pa', [ms+'/pa.h5' for ms in MSs_concat_all.getListStr()],
                            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset',  parset_dir+'/losoto-pa.parset'])

### DONE
########################################################
# 3: find FR
with w.if_todo('cal_fr'):
    if MSs_concat_all.isHBA:
        # Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
        logger.info('Polalign correction...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                       cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
        # Correct beam concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Beam correction...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                           log='$nameMS_beam.log', commandType="DP3")
        logger.info('Iono correction (preliminary)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono.h5 \
                        cor.correction=phase000', log='$nameMS_cor-preIONO.log', commandType="DP3")

        # Smooth data concat_all:CORRECTED_DATA -> SMOOTHED_DATA
        MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')
        # Probably we do not need BLsmoothing since we have long time intervals and smoothnessconstraint?
        logger.info('Converting to circular (DATA -> CIRC_PHASEDIFF_DATA)...')
        # lin2circ conversion TCxx.MS:DATA -> CIRC_PHASEDIFF_DATA # use no dysco here!
        MSs_concat_all.run(
            f'DP3 {parset_dir}/DP3-lin2circ.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CIRC_PHASEDIFF_DATA',
            log='$nameMS_lin2circ.log', commandType="DP3")

        # Get circular phase diff TC.MS: CIRC_PHASEDIFF_DATA -> CIRC_PHASEDIFF_DATA
        logger.info('Get circular phase difference...')
        MSs_concat_all.run('taql "UPDATE $pathMS SET\
             CIRC_PHASEDIFF_DATA[,0]=0.5*EXP(1.0i*(PHASE(CIRC_PHASEDIFF_DATA[,0])-PHASE(CIRC_PHASEDIFF_DATA[,3]))), \
             CIRC_PHASEDIFF_DATA[,3]=CIRC_PHASEDIFF_DATA[,0], \
             CIRC_PHASEDIFF_DATA[,1]=0+0i, \
             CIRC_PHASEDIFF_DATA[,2]=0+0i"', log='$nameMS_taql.log', commandType='general')

        logger.info('Creating MODEL_DATA_FR...')  # take from MODEL_DATA but overwrite
        MSs_concat_all.addcol('MODEL_DATA_FR', 'DATA', usedysco=False)
        MSs_concat_all.run('taql "UPDATE $pathMS SET MODEL_DATA_FR[,0]=0.5+0i, MODEL_DATA_FR[,1]=0.0+0i, MODEL_DATA_FR[,2]=0.0+0i, \
             MODEL_DATA_FR[,3]=0.5+0i"', log='$nameMS_taql.log', commandType='general')

        # Solve TC.MS: CIRC_PHASEDIFF_DATA against MODEL_DATA_FR (only solve - solint=2m nchan=0 as it has the smoothnessconstrain)
        logger.info('Solving circ phase difference ...')

        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5',
                log='$nameMS_solFR.log', commandType="DP3")
        lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_concat_all.getListStr()], [parset_dir + '/losoto-fr-circ.parset'])
    else:
        # Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
        logger.info('Polalign correction...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                       cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
        # Correct beam concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Beam correction...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                           log='$nameMS_beam.log', commandType="DP3")
        logger.info('Iono correction (preliminary)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono.h5 \
                        cor.correction=phase000', log='$nameMS_cor-preIONO.log', commandType="DP3")

        # Smooth data concat_all:CORRECTED_DATA -> SMOOTHED_DATA
        MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

        # Solve concat_all.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating FR...')
        # We solve for rot+diag or rot+scalar here and not just rot since we can have phase offsets from the preliminary iono!!
        # HE: do not use smoothnessconstraint, gives quite bad results here, at least at 1-2 MHz kernel and above!
        MSs_concat_all.run(
            f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 \
            sol.mode=rotation+diagonal sol.rotationdiagonalmode=scalarphase sol.minvisratio=0.1 sol.datause=full \
            sol.solint={small_timestep} sol.uvlambdamin={uvlmin} sol.nchan={int(small_freqstep / 2)}',
            log='$nameMS_solFR.log', commandType="DP3")

        lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_concat_all.getListStr()],
                            [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-plot-rot.parset',
                             parset_dir + '/losoto-fr.parset'])

### DONE
#################################################

# 4: calibrate iono + clock
with w.if_todo('cal_iono'):
    # Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                       log='$nameMS_beam.log', commandType="DP3")
    # Correct FR concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-fr.h5 \
                   cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    logger.info('Calibrating IONO...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.datause=single sol.solint=8 sol.nchan=2 \
                        sol.smoothnessconstraint=5.0e6 sol.smoothness_kernel_truncation=False sol.solutions_per_direction=[8] \
                        sol.antenna_averaging_factors=[CS*:8,RS*:1] sol.antenna_smoothness_factors=[CS*:1,RS*:0.05] \
                        sol.uvlambdamin={uvlmin}', log='$nameMS-preIONO_sol.log',
                       commandType="DP3")
    lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset'])
### DONE
######################################################

# 5: find BP
with w.if_todo('cal_bp'):
    ## Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                           log='$nameMS_beam.log', commandType="DP3")
    # FR corruption concat_all.MS:MODEL_DATA -> MODEL_DATA_FRCOR
    logger.info('Faraday rotation corruption (MODEL_DATA - > MODEL_DATA_FRCOR)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA_FRCOR \
                        cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000 cor.invert=False',
                       log='$nameMS_corFR.log', commandType="DP3")
    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve) against FR-corrupted MODEL_DATA
    # do not use datause=dual since the cross-hands are NOT trivial (FR-corrupted)
    logger.info('Calibrating BP...')
    timestep = round(20 / tint)  # brings down to 20s
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/bp-sub.h5 sol.mode=diagonal sol.datause=full \
                        sol.modeldatacolumns=[MODEL_DATA_FRCOR]  sol.uvlambdamin={uvlmin} sol.minvisratio=0.1 sol.solint={str(timestep)} sol.nchan=1',
                       log='$nameMS_solBP.log', commandType="DP3")

    flag_parset = '/losoto-flag-lessaggressive.parset'
    lib_util.run_losoto(s, 'bp-sub', [ms + '/bp-sub.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + flag_parset])

    # merge the solution with the bandpass before losoto
    s.add('h5_merger.py --h5_out cal-bp.h5 --h5_tables cal-bp-sub.h5 cal-bp-theo.h5 --propagate_flags'
            , log='h5_merger.log', commandType='python')
    s.run(check=True)

    # Custom losoto script for decameter
    lib_util.run_losoto(s, 'cal-bp.h5', 'cal-bp.h5', [parset_dir + '/losoto-bp-decameter.parset'], plots_dir='plots-bp')
### DONE

with w.if_todo('cal_final'):
    # Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                       log='$nameMS_beam.log', commandType="DP3")
    # Correct amp BP concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('amp correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-bp-sub.h5 \
        cor.correction=amplitude000 cor.updateweights=True', log='$nameMS_corBP.log', commandType="DP3")
    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    # FR correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA_FR
    logger.info('FR correction (for imaging)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msout.datacolumn=CORRECTED_DATA_FR \
                       cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
    ###


with w.if_todo('cal_imaging'):
    if not os.path.exists('img'):
        os.makedirs('img')

    logger.info('Cleaning...')
    lib_util.run_wsclean(s, 'wsclean.log', MSs_concat_all.getStrWsclean(), name='img/cal', size=2048, data_column='CORRECTED_DATA_FR',
                         scale='8arcsec', pol='I', auto_mask=5, parallel_gridding=6,
                         weight='uniform', niter=100000, mgain=0.4, no_update_model_required='', minuv_l=30,
                         baseline_averaging='', auto_threshold=2, join_channels='', fit_spectral_pol=4,
                         channels_out=MSs_concat_all.getChout(4e6))
with w.if_todo('cal_imaging_lres'):

    logger.info('Cleaning...')
    lib_util.run_wsclean(s, 'wsclean.log', MSs_concat_all.getStrWsclean(), name='img/cal-lres', size=2048, data_column='CORRECTED_DATA_FR',
                         scale='20arcsec', pol='I', auto_mask=5, parallel_gridding=6,
                         weight='briggs -1.0', maxuvw_m=3000,  niter=100000, mgain=0.4, no_update_model_required='',
                         baseline_averaging='', auto_threshold=2, join_channels='', fit_spectral_pol=4,
                         channels_out=MSs_concat_all.getChout(4e6))
