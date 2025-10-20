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

def debug_imaging(MSs, suffix, column='CORRECTED_DATA'):
    """
    Make an image of the calibrator
    Parameters
    ----------
    MSs
    suffix: name of the image
    """
    if not os.path.exists('img'):
        os.makedirs('img')

    imgsizepix = 2048 if MSs.hasIS else 1024
    scale = '0.3arcsec' if MSs.hasIS else '3arcsec'

    logger.info('Cleaning...')
    imagename = f'img/cal-{suffix}'
    lib_util.run_wsclean(s, 'wsclean.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, data_column=column,
                         scale=scale, pol='I,V', auto_mask=5, # local_rms='', local_rms_method='rms-with-min',
                         weight='briggs -0.3', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.6,
                         baseline_averaging='', auto_threshold=2, join_channels='', fit_spectral_pol=3,
                         channels_out=MSs.getChout(4e6), use_shm=use_shm)

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
for i, msg in enumerate(np.array_split(sorted(glob.glob(data_dir + '/*MS')), 1)):
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
small_timestep = 8 # to 32 s

# relax in case of IS or decamenter
if MSs_concat_all.hasIS:
    small_freqstep //= 2
#    small_timestep //= 2 # ~time variation not that different between RS and IS
if min(MSs_concat_all.getFreqs()) < 20.e6:
    small_freqstep //= 4
    small_timestep //= 2
elif min(MSs_concat_all.getFreqs()) < 40.e6:
    small_freqstep //= 2

uvlambdamin = 50 if min(MSs_concat_all.getFreqs()) < 30e6 else 100 # for Decameter we don't have any data otherwise...

######################################################
# rescale data to expected theoretical bandpass
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
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS ant.baseline=\"{bl2flag}\" \
            aoflagger.strategy={parset_dir}/LBAdefaultwideband.lua',
            log='$nameMS_flag.log', commandType='DP3')
    
    # extend flags on bad time/freq slots
    logger.info('Remove bad time/freq stamps...')
    MSs_concat_all.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
### DONE

# Predict cal concat_all.MS:MODEL_DATA
with w.if_todo('predict_all'):
    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    MSs_concat_all.run(
        f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}",
        log="$nameMS_pre.log", commandType="DP3")

if use_GNSS:
    with w.if_todo('get_gps_tec_rm'):
        lib_util.check_rm('cal-gps*.h5')
        # Get tec h5 parm from GPS data using spinifex (https://git.astron.nl/RD/spinifex).
        logger.info('Get RM from GPS data (spinifex)...')
        MSs_concat_all.run('spinifex get_rm_h5parm_from_ms $pathMS -o cal-gps-rm.h5',
                           log='spinifex_gps_rm.log', commandType='general')
        lib_util.run_losoto(s, 'cal-gps-rm.h5', ['cal-gps-rm.h5'], [parset_dir + '/losoto-plot-rm.parset'], plots_dir='plots-gps-rm')
        logger.info('Get TEC from GPS data (spinifex)...')
        MSs_concat_all.run('spinifex get_tec_h5parm_from_ms $pathMS -o cal-gps-tec.h5',
                           log='spinifex_gps_tec.log', commandType='general')
        # smooth gps TEC. (fitting works better on smoothed data)
        s.add("smooth_gps_tec.py cal-gps-tec.h5 tec", log='smooth_gps_tec.log', commandType='python')
        s.run()    
        
        lib_util.run_losoto(s, 'cal-gps-tec.h5', ['cal-gps-tec.h5'], [parset_dir + '/losoto-plot-tec.parset'], plots_dir='plots-gps-tec')

# if develop:
#     # Smooth data concat_all-all DATA -> SMOOTHED_DATA (BL-based smoothing)
#     MSs_concat_all.run_Blsmooth(nofreq=True, logstr='smooth3')
#
#     # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#     logger.info('Calibrating test...')
#     MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
#             sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
#             sol.solint=1 sol.nchan=1', \
#             log='$nameMS_solTEST.log', commandType="DP3")
#     lib_util.run_losoto(s, 'test-init', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
#                         [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

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

    if use_GNSS:
        # Preliminary rm correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('pre-correcion RM from GPS...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-rm.h5 \
                    cor.correction=rotationmeasure000', log='$nameMS_cor-gps-rm.log', commandType="DP3")
        # Preliminary tec correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('pre-correcion TEC from GPS...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-tec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-gps-tec.log', commandType="DP3")

    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # Solve concat_all.MS:SMOOTHED_DATA (only solve on CS-CS BL)
    # not more smoothing since in rare case some CS have strong delay!
    # solint = 8 * 4s = 32s
    logger.info('Calibrating IONO (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        sol.h5parm=$pathMS/preiono.h5 sol.mode=scalarphase sol.datause=single sol.solint=8 sol.nchan=1 msin.baseline="CS*&CS*" \
                        sol.smoothnessconstraint=0.5e6 sol.uvlambdamin={uvlambdamin}', log='$nameMS-preIONO_sol.log',
                       commandType="DP3")

    lib_util.run_losoto(s, 'preiono-cs', [ms + '/preiono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-ionoCS.parset'])

    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA (unit correction for RS)
    logger.info('Iono correction (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono-cs.h5 \
                cor.correction=phase000', log='$nameMS-preIONO_cor.log', commandType="DP3")

    # Phasing up the cose stations concat_all.MS:SMOOTHED_DATA -> concat_all-phaseup-preIONO.MS:DATA
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseup-preIONO.MS')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-phaseup.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        msout=concat_all-phaseup-preIONO.MS', log='$nameMS_phaseup.log', commandType="DP3")

    MSs_concat_phaseupIONO = lib_ms.AllMSs(['concat_all-phaseup-preIONO.MS'], s, check_flags=False)

    # Predict cal concat_all.MS:MODEL_DATA
    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    MSs_concat_phaseupIONO.run(
        f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}",
        log="$nameMS_pre.log", commandType="DP3")
    # TODO RS and IS need very different smoothing - could do an extra iteration for the IS?
    # Equation for allowed smoothing, assuming ONLY TEC and kernel is one sixth of the bandwidth it takes for one wrap at 54 MHz
    # kernelsize_smoothnessconstraint [MHz] = 0.3 / dTEC [TECU]
    # For RS: Expect up to 0.5 TECU, kernel should be ~0.6 MHz (smaller since we also have clock)
    # FOR IS: Expect up to 5 TECU, kernel should be ~0.1 MHz
    # Smooth data concat_all.MS:DATA -> SMOOTHED_DATA
    MSs_concat_phaseupIONO.run_Blsmooth(incol='DATA', logstr='smooth')
    # Solve concat_all-phaseupIONO.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating IONO (distant stations)...')
    smoothnessconstraint = '0.1e6' if MSs_concat_all.hasIS else '0.5e6'
    MSs_concat_phaseupIONO.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                           sol.h5parm=$pathMS/preiono.h5 sol.mode=scalarphase  sol.datause=single \
                           sol.solint=1 sol.nchan=1 sol.smoothnessconstraint={smoothnessconstraint} sol.smoothnessreffrequency=54e6', \
                           log='$nameMS_sol.log', commandType="DP3")

    lib_util.run_losoto(s, 'preiono', [ms + '/preiono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                        [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset'])
    if use_GNSS:
        lib_util.check_rm('cal-dtec.h5')
        logger.info('fit residual dTEC...')
        s.add("dtec_finder.py --gps_corrected cal-preiono.h5", log='dtec_finder.log', commandType='python')
        s.run(check=True)
        lib_util.run_losoto(s, 'cal-dtec.h5', ['cal-dtec.h5'], [parset_dir + '/losoto-plot-tec.parset'], plots_dir='plots-dtec-finder')

### DONE
########################################################

# 2: find PA
with w.if_todo('cal_pa'):
    if use_GNSS:
        # Correct gps-tec concat_all:DATA -> CORRECTED_DATA
        logger.info('TEC correction (GPS)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-tec.h5 msin.datacolumn=DATA\
                    cor.correction=tec000', log='$nameMS_cor-gps-tec.log', commandType="DP3")
    else:
        # Correct pre-iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Iono correction (preliminary)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono-cs.h5 \
                    cor.correction=phase000', log='$nameMS_cor-preIONO.log', commandType="DP3")
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono.h5 \
                    cor.correction=phase000', log='$nameMS_cor-preIONO.log', commandType="DP3")
    # Smooth data concat_all:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    if MSs_concat_all.hasIS:
        # For IS data speedup: Concat+avg SMOOTHED_DATA -> concat_pa.MS:DATA
        # average a lot to speed up and increase S/N -> fast scalarphase is gone, so this should be fine.
        # the averaging here needs to be sufficient to catch the FR
        logger.info('Create small concatenated MS (averaging)...')
        lib_util.check_rm('concat_pa.MS')
        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin=concat_all.MS msout=concat_pa.MS msin.datacolumn=SMOOTHED_DATA \
                   avg.timestep={small_timestep} avg.freqstep={small_freqstep}', log='concat_pa.log', commandType='DP3')
        s.run(check=True)
        MSs_pa = lib_ms.AllMSs(['concat_pa.MS'], s, check_flags=False)

        # predict the element-corrupted model
        logger.info(f'Add beam-corrupted model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
        MSs_pa.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname} \
                  pre.usebeammodel=True, pre.beammode=element pre.beam_interval=120", log="$nameMS_pre.log", commandType="DP3")
        if use_GNSS:
            # FR corruption concat_pa.MS:MODEL_DATA -> MODEL_DATA
            logger.info('Faraday rotation corruption (MODEL_DATA - > MODEL_DATA)...')
            MSs_pa.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                                cor.parmdb=cal-gps-rm.h5 cor.correction=rotationmeasure000 cor.invert=False',
                       log='$nameMS_corGPSFR.log', commandType="DP3")
        # Beam corruption concat_pa.MS:MODEL_DATA -> MODEL_DATA
        logger.info(f'Beam model corruption (MODEL_DATA - > MODEL_DATA)...')
        phase_center = MSs_pa.getListObj()[0].getPhaseCentre()
        MSs_pa.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA setbeam.beammode=Element corrbeam.updateweights=False corrbeam.invert=False',
            log='$nameMS_beam.log', commandType="DP3")
        # HE: sol.rotationdiagonalmode diagonalphase seemes to give more stable results and surpresses the ~60 MHz bump weirdness
        # HE: do not use smoothnessconstraint, gives quite bad results here, at least at 1-2 MHz kernel and above!
        # Solve concat_pa.MS:DATA (only solve)
        logger.info('Calibrating PA...')
        MSs_pa.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/pa.h5 \
               sol.mode=rotation+diagonal sol.rotationdiagonalmode=diagonalphase sol.datause=full \
               sol.solint=1 sol.nchan=1', log='$nameMS_solPA.log', commandType="DP3")

        lib_util.run_losoto(s, 'pa', [ms+'/pa.h5' for ms in MSs_pa.getListStr()],
                [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset',  parset_dir+'/losoto-pa.parset'])

    else:
        # predict the element-corrupted model
        logger.info(f'Add beam-corrupted model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA_BEAMCOR...')
        MSs_concat_all.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_BEAMCOR \
                           pre.sourcedb={skymodel} pre.sources={calname} pre.beam_interval=120 \
                           pre.usebeammodel=True, pre.beammode=element", log="$nameMS_pre.log", commandType="DP3")
        if use_GNSS:
            # FR corruption concat_concat_all.MS:MODEL_DATA_BEAMCOR -> MODEL_DATA_BEAMCOR
            MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_BEAMCOR msout.datacolumn=MODEL_DATA_BEAMCOR \
                                cor.parmdb=cal-gps-rm.h5 cor.correction=rotationmeasure000 cor.invert=False',
                            log='$nameMS_corGPSFR.log', commandType="DP3")

        # Solve concat_all.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating PA...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/pa.h5 sol.modeldatacolumns=[MODEL_DATA_BEAMCOR] \
                   sol.mode=rotation+diagonal sol.rotationdiagonalmode=diagonalphase sol.datause=full \
                   sol.solint={small_timestep} sol.nchan={small_freqstep}', log='$nameMS_solPA.log', commandType="DP3")

        lib_util.run_losoto(s, 'pa', [ms+'/pa.h5' for ms in MSs_concat_all.getListStr()],
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset',  parset_dir+'/losoto-pa.parset'])
        
### DONE
########################################################
# 3: find FR
with w.if_todo('cal_fr'):
    # Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                       log='$nameMS_beam.log', commandType="DP3")
    if use_GNSS:
        # Preliminary RM correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('FR pre-correction (GPS)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-rm.h5 \
                    cor.correction=rotationmeasure000', log='$nameMS_cor-gps-rm.log', commandType="DP3")
        # Correct gps-tec concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('TEC pre-correction (GPS)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-tec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-gps-tec.log', commandType="DP3")
        # Correct TEC concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('dTEC correction (fitted)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-dtec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-dtec.log', commandType="DP3")
    else:
        logger.info('Iono correction (preliminary)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono-cs.h5 \
                    cor.correction=phase000', log='$nameMS_cor-preIONO.log', commandType="DP3")
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
        sol.mode=rotation+diagonal sol.rotationdiagonalmode=scalarphase sol.datause=full \
        sol.solint={small_timestep} sol.nchan={int(small_freqstep / 2)}', 
        log='$nameMS_solFR.log', commandType="DP3")

    if (min(MSs_concat_all.getFreqs()) < 30.e6):
        lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_concat_all.getListStr()],
                            [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-plot-rot.parset',
                             parset_dir + '/losoto-fr-low.parset'])
    elif MSs_concat_all.hasIS:
        lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_concat_all.getListStr()],
                            [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-plot-rot.parset',
                             parset_dir + '/losoto-fr-IS.parset'])
        # workaround to remove all flags from cal-fr.h5
        #logger.info('unflag cal-fr.h5...')
        #s.add("h5_remove_flags.py cal-fr.h5 rotationmeasure", log='h5_remove_flag.log', commandType='python')
        #s.run()
        #os.system('mv cal-fr.h5 cal-fr-original.h5')
        #os.system('mv cal-fr-noflag.h5 cal-fr.h5')
    else:
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
    if use_GNSS:
        # Correct gps-tec concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('TEC correction (GPS)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-tec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-gps-tec.log', commandType="DP3")
        # Correct TEC concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('dTEC correction (fitted)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-dtec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-dtec.log', commandType="DP3")
        # Correct FR concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation pre-correction (GPS)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-gps-rm.h5 \
                      cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
    # Correct FR concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-fr.h5 \
                   cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # Solve concat_all.MS:DATA (only solve on CS-CS BL)
    # not more smoothing since in rare case some CS have strong delay!
    logger.info('Calibrating IONO (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.datause=single sol.solint=8 sol.nchan=1 msin.baseline="CS*&CS*" \
                        sol.smoothnessconstraint=0.5e6 sol.uvlambdamin={uvlambdamin}', log='$nameMS_solIONO.log',
                       commandType="DP3")

    lib_util.run_losoto(s, 'iono-cs', [ms + '/iono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-ionoCS.parset'])

    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA (unit correction for others)
    logger.info('Iono correction (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")

    # Phasing up the cose stations CORRECTED_DATA -> concat_all-phaseup-IONO.MS:DATA
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseup-IONO.MS')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-phaseup.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        msout=concat_all-phaseup-IONO.MS', log='$nameMS_phaseup.log', commandType="DP3")

    MSs_concat_phaseupIONO = lib_ms.AllMSs(['concat_all-phaseup-IONO.MS'], s, check_flags=False)

    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    os.system('cp -r %s %s' % (skymodel, MSs_concat_phaseupIONO.getListObj()[0].pathMS))
    MSs_concat_phaseupIONO.run(
        f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}',
        log="$nameMS_pre.log", commandType="DP3")

    # Smooth data concat_all-phaseup-IONO.MS:DATA -> SMOOTHED_DATA
    MSs_concat_phaseupIONO.run_Blsmooth(incol='DATA', nofreq=True, logstr='smooth')
    # Solve concat_all-phaseup-IONO.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating IONO (distant stations)...')
    smoothnessconstraint = '0.1e6' if MSs_concat_all.hasIS else '0.5e6'
    MSs_concat_phaseupIONO.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS \
                           sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.datause=single \
                           sol.solint=1 sol.nchan=1 sol.smoothnessconstraint={smoothnessconstraint} sol.smoothnessreffrequency=54e6', \
                           log='$nameMS_sol.log', commandType="DP3")

    if (min(MSs_concat_phaseupIONO.getFreqs()) < 35.e6):
        lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset',
                             parset_dir + '/losoto-iono3rd.parset'])
    else:
        lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset',
                             parset_dir + '/losoto-iono.parset'])
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
    if use_GNSS:
        # Correct gps-TEC concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('TEC correction (GPS)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-gps-tec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-gps-tec.log', commandType="DP3")
        # Correct TEC concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('dTEC correction (fitted)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-dtec.h5 \
                    cor.correction=tec000', log='$nameMS_cor-dtec.log', commandType="DP3")
        # FR - prepcorruption concat_all.MS:MODEL_DATA_FRCOR -> MODEL_DATA_FRCOR
        logger.info('Faraday rotation pre-corruption (GPS) (MODEL_DATA_FRCOR - > MODEL_DATA_FRCOR)...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_FRCOR msout.datacolumn=MODEL_DATA_FRCOR \
                            cor.parmdb=cal-gps-rm.h5 cor.correction=rotationmeasure000 cor.invert=False',
                           log='$nameMS_corFR.log', commandType="DP3")
    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve) against FR-corrupted MODEL_DATA
    # do not use datause=dual since the cross-hands are NOT trivial (FR-corrupted)
    logger.info('Calibrating BP...')
    timestep = round(20 / tint)  # brings down to 20s
    if min(MSs_concat_all.getFreqs()) < 20.e6:
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/bp-sub.h5 sol.mode=diagonal sol.datause=full \
                            sol.modeldatacolumns=[MODEL_DATA_FRCOR] sol.solint={str(timestep)} sol.nchan=1 sol.smoothnessconstraint=1e6 \
                            sol.smoothnessreffrequency=0. ',
                           log='$nameMS_solBP.log', commandType="DP3")
    else:
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/bp-sub.h5 sol.mode=diagonal sol.datause=full \
                            sol.modeldatacolumns=[MODEL_DATA_FRCOR] sol.solint={str(timestep)} sol.nchan=1',
                           log='$nameMS_solBP.log', commandType="DP3")

    flag_parset = '/losoto-flag-lessaggressive.parset' if less_aggressive_flag else '/losoto-flag.parset'
    lib_util.run_losoto(s, 'bp-sub', [ms + '/bp-sub.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + flag_parset])

    # merge the solution with the bandpass before losoto
    s.add('h5_merger.py --h5_out cal-bp.h5 --h5_tables cal-bp-sub.h5 cal-bp-theo.h5 --propagate_flags'
            , log='h5_merger.log', commandType='python')
    s.run(check=True)

    # Custom losoto script for decameter
    if min(MSs_concat_all.getFreqs()) < 20.e6:
        lib_util.run_losoto(s, 'cal-bp.h5', 'cal-bp.h5', [parset_dir + '/losoto-bp-decameter.parset'], plots_dir='plots-bp')
    else:
        lib_util.run_losoto(s, 'cal-bp.h5', 'cal-bp.h5', [parset_dir + '/losoto-bp.parset'], plots_dir='plots-bp')
### DONE

if develop:
    # Restore original DATA
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=DATA '
                       f'cor.parmdb=cal-bp-theo.h5 cor.correction=amplitude000 cor.invert=False',
                       log="$nameMS_bpscaleTEST.log", commandType="DP3")

    # Pol align correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 \
                   cor.correction=polalign', log='$nameMS_corPATEST.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pa', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False',
                           log='$nameMS_beamTEST.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeam', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('BP correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-bp.h5 \
                cor.correction=amplitudeSmooth cor.updateweights=False',
                log='$nameMS_corBPTEST.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambp', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 \
                   cor.correction=rotationmeasure000', log='$nameMS_corFR-TEST.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambpfr', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Correct iono CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                   cor.correction=phase000', log='$nameMS_corIONO_CS-TEST.log', commandType="DP3")
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                   cor.correction=phase000', log='$nameMS_corIONO-TEST.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambpfriono', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

if not develop:
    with w.if_todo('compressing_h5'):
        logger.info('Compressing caltables...')
        # os.system('cp cal-pa.h5 fullcal-pa.h5')
        # os.system('cp cal-fr.h5 fullcal-fr.h5') # no need to keep orig
        # os.system('cp cal-bp.h5 fullcal-bp.h5')
        # os.system('cp cal-iono.h5 fullcal-iono.h5')
        s.add('losoto -d sol000/phaseOrig000 cal-pa.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/phase000 cal-pa.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/rotation000 cal-pa.h5', log='losoto-compress.log', commandType="python")

        s.add('losoto -d sol000/phase000 cal-fr.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/rotation000 cal-fr.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/rotationResid000 cal-fr.h5', log='losoto-compress.log', commandType="python")

        s.add('losoto -d sol000/amplitude000 cal-bp.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/phase000 cal-bp.h5', log='losoto-compress.log', commandType="python")
        #s.add('losoto -d sol000/amplitudeRes cal-bp.h5', log='losoto-compress.log', commandType="python") # keep for quality evaluation

        s.add('losoto -d sol000/phase_offset000 cal-iono-cs.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/phaseResid000 cal-iono-cs.h5', log='losoto-compress.log', commandType="python")

        s.add('losoto -d sol000/phase_offset000 cal-iono.h5', log='losoto-compress.log', commandType="python")
        s.add('losoto -d sol000/phaseResid000 cal-iono.h5', log='losoto-compress.log', commandType="python")

        s.run(maxProcs=1, check=True) # final check on losoto-compress.log

        os.system('h5repack cal-pa.h5 cal-pa-compressed.h5; mv cal-pa-compressed.h5 cal-pa.h5')
        os.system('h5repack cal-fr.h5 cal-fr-compressed.h5; mv cal-fr-compressed.h5 cal-fr.h5')
        os.system('h5repack cal-bp.h5 cal-bp-compressed.h5; mv cal-bp-compressed.h5 cal-bp.h5')
        os.system('h5repack cal-iono-cs.h5 cal-iono-cs-compressed.h5; mv cal-iono-cs-compressed.h5 cal-iono-cs.h5')
        os.system('h5repack cal-iono.h5 cal-iono-compressed.h5; mv cal-iono-compressed.h5 cal-iono.h5')
### DONE

# a debug image
if imaging:
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
        logger.info('BP correction...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-bp-sub.h5 \
            cor.correction=amplitude000 cor.updateweights=True', log='$nameMS_corBP.log', commandType="DP3")
        # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Iono correction...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                    cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                    cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")

        
        # check weights
        #MSs_concat_all.run('reweight.py $pathMS -v -p -a %s' % (MSs_concat_all.getListObj()[0].getAntennas()[0]),
        #        log='$nameMS_weights.log', commandType='python')
        #os.system('mv *png plots-weights/postbp.png')

        ### DEBUG
        # FR correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA_FR
        #logger.info('FR correction (for imaging)...')
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msout.datacolumn=CORRECTED_DATA_FR \
        #                    cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
        #debug_imaging(MSs_concat_all, 'afterbp', column='CORRECTED_DATA_FR')
        ###

        # Here we solve for fast scalaramp - this should capture scintillations. BLsmooth might smear the scintillations for the CS, but seems to give better results in image space
        # Solve concat_all.MS:CORRECTED_DATA (only solve) against FR-corrupted MODEL_DATA
        #logger.info('Calibrating amp scintillations...')
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/amp.h5 sol.mode=scalaramplitude sol.datause=full \
        #                    sol.modeldatacolumns=[MODEL_DATA_FRCOR] sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=5e6',
        #                   log='$nameMS_solBP.log', commandType="DP3")

        #lib_util.run_losoto(s, 'amp', [ms + '/amp.h5' for ms in MSs_concat_all.getListStr()],
        #                    [parset_dir + '/losoto-plot-scalaramp.parset'])

        # Fast amp correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Amp scintillation correction...')
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-amp.h5 \
        #    cor.correction=amplitude000  cor.updateweights=True', log='$nameMS_corBP.log', commandType="DP3")

        # check weights
        #MSs_concat_all.run('reweight.py $pathMS -v -p -a %s' % (MSs_concat_all.getListObj()[0].getAntennas()[0]),
        #        log='$nameMS_weights.log', commandType='python')
        #os.system('mv *png plots-weights/postfastamp.png')

        ### DEBUG
        # FR correction concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA_FR
        #logger.info('FR correction (for imaging)...')
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msout.datacolumn=CORRECTED_DATA_FR \
        #                    cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
        #debug_imaging(MSs_concat_all, 'afteramp', column='CORRECTED_DATA_FR')
        ###

        # Solve concat_all.MS:CORRECTED_DATA (only solve)
        #logger.info('Calibrating slow leakage...')
        #timestep = 20*round(120 / tint)  # brings down to 20min
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
        #                   sol.h5parm=$pathMS/fj.h5 sol.mode=fulljones sol.minvisratio=0 \
        #                   sol.modeldatacolumns=[MODEL_DATA_FRCOR] sol.solint={str(timestep)} sol.smoothnessconstraint=5e6',
        #                   log='$nameMS_solFJ.log', commandType="DP3")

        #lib_util.run_losoto(s, 'fj', [ms + '/fj.h5' for ms in MSs_concat_all.getListStr()],
        #                    [parset_dir + '/losoto-plot-fullj.parset'])
        
        # Correct amp BP concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Leakage correction...')
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fj.h5 \
        #                    cor.correction=fulljones cor.soltab=[amplitude000,phase000] cor.updateweights=False',
        #                    log='$nameMS_corBP.log', commandType="DP3")

        # Correct FR concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Faraday rotation correction...')
        #MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 \
        #               cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")

    with w.if_todo('cal_imaging'):
        debug_imaging(MSs_concat_all, 'final')
        logger.info('Source to peel: set SUBTRACTED_CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_FRCOR')
        MSs_concat_all.addcol('SUBTRACTED_CORRECTED_DATA', 'DATA', log='$nameMS_addcol.log')
        MSs_concat_all.run('taql "update $pathMS set SUBTRACTED_CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_FRCOR"', \
                        log='$nameMS_taql.log', commandType='general')
        debug_imaging(MSs_concat_all, 'subtractedcorr', column='SUBTRACTED_CORRECTED_DATA')
    ### DONE

    with w.if_todo('model_corruption'):
        MSs_concat_all.addcol('MODEL_DATA_CORRUPT', 'MODEL_DATA', log='$nameMS_addcol.log')
        # Correct FR concat_all.MS:MODEL_DATA -> MODEL_DATA_CORRUPT
        logger.info('Faraday rotation corruption...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA_CORRUPT\
                       cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")
        # Correct iono concat_all:MODEL_DATA_CORRUPT -> MODEL_DATA_CORRUPT
        logger.info('Iono corruption...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_CORRUPT msout.datacolumn=MODEL_DATA_CORRUPT \
                           cor.parmdb=cal-iono-cs.h5 cor.correction=phase000  cor.invert=False', log='$nameMS_corIONO_CS.log', commandType="DP3")
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_CORRUPT msout.datacolumn=MODEL_DATA_CORRUPT \
                           cor.parmdb=cal-iono.h5 cor.correction=phase000 cor.invert=False', log='$nameMS_corIONO.log', commandType="DP3")
        # Correct amp BP concat_all.MS:MODEL_DATA_CORRUPT -> MODEL_DATA_CORRUPT
        logger.info('BP corruption...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_CORRUPT msout.datacolumn=MODEL_DATA_CORRUPT \
                           cor.parmdb=cal-bp-sub.h5 cor.correction=amplitude000 cor.updateweights=False cor.invert=False', log='$nameMS_corBP.log', commandType="DP3")
        # Correct beam concat_all.MS:MODEL_DATA_CORRUPT -> MODEL_DATA_CORRUPT
        logger.info('Beam corruption...')
        MSs_concat_all.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS msin.datacolumn=MODEL_DATA_CORRUPT msout.datacolumn=MODEL_DATA_CORRUPT \
                        setbeam.beammode=element corrbeam.invert=False', log='$nameMS_beam.log', commandType='DP3')
        # Pol align correction concat_all.MS:MODEL_DATA_CORRUPT -> MODEL_DATA_CORRUPT
        logger.info('Polalign corruption...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_CORRUPT msout.datacolumn=MODEL_DATA_CORRUPT \
                       cor.parmdb=cal-pa.h5 cor.correction=polalign cor.invert=False', log='$nameMS_corPA.log', commandType="DP3")
        
        logger.info('Source to peel: set SUBTRACTED_DATA = DATA - MODEL_DATA_CORRUPT')
        MSs_concat_all.addcol('SUBTRACTED_DATA', 'DATA', log='$nameMS_addcol.log')
        MSs_concat_all.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA_CORRUPT"', \
                        log='$nameMS_taql.log', commandType='general')
        debug_imaging(MSs_concat_all, 'subtracted', column='SUBTRACTED_DATA')


if not develop:
    logger.info('Cleaning up...')
    lib_util.check_rm('cal-preiono.h5 cal-preiono-cs.h5 cal-bp-sub.h5 cal-bp-theo.h5')
    os.system('rm -r *MS')

w.alldone()
