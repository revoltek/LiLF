#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation .
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

# TODO: debugplots currently broken
# TODO reset minvisration to 0.3 once bug fixed

import sys, os, glob, re
import casacore.tables as pt
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_h5

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
sparse_sb = parset.getboolean('LOFAR_cal', 'sparse_sb') # change flagging to hande data that uses only alternating sb
bl2flag = parset.get('flag', 'stations')
debugplots = False

#############################################################

def debug_imaging(MSs, suffix):
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
    lib_util.run_wsclean(s, 'wsclean.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix,
                         scale=scale, auto_mask=5, # local_rms='', local_rms_method='rms-with-min',
                         weight='briggs -0.3', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.6,
                         baseline_averaging='', auto_threshold=2, join_channels='', fit_spectral_pol=5,
                         channels_out=MSs.getChout(4e6))

#############################################################

# unpack tar files if present
for tarfile in glob.glob(data_dir + '/*tar'):
    if not os.path.exists(tarfile.replace('.tar','')):
        s.add(f'tar xf {tarfile} --one-top-level={data_dir}', log='tar.log', commandType='general')
if len(s.action_list) > 0:
    logger.info('Untar files...')
    s.run(check=True, maxThreads=5)

MSs = lib_ms.AllMSs(glob.glob(data_dir + '/*MS'), s, check_flags=False)

### This part is done so that missing subbands get concatenated correctly.
for i, msg in enumerate(np.array_split(sorted(glob.glob(data_dir+'/*MS')), 1)):
    if fillmissingedges:
        min_nu = pt.table(MSs.getListStr()[0], ack=False).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MIN']
        max_nu = pt.table(MSs.getListStr()[0], ack=False).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MAX']
    else:
        min_nu = min(MSs.getFreqs())/1e6
        max_nu = max(MSs.getFreqs())/1e6
    num_init = lib_util.lofar_nu2num(min_nu) + 1  # +1 because FREQ_MIN/MAX somewhat have the lowest edge of the SB freq
    num_fin = lib_util.lofar_nu2num(max_nu) + 1
    prefix = re.sub('SB[0-9]*.MS', '', msg[0])
    msg = []
    for j in range(num_init, num_fin + 1):
        msg.append(prefix + 'SB%03i.MS' % j)

#msg = sorted(glob.glob(data_dir+'/*MS')) # does not fill gaps

if skymodel == '':  # default case
    if MSs.hasIS:
        skymodel = os.path.dirname(__file__) + '/../models/calib-highres.skymodel'
    else:
        skymodel = os.path.dirname(__file__) + '/../models/calib-simple.skymodel'

calname = MSs.getListObj()[0].getNameField()
nchan = MSs.mssListObj[0].getNchan()
tint = MSs.mssListObj[0].getTimeInt()
freqres = MSs.mssListObj[0].getChanband()

logger.info(f"Initial time res: {tint:.1f}, nchan: {nchan}")

with w.if_todo('concat_all'):
    freqstep = 1  # keep all channels
    timestep = int(np.rint(4 / tint))  # brings down to 4s
    # concat all SBs
    # SB.MS:DATA -> concat.MS:DATA
    logger.info('Concatenating data all (avg time %i)...' % timestep)
    lib_util.check_rm('concat_all.MS')
    s.add('DP3 ' + parset_dir + '/DP3-concat.parset msin="[' + ','.join(msg) + ']" msout=concat_all.MS \
              msin.baseline="*&" avg.freqstep=' + str(freqstep) + ' avg.timestep=' + str(timestep),
          log='concat.log', commandType='DP3')
    s.run(check=True)
### DONE

MSs_concat_all = lib_ms.AllMSs(['concat_all.MS'], s, check_flags=False)
MSs_concat_all.print_HAcov()
# for averaged data data set - the only purpose is to speed up the PA and FR solves
#small_freqres = 781248 # four SB -390624 # two SB - 195312 # one SB
small_freqstep = int(np.rint(781248/freqres))
small_timestep = 8 # to 32 s

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
# flag bad stations, flags will propagate
with w.if_todo('flag'):
    logger.info("Flagging...")
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS ant.baseline=\"{bl2flag}\"',
                       log="$nameMS_flag.log", commandType="DP3")
    # extend flags
    logger.info('Remove bad time/freq stamps...')
    MSs_concat_all.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
### DONE

# Predict cal
with w.if_todo('predict_all'):
    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    MSs_concat_all.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}",
                       log="$nameMS_pre.log", commandType="DP3")

if debugplots:
    # Smooth data concat_all-all DATA -> SMOOTHED_DATA (BL-based smoothing)
    MSs_concat_all.run_Blsmooth(nofreq=True, logstr='smooth3')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS\
            sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
            sol.solint=1 sol.nchan=1', \
            log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-init', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

###################################################

# 1: PRE-calibration: remove the fast-wrapping scalarphase (clock+disp. delay + 3rd order).
# This is necessary to increase the solution intervals/channels for the PA rot+diag step, that otherwise becomes
# too noisy for international stations. Do this on FR-corrected data to reduce glitches from FR+PA interplay.
with w.if_todo('pre_cal'):
    # Correct beam concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn=DATA corrbeam.updateweights=True',
                       log='$nameMS_beam.log', commandType="DP3")
    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # Solve concat_all.MS:DATA (only solve on CS-CS BL)
    # not more smoothing since in rare case some CS have strong delay!
    # solint = 15 * 4s = 60s
    # HE: we could also experiment with moving to rotation+scalar on a somewhat lower time resolution here
    logger.info('Calibrating IONO (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.minvisratio=0.1 \
                        sol.h5parm=$pathMS/preiono.h5 sol.mode=scalarphase sol.datause=single sol.solint=15 sol.nchan=1 msin.baseline="CS*&CS*" \
                        sol.smoothnessconstraint=0.5e6 sol.uvlambdamin={uvlambdamin}', log='$nameMS_solIONO_CS.log',
                       commandType="DP3")

    lib_util.run_losoto(s, 'preiono-cs', [ms + '/preiono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-ionoCS.parset'])

    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA (unit correction for RS)
    logger.info('Iono correction (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono-cs.h5 \
                cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")

    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # Phasing up the cose stations concat_all.MS:SMOOTHED_DATA -> concat_all-phaseupIONO.MS:DATA
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseupIONO.MS')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-phaseup.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        msout=concat_all-phaseupIONO.MS', log='$nameMS_phaseup.log', commandType="DP3")

    MSs_concat_phaseupIONO = lib_ms.AllMSs(['concat_all-phaseupIONO.MS'], s, check_flags=False)

    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    MSs_concat_phaseupIONO.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}",
                       log="$nameMS_pre.log", commandType="DP3")

    # Solve concat_all-phaseupIONO.MS:DATA (only solve)
    logger.info('Calibrating IONO (distant stations)...')
    MSs_concat_phaseupIONO.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=DATA \
                           sol.h5parm=$pathMS/preiono.h5 sol.mode=scalarphase sol.datause=single sol.minvisratio=0.1 \
                           sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=0.1e6 sol.smoothnessreffrequency=54e6', \
                           log='$nameMS_solIONO.log', commandType="DP3")

    if min(MSs_concat_all.getFreqs()) < 35.e6:
        lib_util.run_losoto(s, 'preiono', [ms + '/preiono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset',
                             parset_dir + '/losoto-iono3rd.parset'])
    else:
        lib_util.run_losoto(s, 'preiono', [ms + '/preiono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset',
                             parset_dir + '/losoto-iono.parset'])
### DONE
########################################################

# 2: find PA
with w.if_todo('cal_pa'):
    # Correct iono concat_all:DATA -> CORRECTED_DATA
    logger.info('Iono correction (preliminary)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono-cs.h5 msin.datacolumn=DATA \
                cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono.h5 \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    # Smooth data concat_all:DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    if MSs_concat_all.hasIS:
        # For IS data speedup: Concat+avg SMOOTHED_DATA -> concat_pa.MS:DATA
        # average a lot to speed up and increase S/N -> fast scalarphase is gone, so this should be fine.
        # the averaging here needs to be sufficient to catch the FR
        logger.info('Create small concatenated MS (averaging)...')
        lib_util.check_rm('concat_pa.MS')
        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin=concat_all.MS msout=concat_pa.MS msin.datacolumn=SMOOTHED_DATA \
                   avg.timestep={small_timestep} avg.freqstep={small_freqstep}', log='concat_pa_model.log', commandType='DP3')
        s.run(check=True)
        MSs_pa = lib_ms.AllMSs(['concat_pa.MS'], s, check_flags=False)

        # predict the element-corrupted model
        logger.info('Add beam-corrupted model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
        MSs_pa.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname} \
                  pre.usebeammodel=True, pre.beammode=element", log="$nameMS_pre.log", commandType="DP3")

        # HE: sol.rotationdiagonalmode diagonalphase seemes to give more stable results and surpresses the ~60 MHz bump weirdness
        # Solve concat_pa.MS:DATA (only solve)
        logger.info(f'Calibrating PA...')
        MSs_pa.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal sol.rotationdiagonalmode=diagonalphase \
               sol.solint=1 sol.nchan=1', log='$nameMS_solPA.log', commandType="DP3")

    else:
        # predict the element-corrupted model
        logger.info(f'Add beam-corrupted model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA_BEAMCOR...')
        MSs_concat_all.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_BEAMCOR \
                           pre.sourcedb={skymodel} pre.sources={calname} \
                           pre.usebeammodel=True, pre.beammode=element", log="$nameMS_pre.log", commandType="DP3")

        # Solve concat_pa.MS:DATA (only solve)
        logger.info('Calibrating PA...')
        MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.modeldatacolumns=[MODEL_DATA_BEAMCOR] \
                   sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal \
                   sol.solint={small_timestep} sol.nchan={small_freqstep}', log='$nameMS_solPA.log', commandType="DP3")

    lib_util.run_losoto(s, 'pa-noavg', [ms+'/pa.h5' for ms in MSs_concat_all.getListStr()],
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
    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction (preliminary)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono-cs.h5 \
                    cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-preiono.h5 \
                    cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    # Smooth data concat_all:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')

    # Concat+avg SMOOTHED_DATA -> concat_fr.MS:DATA
    # average a lot to speed up and increase S/N -> fast scalarphase is gone, so this should be fine.
    #logger.info('Create small concatenated MS (averaging)...')
    #lib_util.check_rm('concat_fr.MS')
    #s.add(f'DP3 {parset_dir}/DP3-avg.parset msin=concat_all.MS msout=concat_fr.MS msin.datacolumn=SMOOTHED_DATA \
    #                avg.timestep={int(small_timestep)} avg.freqresolution={small_freqres}', log='concat_fr_avg.log', commandType='DP3')
    #s.run(check=True)
    #MSs_fr = lib_ms.AllMSs(['concat_fr.MS'], s, check_flags=False)
    # predict the MODEL_DATA
    #logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
    #MSs_fr.run(f"DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}",
    #           log="$nameMS_pre.log", commandType="DP3")

    # Solve concat_fr.MS:DATA (only solve)
    # TODO test sol.rotationdiagonalmode=scalarphase
    #logger.info(f'Calibrating FR...')
    # Need to solve for rot+diag and not just rot since we can have phase offsets from the preliminary iono!
    #MSs_fr.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/fr.h5 \
    #           sol.mode=rotation+diagonal sol.rotationdiagonalmode=scalar\
    #           sol.solint=1 sol.nchan=1', log='$nameMS_solFR.log', commandType="DP3")
    #lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_fr.getListStr()],
    #                    [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-plot-rot.parset',
    #                     parset_dir + '/losoto-plot-scalaramp.parset', parset_dir + '/losoto-fr.parset'])
    #lib_util.check_rm('concat_fr.MS')

    # Solve concat_all.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating FR...')
    # We solve for rot+diag or rot+scalar here and not just rot since we can have phase offsets from the preliminary iono!!
    # TODO can we add smoothnessconstraint here?
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 \
               sol.mode=rotation+diagonal sol.rotationdiagonalmode=scalarphase\
               sol.solint={small_timestep} sol.nchan={small_freqstep}', log='$nameMS_solFR.log', commandType="DP3")

    # TODO add residual rotation plot after FR fit as soon as this option is present in LoSoTo!
    lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-plot-rot.parset',
                         parset_dir + '/losoto-fr.parset'])

### DONE

#################################################

# 3: calibrate iono + clock
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

    # Solve concat_all.MS:DATA (only solve on CS-CS BL)
    # not more smoothing since in rare case some CS have strong delay!
    logger.info('Calibrating IONO (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.minvisratio=0.1 \
                        sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.datause=single sol.solint=15 sol.nchan=1 msin.baseline="CS*&CS*" \
                        sol.smoothnessconstraint=0.5e6 sol.uvlambdamin={uvlambdamin}', log='$nameMS_solIONO_CS.log',
                       commandType="DP3")

    lib_util.run_losoto(s, 'iono-cs', [ms + '/iono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-ionoCS.parset'])

    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA (unit correction for others)
    logger.info('Iono correction (Core Stations)...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")

    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth3')

    # Phasing up the cose stations SMOOTHED_DATA -> concat_all-phaseupIONO.MS:DATA
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseupIONO.MS')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-phaseup.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        msout=concat_all-phaseupIONO.MS', log='$nameMS_phaseup.log', commandType="DP3")

    MSs_concat_phaseupIONO = lib_ms.AllMSs(['concat_all-phaseupIONO.MS'], s, check_flags=False)

    logger.info(f'Add model of {calname} from {os.path.basename(skymodel)} to MODEL_DATA...')
    os.system('cp -r %s %s' % (skymodel, MSs_concat_phaseupIONO.getListObj()[0].pathMS))
    MSs_concat_phaseupIONO.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={skymodel} pre.sources={calname}',
                               log="$nameMS_pre.log", commandType="DP3")

    # Solve concat_all-phaseupIONO.MS:DATA (only solve)
    logger.info('Calibrating IONO (distant stations)...')
    MSs_concat_phaseupIONO.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS msin.datacolumn=DATA \
                           sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.datause=single sol.minvisratio \
                           sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=0.1e6 sol.smoothnessreffrequency=54e6', \
                           log='$nameMS_solIONO.log', commandType="DP3")
   
    if (min(MSs_concat_phaseupIONO.getFreqs()) < 35.e6):
        lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-iono3rd.parset'])
    else:
        lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-iono.parset'])
### DONE

######################################################
# 4: find BP
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
                        cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")
    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                cor.correction=phase000', log='$nameMS_corIONO_CS.log', commandType="DP3")
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve) against FR-corrupted MODEL_DATA
    logger.info('Calibrating BP...')
    timestep = int(np.rint(60 / tint))  # brings down to 60s
    MSs_concat_all.run(f'DP3 {parset_dir}/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/bp.h5 sol.mode=diagonal sol.datause=dual \
                        sol.modeldatacolumns=[MODEL_DATA_FRCOR] sol.solint={str(timestep)} sol.nchan=1',
                       log='$nameMS_solBP.log', commandType="DP3")

    flag_parset = '/losoto-flag-sparse.parset' if sparse_sb else '/losoto-flag.parset'
    lib_util.run_losoto(s, 'bp', [ms + '/bp.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + flag_parset, parset_dir + '/losoto-bp.parset'])
### DONE

if debugplots:
    # Pol align correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 \
                   cor.correction=polalign', log='$nameMS_corPA2.log', commandType="DP3")
    
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth3')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pa', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', 
                           log='$nameMS_beam2.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth3')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeam', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('BP correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-bp.h5 \
            cor.correction=fulljones cor.soltab=[amplitudeSmooth,phaseNull] cor.updateweights=False',
                           log='$nameMS_corBP.log', commandType="DP3")
    
    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth3')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambp', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 \
                   cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth3')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambpfr', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Correct iono CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                   cor.correction=phase000', log='$nameMS_corIONO2_CS.log', commandType="DP3")
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                   cor.correction=phase000', log='$nameMS_corIONO2.log', commandType="DP3")

    MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', nofreq=True, logstr='smooth3')

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambpfriono', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

with w.if_todo('compressing_h5'):
    logger.info('Compressing caltables...')
    # os.system('cp cal-pa.h5 fullcal-pa.h5')
    # os.system('cp cal-fr.h5 fullcal-fr.h5') # no need to keep orig
    # os.system('cp cal-bp.h5 fullcal-bp.h5')
    # os.system('cp cal-iono.h5 fullcal-iono.h5')
    s.add('losoto -d sol000/phase000 cal-pa.h5', log='losoto-final.log', commandType="python")
    # s.add('losoto -d sol000/amplitude000 cal-pa.h5', log='losoto-final.log', commandType="python")
    s.add('losoto -d sol000/rotation000 cal-pa.h5', log='losoto-final.log', commandType="python")
    #s.add('losoto -d sol000/phaseResid000 cal-pa.h5', log='losoto-final.log', commandType="python")

    s.add('losoto -d sol000/phase000 cal-fr.h5', log='losoto-final.log', commandType="python")
    #s.add('losoto -d sol000/phaseResid000 cal-fr.h5', log='losoto-final.log', commandType="python")

    s.add('losoto -d sol000/amplitude000 cal-bp.h5', log='losoto-final.log', commandType="python")
    s.add('losoto -d sol000/phase000 cal-bp.h5', log='losoto-final.log', commandType="python")
    s.add('losoto -d sol000/amplitudeRes cal-bp.h5', log='losoto-final.log', commandType="python")

    s.add('losoto -d sol000/phase_offset000 cal-iono-cs.h5', log='losoto-final.log', commandType="python")
    s.add('losoto -d sol000/phaseResid000 cal-iono-cs.h5', log='losoto-final.log', commandType="python")

    s.add('losoto -d sol000/phase_offset000 cal-iono.h5', log='losoto-final.log', commandType="python")
    s.add('losoto -d sol000/phaseResid000 cal-iono.h5', log='losoto-final.log', commandType="python")
    
    s.run(maxThreads=1, check=True) # final check on losoto-final.log

    os.system('h5repack cal-pa.h5 cal-pa-compressed.h5; mv cal-pa-compressed.h5 cal-pa.h5')
    os.system('h5repack cal-fr.h5 cal-fr-compressed.h5; mv cal-fr-compressed.h5 cal-fr.h5')
    os.system('h5repack cal-bp.h5 cal-bp-compressed.h5; mv cal-bp-compressed.h5 cal-bp.h5')
    os.system('h5repack cal-iono-cs.h5 cal-iono-cs-compressed.h5; mv cal-iono-cs-compressed.h5 cal-iono-cs.h5')
    os.system('h5repack cal-iono.h5 cal-iono-compressed.h5; mv cal-iono-compressed.h5 cal-iono.h5')

    # remove unnecessary tables
    lib_util.check_rm('cal-preiono.h5 cal-preiono-cs.h5')

### DONE

# a debug image
if imaging:
    with w.if_todo('cal_leakage'):
        # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
        logger.info('BP correction...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-bp.h5 \
            cor.correction=amplitudeSmooth  cor.updateweights=True', log='$nameMS_corBP.log', commandType="DP3")

        debug_imaging(MSs_concat_all, 'beforefj')

        # Smooth data concat_all.MS:CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
        # unclear if we should smooth here
        MSs_concat_all.run_Blsmooth(incol='CORRECTED_DATA', logstr='smooth')  # now we can also smooth in time

        # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating Leakage...')
        timestep = int(np.rint(60 / tint))  # brings down to 60s
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS sol.h5parm=$pathMS/fj.h5 sol.mode=fulljones \
                            sol.modeldatacolumns=[MODEL_DATA_FRCOR] sol.solint=' + str(timestep) + ' sol.nchan=' + str(nchan),
                           log='$nameMS_solFJ.log', commandType="DP3")

        lib_util.run_losoto(s, 'fj', [ms + '/fj.h5' for ms in MSs_concat_all.getListStr()],
                            [parset_dir + '/losoto-plot-fullj.parset'])

    with w.if_todo('cal_imaging'):

        # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Leakage correction...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fj.h5 \
            cor.correction=fulljones cor.soltab=[amplitude000,phase000] cor.updateweights=False',
                           log='$nameMS_corBP.log', commandType="DP3")

        # Correct FR concat_all-phaseup.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cat cor.parmdb=cal-fr.h5 \
                       cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")

        debug_imaging(MSs_concat_all, 'final')

    ### DONE

logger.info('Cleaning up...')
os.system('rm -r *MS')

w.alldone()
