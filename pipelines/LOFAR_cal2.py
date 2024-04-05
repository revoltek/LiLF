#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation.
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import sys, os, glob, re
import casacore.tables as pt
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_h5

logger_obj = lib_log.Logger('pipeline-cal2')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry=False)
w = lib_util.Walker('pipeline-cal2.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: ' + str(dict(parset['LOFAR_cal2'])))
parset_dir = parset.get('LOFAR_cal2', 'parset_dir')
data_dir = parset.get('LOFAR_cal2', 'data_dir')
skymodel = parset.get('LOFAR_cal2', 'skymodel')
imaging = parset.getboolean('LOFAR_cal2', 'imaging')
fillmissingedges = parset.getboolean('LOFAR_cal2', 'fillmissingedges')
bl2flag = parset.get('flag', 'stations')
debugplots = False

#############################################################

def debug_imaging(MSs, suffix):
    if not os.path.exists('img'):
        os.makedirs('img')

    imgsizepix = 1024

    logger.info('Cleaning...')
    imagename = f'img/cal-{suffix}'
    lib_util.run_wsclean(s, 'wsclean.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix,
                         scale='3arcsec', auto_mask=5, # local_rms='', local_rms_method='rms-with-min',
                         weight='briggs -0.3', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.5,
                         baseline_averaging='', auto_threshold=2, join_channels='', fit_spectral_pol=5,
                         channels_out=12)

MSs = lib_ms.AllMSs(glob.glob(data_dir + '/*MS'), s, check_flags=False)

### This part is done so that missing subbands get concatenated correctly.
for i, msg in enumerate(np.array_split(sorted(glob.glob(data_dir+'/*MS')), 1)):
    if fillmissingedges:
        min_nu = pt.table(MSs.getListStr()[0]).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MIN']
        max_nu = pt.table(MSs.getListStr()[0]).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MAX']
    else:
        min_nu = min(MSs.getFreqs())/1e6
        max_nu = max(MSs.getFreqs())/1e6
    print(min_nu,max_nu)
    num_init = lib_util.lofar_nu2num(min_nu) + 1  # +1 because FREQ_MIN/MAX somewhat have the lowest edge of the SB freq
    num_fin = lib_util.lofar_nu2num(max_nu) + 1
    prefix = re.sub('SB[0-9]*.MS', '', msg[0])
    msg = []
    for j in range(num_init, num_fin + 1):
        msg.append(prefix + 'SB%03i.MS' % j)

if min(MSs.getFreqs()) < 35.e6:
    iono3rd = True
    logger.debug('Include iono 3rd order.')
else:
    iono3rd = False

if skymodel == '':  # default case
    if MSs.hasIS:
        skymodel = os.path.dirname(__file__) + '/../models/calib-highres.skydb'
    else:
        skymodel = os.path.dirname(__file__) + '/../models/calib-simple.skydb'

calname = MSs.getListObj()[0].getNameField()
nchan = MSs.mssListObj[0].getNchan()
tint = MSs.mssListObj[0].getTimeInt()

#with w.if_todo('concat_core'):
#    freqstep = nchan  # brings down to 1ch/sb for
#    timestep = int(np.rint(8 / tint))  # brings down to 8s
#    # concat all SBs
#    # SB.MS:DATA -> concat.MS:DATA
#    logger.info('Concatenating data core...')
#    lib_util.check_rm('concat_core.MS')
#    s.add('DP3 ' + parset_dir + '/DP3-concat.parset msin="[' + ','.join(msg) + ']" msout=concat_core.MS \
#                      msin.baseline="CS*&" avg.freqstep=' + str(freqstep) + ' avg.timestep=' + str(timestep),
#          log='concat.log', commandType='DP3')
#
#    s.run(check=True)
### DONE
#MSs_concat_core = lib_ms.AllMSs(['concat_core.MS'], s, check_flags=False)

#with w.if_todo('predict_core'):
#    # predict to save time ms:MODEL_DATA
#    logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
#    os.system('cp -r %s %s' % (skymodel, MSs_concat_core.getListObj()[0].pathMS))
#    MSs_concat_core.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/"
#                        + os.path.basename(skymodel) + " pre.sources=" + calname, log="$nameMS_pre.log",
#                        commandType="DP3")
### DONE

with w.if_todo('concat_all'):
    freqstep = 1  # keep all channels
    timestep = int(np.rint(4 / tint))  # brings down to 4s
    # concat all SBs
    # SB.MS:DATA -> concat.MS:DATA
    logger.info('Concatenating data all...')
    lib_util.check_rm('concat_all.MS')
    s.add('DP3 ' + parset_dir + '/DP3-concat.parset msin="[' + ','.join(msg) + ']" msout=concat_all.MS \
              msin.baseline="*&" avg.freqstep=' + str(freqstep) + ' avg.timestep=' + str(timestep),
          log='concat.log', commandType='DP3')
    s.run(check=True)
### DONE
MSs_concat_all = lib_ms.AllMSs(['concat_all.MS'], s, check_flags=False)

######################################################
# flag bad stations, flags will propagate
with w.if_todo('flag'):
    logger.info("Flagging...")
    #MSs_concat_core.run("DP3 " + parset_dir + "/DP3-flag.parset msin=$pathMS ant.baseline=\"" + bl2flag + "\"",
    #                    log="$nameMS_flag.log", commandType="DP3")
    MSs_concat_all.run("DP3 " + parset_dir + "/DP3-flag.parset msin=$pathMS ant.baseline=\"" + bl2flag + "\"",
                       log="$nameMS_flag.log", commandType="DP3")
    # extend flags
    logger.info('Remove bad time/freq stamps...')
    #MSs_concat_core.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
    MSs_concat_all.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

### DONE

#with w.if_todo('phaseupcore'):
#    # Phasing up the cose stations
#    logger.info('Phasing up Core Stations...')
#    lib_util.check_rm('concat_all-phaseup.MS')
#    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-phaseup.parset msin=$pathMS msin.datacolumn=DATA msout=concat_all-phaseup.MS',
#                       log='$nameMS_phaseup.log', commandType="DP3")
#
### DONE
#
#MSs_concat_phaseup = lib_ms.AllMSs(['concat_all-phaseup.MS'], s, check_flags=False)
#
#with w.if_todo('predict_all-phaseup'):
#    # predict to save time ms:MODEL_DATA
#    logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
#    os.system('cp -r %s %s' % (skymodel, MSs_concat_phaseup.getListObj()[0].pathMS))
#    MSs_concat_phaseup.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/"
#                           + os.path.basename(skymodel) + " pre.sources=" + calname, log="$nameMS_pre.log",
#                           commandType="DP3")

### DONE

if debugplots:
    # predict to save time ms:MODEL_DATA
    logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
    os.system('cp -r %s %s' % (skymodel, MSs_concat_all.getListObj()[0].pathMS))
    MSs_concat_all.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/"
                           + os.path.basename(skymodel) + " pre.sources=" + calname, log="$nameMS_pre.log",
                           commandType="DP3")    
    logger.info('BL-smooth...')
    MSs_concat_all.run(
        f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth3.log', commandType='python', maxThreads=1)
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
            sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
            sol.solint=1 sol.nchan=1', \
            log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-init', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

###################################################
# 1: find PA

with w.if_todo('cal_pa'):
  
    # Smooth data concat_all-all DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs_concat_all.run(
        f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth1.log', commandType='python', maxThreads=1)

    # Get phase diff SMOOTHED_DATA -> SMOOTHED_DATA
    logger.info('Get phase difference...')
    MSs_concat_all.run('taql "UPDATE $pathMS SET\
                             SMOOTHED_DATA[,0]=0.5*EXP(1.0i*(PHASE(SMOOTHED_DATA[,0])-PHASE(SMOOTHED_DATA[,3]))), \
                             SMOOTHED_DATA[,3]=SMOOTHED_DATA[,0], SMOOTHED_DATA[,1]=0+0i, SMOOTHED_DATA[,2]=0+0i"', 
                             log='$nameMS_taql_phdiff.log', commandType='general')
    
    logger.info('Creating MODEL_DATA...')
    MSs_concat_all.addcol('MODEL_DATA', 'DATA', usedysco=False)  # need this to make sure no dysco, if we have dyso we cannot set values to zero
    MSs_concat_all.run('taql "UPDATE $pathMS SET MODEL_DATA[,0]=0.5+0i, MODEL_DATA[,1]=0.0+0i, \
                             MODEL_DATA[,2]=0.0+0i, MODEL_DATA[,3]=0.5+0i"', 
                             log='$nameMS_taql_model.log', commandType='general')
    
    # Solve cal_SB.MS:DATA against MODEL_DATA (only solve)
    logger.info('Calibrating PA...')
    freqstep = nchan  # brings down to 1ch/sb for
    timestep = int(np.rint(600 / tint))  # brings down to 10m
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS \
                             sol.h5parm=$pathMS/pa.h5 sol.mode=phaseonly sol.solint='+str(timestep)+' sol.nchan='+str(freqstep)+' \
                             sol.smoothnessconstraint=2e6',
                             log='$nameMS_solPA.log', commandType="DP3")

    lib_util.run_losoto(s, 'pa', [ms + '/pa.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-pa.parset'])
    
### DONE

########################################################
# 2: find FR
with w.if_todo('cal_fr'):

    # Pol align correction concat_all.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', 
                       log='$nameMS_beam.log', commandType="DP3")

    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs_concat_all.run(
        f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth2.log', commandType='python', maxThreads=1)
    
    # Convert to circular SMOOTHED_DATA -> SMOOTHED_DATA
    logger.info('Converting to circular...')
    MSs_concat_all.run('mslin2circ.py -s -i $pathMS:SMOOTHED_DATA -o $pathMS:SMOOTHED_DATA',
                           log='$nameMS_lincirc.log', commandType='python', maxThreads=1)
    
    # Phasing up the cose stations SMOOTHED_DATA -> concat_all-phaseupFR.MS:DATA
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseupFR.MS')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-phaseup.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        msout=concat_all-phaseupFR.MS', log='$nameMS_phaseup.log', commandType="DP3")

    MSs_concat_phaseupFR = lib_ms.AllMSs(['concat_all-phaseupFR.MS'], s, check_flags=False)

    # Get circular phase diff DATA -> DATA
    logger.info('Get circular phase difference...')
    MSs_concat_phaseupFR.run('taql "UPDATE $pathMS SET\
                             DATA[,0]=0.5*EXP(1.0i*(PHASE(DATA[,0])-PHASE(DATA[,3]))), \
                             DATA[,3]=DATA[,0], DATA[,1]=0+0i, DATA[,2]=0+0i"', 
                             log='$nameMS_taql_phdiff.log', commandType='general')
    logger.info('Creating MODEL_DATA...')  # take from MODEL_DATA but overwrite
    MSs_concat_phaseupFR.addcol('MODEL_DATA', 'DATA', usedysco=False)  # need this to make sure no dysco, if we have dyso we cannot set values to zero
    MSs_concat_phaseupFR.run('taql "UPDATE $pathMS SET MODEL_DATA[,0]=0.5+0i, MODEL_DATA[,1]=0.0+0i, \
                             MODEL_DATA[,2]=0.0+0i, MODEL_DATA[,3]=0.5+0i"', 
                             log='$nameMS_taql_model.log', commandType='general')
    # Solve cal_SB.MS:DATA against MODEL_DATA (only solve)
    logger.info('Calibrating FR...')
    MSs_concat_phaseupFR.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA \
                             sol.h5parm=$pathMS/fr.h5 sol.mode=phaseonly sol.solint=4 sol.nchan=4 \
                             sol.smoothnessconstraint=2e6 sol.smoothnessreffrequency=54e6',
                             log='$nameMS_solFR.log', commandType="DP3")

    lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs_concat_phaseupFR.getListStr()],
                        [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-fr.parset'])

### DONE 

#################################################
with w.if_todo('predict_all'):
    logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
    os.system('cp -r %s %s' % (skymodel, MSs_concat_all.getListObj()[0].pathMS))
    MSs_concat_all.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/"
                       + os.path.basename(skymodel) + " pre.sources=" + calname, log="$nameMS_pre.log",
                       commandType="DP3")

# 3: find iono
with w.if_todo('cal_iono'):

    # Pol align correction concat_all-phaseup.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam concat_all-phaseup.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', 
                           log='$nameMS_beam.log', commandType="DP3")
    # Correct FR concat_all-phaseup.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-fr.h5 \
                   cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")

    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs_concat_all.run(
        f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth2.log', commandType='python', maxThreads=1)

    # Solve cal_SB.MS:DATA CS-CS baselines(only solve)
    logger.info('Calibrating IONO (Core Stations)...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.solint=16 sol.nchan=1 msin.baseline="CS*&CS*" \
                        sol.smoothnessconstraint=4e6 sol.smoothnessreffrequency=54e6', log='$nameMS_solIONO_CS.log',
                       commandType="DP3")

    lib_util.run_losoto(s, 'iono-cs', [ms + '/iono.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-iono.parset'])

    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA (unit correction for others)
    logger.info('Iono correction (Core Stations)...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                cor.correction=phaseOrig000', log='$nameMS_corIONO_CS.log', commandType="DP3")

    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs_concat_all.run(
        f'BLsmooth.py -r -q -c {16 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth2.log', commandType='python', maxThreads=1)

    # Phasing up the cose stations SMOOTHED_DATA -> concat_all-phaseupFR.MS:DATA
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseupIONO.MS')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-phaseup.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA \
                        msout=concat_all-phaseupIONO.MS', log='$nameMS_phaseup.log', commandType="DP3")
    
    MSs_concat_phaseupIONO = lib_ms.AllMSs(['concat_all-phaseupIONO.MS'], s, check_flags=False)

    logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
    os.system('cp -r %s %s' % (skymodel, MSs_concat_phaseupIONO.getListObj()[0].pathMS))
    MSs_concat_phaseupIONO.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/"
                           + os.path.basename(skymodel) + " pre.sources=" + calname, log="$nameMS_pre.log",
                           commandType="DP3")   

    # Solve cal_SB.MS:DATA (only solve)
    logger.info('Calibrating IONO (distant stations)...')
    MSs_concat_phaseupIONO.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA \
                           sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase \
                           sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=0.1e6 sol.smoothnessreffrequency=54e6', \
                           log='$nameMS_solIONO.log', commandType="DP3")
   
    if iono3rd:
        lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-iono3rd.parset'])
    else:
        lib_util.run_losoto(s, 'iono', [ms + '/iono.h5' for ms in MSs_concat_phaseupIONO.getListStr()],
                            [parset_dir + '/losoto-ref-ph.parset', parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-iono.parset'])
### DONE

# TEST
#MSs_concat_phaseup = MSs_concat_all

######################################################
# 4: find BP
with w.if_todo('cal_bp'):
    ## Pol align correction concat_all-phaseup.MS:DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                   cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")
    # Correct beam concat_all-phaseup.MS:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', 
                           log='$nameMS_beam.log', commandType="DP3")
    # # Correct FR concat_all-phaseup.MS:CORRECTED_DATA -> CORRECTED_DATA
    # logger.info('Faraday rotation correction...')
    # MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cat cor.parmdb=cal-fr.h5 \
    #                cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
    # Corrupt FR concat_all-phaseup.MS:MODEL_DATA -> MODEL_DATA
    logger.info('Faraday rotation corruption (MODEL_DATA)...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")

    # Correct iono concat_all:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                cor.correction=phaseOrig000', log='$nameMS_corIONO_CS.log', commandType="DP3")
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                cor.correction=phaseOrig000', log='$nameMS_corIONO.log', commandType="DP3")
    
    # Smooth data concat_all-phaseup.MS:CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs_concat_all.run(
        f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth3.log', commandType='python', maxThreads=1)

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating BP...')
    freqstep = nchan  # brings down to 1ch/sb for
    timestep = int(np.rint(60 / tint))  # brings down to 60s
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
            sol.h5parm=$pathMS/bp.h5 sol.mode=fulljones \
            sol.solint='+str(timestep)+' sol.nchan='+str(freqstep),
            log='$nameMS_solBP.log', commandType="DP3")
    
    lib_util.run_losoto(s, 'bp', [ms + '/bp.h5' for ms in MSs_concat_all.getListStr()],
                        [parset_dir + '/losoto-flag.parset', parset_dir + '/losoto-plot-fullj.parset',
                         parset_dir + '/losoto-bp.parset', parset_dir + '/losoto-flagstations.parset'])
### DONE
    
if debugplots:
    # we need to predict again since we corrupted the MODEL_DATA (we could also correct the model again or give the corrupted model a new name)
    logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
    MSs_concat_all.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/"
                       + os.path.basename(skymodel) + " pre.sources=" + calname, log="$nameMS_pre.log",
                       commandType="DP3")
    # Pol align correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 \
                   cor.correction=polalign', log='$nameMS_corPA2.log', commandType="DP3")
    
    logger.info('BL-smooth...')
    MSs_concat_all.run(
            f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
            log='$nameMS_smooth3.log', commandType='python', maxThreads=1)
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pa', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs_concat_all.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', 
                           log='$nameMS_beam2.log', commandType="DP3")

    logger.info('BL-smooth...')
    MSs_concat_all.run(
            f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
            log='$nameMS_smooth3.log', commandType='python', maxThreads=1)
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
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
    
    logger.info('BL-smooth...')
    MSs_concat_all.run(
            f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
            log='$nameMS_smooth3.log', commandType='python', maxThreads=1)
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambp', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 \
                   cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DP3")

    logger.info('BL-smooth...')
    MSs_concat_all.run(
            f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
            log='$nameMS_smooth3.log', commandType='python', maxThreads=1)
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
                sol.h5parm=$pathMS/test.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1', \
                log='$nameMS_solTEST.log', commandType="DP3")
    lib_util.run_losoto(s, 'test-pabeambpfr', [ms + '/test.h5' for ms in MSs_concat_all.getListStr()],
                [parset_dir + '/losoto-plot-fullj.parset', parset_dir + '/losoto-bp.parset'])

    # Correct iono CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Iono correction...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono-cs.h5 \
                   cor.correction=phaseOrig000', log='$nameMS_corIONO2_CS.log', commandType="DP3")
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 \
                   cor.correction=phaseOrig000', log='$nameMS_corIONO2.log', commandType="DP3")

    logger.info('BL-smooth...')
    MSs_concat_all.run(
            f'BLsmooth.py -r -q -c {20 if MSs.hasIS else 1} -n {s.max_processors} -f {.2e-3 if MSs.hasIS else 1e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
            log='$nameMS_smooth3.log', commandType='python', maxThreads=1)
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating test...')
    MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
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
    s.add('losoto -d sol000/amplitude000 cal-pa.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/phase000 cal-pa.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/phaseOrig000 cal-pa.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-pa.h5 cal-pa-compressed.h5; mv cal-pa-compressed.h5 cal-pa.h5')

    s.add('losoto -d sol000/phase000 cal-fr.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/phaseOrig000 cal-fr.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-fr.h5 cal-fr-compressed.h5; mv cal-fr-compressed.h5 cal-fr.h5')

    s.add('losoto -d sol000/amplitude000 cal-bp.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/amplitudeRes cal-bp.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/phase000 cal-bp.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-bp.h5 cal-bp-compressed.h5; mv cal-bp-compressed.h5 cal-bp.h5')

    s.add('losoto -d sol000/phase_offset000 cal-iono-cs.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-iono-cs.h5 cal-iono-cs-compressed.h5; mv cal-iono-cs-compressed.h5 cal-iono-cs.h5')
    s.add('losoto -d sol000/phase_offset000 cal-iono.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-iono.h5 cal-iono-compressed.h5; mv cal-iono-compressed.h5 cal-iono.h5')

### DONE

# a debug image
if imaging:

    with w.if_todo('cal_imaging'):

        # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
        logger.info('BP correction...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-bp.h5 \
            cor.correction=fulljones cor.soltab=[amplitudeSmooth,phaseNull] cor.updateweights=True',
                           log='$nameMS_corBP.log', commandType="DP3")

        # Correct FR concat_all-phaseup.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        MSs_concat_all.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cat cor.parmdb=cal-fr.h5 \
                       cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")

        debug_imaging(MSs_concat_all, 'final')

    ### DONE

logger.info("Done.")
