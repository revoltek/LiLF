#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation.
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import sys, os, glob, re
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_h5
logger_obj = lib_log.Logger('pipeline-cal')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-cal.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_cal2'])))
parset_dir = parset.get('LOFAR_cal2','parset_dir')
data_dir = parset.get('LOFAR_cal2','data_dir')
skymodel = parset.get('LOFAR_cal2','skymodel')
imaging = parset.getboolean('LOFAR_cal2','imaging')
bl2flag = parset.get('flag','stations')

#############################################################

MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s, check_flags=False)

if min(MSs.getFreqs()) < 35.e6:
    iono3rd = True
    logger.debug('Include iono 3rd order.')
else: iono3rd = False

if skymodel == '':  # default case
    if MSs.isLBA:
        if MSs.hasIS:
            skymodel = os.path.dirname(__file__) + '/../models/calib-highres.skydb'
        else:
            skymodel = os.path.dirname(__file__) + '/../models/calib-simple.skydb'
    elif MSs.isHBA:
        skymodel = os.path.dirname(__file__) + '/../models/calib-hba.skydb'

calname = MSs.getListObj()[0].getNameField()
nchan = MSs.mssListObj[0].getNchan()
tint = MSs.mssListObj[0].getTimeInt()

#if nchan > 4:
#    base_nchan = int(np.rint(nchan / 4)) # brings down to 4ch/sb for 
#else: base_nchan = 1
#if tint < 4:
#    base_solint = int(np.rint(4/tint)) # this is 1 for dutch observations and 2 for IS observations
#else: base_solint = 1

for c in ['core','all']:

    with w.if_todo('concat_%s' % c):
        
        if c == "core":
            freqstep = nchan # brings down to 1ch/sb for 
            timestep = int(np.rint(8/tint)) # brings down to 8s
            # concat all SBs
            # SB.MS:DATA -> concat.MS:DATA
            logger.info('Concatenate data...')
            lib_util.check_rm('concat_core.MS')
            s.add('DP3 '+parset_dir+'/DP3-concat.parset msin="['+','.join(MSs.getListStr())+']" msout=concat_core.MS \
                  msin.baseline="CS*&" avg.freqstep='+str(freqstep)+' avg.timestep='+str(timestep), 
                  log='concat.log', commandType='DP3')
            s.run(check=True)
        else:
            freqstep = 1 # keep all channels 
            timestep = int(np.rint(4/tint)) # brings down to 4s
            # concat all SBs
            # SB.MS:DATA -> concat.MS:DATA
            logger.info('Concatenate data...')
            lib_util.check_rm('concat_all.MS')
            s.add('DP3 '+parset_dir+'/DP3-concat.parset msin="['+','.join(MSs.getListStr())+']" msout=concat_all.MS \
                  msin.baseline="*&" avg.freqstep='+str(freqstep)+' avg.timestep='+str(timestep), 
                  log='concat.log', commandType='DP3')
            s.run(check=True)

    ### DONE

    if c == "core": MSs_concat = lib_ms.AllMSs( ['concat_core.MS'], s, check_flags=False )
    else: MSs_concat = lib_ms.AllMSs( ['concat_all.MS'], s, check_flags=False )
    os.system('cp -r %s %s' % (skymodel, MSs_concat.getListObj()[0].pathMS))

    if c == "all":
        with w.if_todo('correctCS'):
            # Correct CS
            logger.info('Correcting CS...')
            # Pol align correction DATA -> CORRECTED_DATA
            logger.info('Polalign correction...')
            MSs_concat.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
                           cor.parmdb=cal-pa_core.h5 cor.correction=polalign', log='$nameMS_corPA2.log', commandType="DP3")

            # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
            logger.info('AmpBP correction...')
            MSs_concat.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.parmdb=cal-amp_core.h5 \
                cor.correction=amplitudeSmooth cor.updateweights=True', log='$nameMS_corAMP.log', commandType="DP3")

            # Beam correction CORRECTED_DATA -> CORRECTED_DATA
            #logger.info('Beam correction...')
            #MSs_concat.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', 
            #               log='$nameMS_beam2.log', commandType="DP3")

            # Correct FR CORRECTED_DATA -> CORRECTED_DATA
            logger.info('Faraday rotation correction...')
            MSs_concat.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr_core.h5 \
                           cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DP3")

            # Correct Iono CORRECTED_DATA -> CORRECTED_DATA
            logger.info('IONO correction...')
            MSs_concat.run("DP3 " + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono_core.h5 \
                           cor.correction=phaseOrig000', log='$nameMS_corIONO.log', commandType="DP3")
    
        ### DONE

        with w.if_todo('phaseupCS'):
            # Phasing up the cose stations
            logger.info('Phasing up Core Stations...')
            lib_util.check_rm('concat_all-phaseup.MS')
            MSs_concat.run('DP3 ' + parset_dir + '/DP3-phaseup.parset msin=$pathMS msout=concat_all-phaseup.MS', log='$nameMS_phaseup.log', commandType="DP3")
            
        ### DONE
        
        MSs_concat = lib_ms.AllMSs( ['concat_all-phaseup.MS'], s, check_flags=False )

    ######################################################
    # flag bad stations, flags will propagate
    with w.if_todo('flag_%s' % c):
        logger.info("Flagging...")
        MSs_concat.run("DP3 " + parset_dir + "/DP3-flag.parset msin=$pathMS ant.baseline=\"" + bl2flag+"\"", 
                    log="$nameMS_flag.log", commandType="DP3")

        #logger.info("Flagging...")
        #flag_strat = '/HBAdefaultwideband.lua' if MSs.isHBA else '/LBAdefaultwideband.lua'
        #MSs_concat.run("DP3 " + parset_dir + "/DP3-flag2.parset msin=$pathMS aoflagger.strategy=" + parset_dir + flag_strat, 
        #            log="$nameMS_flag.log", commandType="DP3")

        # extend flags
        logger.info('Remove bad time/freq stamps...')
        MSs_concat.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    ### DONE

    with w.if_todo('predict_%s' % c):
        # predict to save time ms:MODEL_DATA
        logger.info('Add model of %s from %s to MODEL_DATA...' % (calname, os.path.basename(skymodel)))
        MSs_concat.run("DP3 " + parset_dir + "/DP3-predict.parset msin=$pathMS pre.sourcedb=$pathMS/" + os.path.basename(skymodel) + " pre.sources=" + calname, \
                log="$nameMS_pre.log", commandType="DP3")

    ### DONE

    ###################################################
    # 1: find PA

    with w.if_todo('cal_pa_%s' % c):
        # Smooth data DATA -> SMOOTHED_DATA (BL-based smoothing)
        logger.info('BL-smooth...')
        MSs_concat.run(f'BLsmooth.py -r -q -c 1 -n 8 -f {1e-2 if MSs.isLBA else .2e-3} -i DATA -o SMOOTHED_DATA $pathMS', 
                       log='$nameMS_smooth1.log', commandType='python', maxThreads=8)

        # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating PA...')
        if c == 'core': smoothnessconstraint = 1e6
        else: smoothnessconstraint = 0.2e6
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/pa_'+c+'.h5 sol.mode=rotation+diagonal \
                sol.solint=1 sol.smoothnessconstraint='+str(smoothnessconstraint)+' sol.smoothnessreffrequency=54e6',
                log='$nameMS_solPA.log', commandType="DP3")

        lib_util.run_losoto(s, 'pa_'+c, [ms+'/pa_'+c+'.h5' for ms in MSs_concat.getListStr()],
                [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

        # Pol align correction DATA -> CORRECTED_DATA
        logger.info('Polalign correction...')
        MSs_concat.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa_'+c+'.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DP3")

    ### DONE

    ########################################################
    # 2: find FR
    with w.if_todo('cal_fr_%s' % c):
        # Beam correction CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Beam correction...')
        MSs_concat.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log',
                commandType="DP3")

        # Smooth data CORRECTED_DATA -> CIRC_PHASEDIFF_DATA (BL-based smoothing)
        logger.info('BL-smooth...')
        MSs_concat.addcol('CIRC_PHASEDIFF_DATA', 'CORRECTED_DATA', usedysco=False) # need this to make sure no dysco, if we have dyso we cannot set values to zero
        MSs_concat.run(f'BLsmooth.py -r -q -c 1 -n 8 -f {1e-2 if MSs.isLBA else .2e-3} -i CIRC_PHASEDIFF_DATA -o CIRC_PHASEDIFF_DATA $pathMS', log='$nameMS_smooth2.log',
                commandType='python', maxThreads=8)

        logger.info('Converting to circular...')
        MSs_concat.run('mslin2circ.py -s -i $pathMS:CIRC_PHASEDIFF_DATA -o $pathMS:CIRC_PHASEDIFF_DATA', log='$nameMS_lincirc.log',
                commandType='python', maxThreads=2)
        # Get circular phase diff CIRC_PHASEDIFF_DATA -> CIRC_PHASEDIFF_DATA
        logger.info('Get circular phase difference...')
        MSs_concat.run('taql "UPDATE $pathMS SET\
        CIRC_PHASEDIFF_DATA[,0]=0.5*EXP(1.0i*(PHASE(CIRC_PHASEDIFF_DATA[,0])-PHASE(CIRC_PHASEDIFF_DATA[,3]))), \
        CIRC_PHASEDIFF_DATA[,3]=CIRC_PHASEDIFF_DATA[,0], \
        CIRC_PHASEDIFF_DATA[,1]=0+0i, \
        CIRC_PHASEDIFF_DATA[,2]=0+0i"',log='$nameMS_taql_phdiff.log', commandType='general')

        logger.info('Creating FR_MODEL_DATA...') # take from MODEL_DATA but overwrite
        MSs_concat.addcol('FR_MODEL_DATA', 'DATA', usedysco=False) # need this to make sure no dysco, if we have dyso we cannot set values to zero
        MSs_concat.run('taql "UPDATE $pathMS SET FR_MODEL_DATA[,0]=0.5+0i, FR_MODEL_DATA[,1]=0.0+0i, FR_MODEL_DATA[,2]=0.0+0i, \
        FR_MODEL_DATA[,3]=0.5+0i"', log='$nameMS_taql_frmodel.log', commandType='general')

        # Solve cal_SB.MS:CIRC_PHASEDIFF_DATA against FR_MODEL_DATA (only solve)
        logger.info('Calibrating FR...')
        if c == 'core': smoothnessconstraint = 5e6
        else: smoothnessconstraint = 2e6
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA \
        sol.modeldatacolumns=[FR_MODEL_DATA] sol.h5parm=$pathMS/fr_'+c+'.h5 sol.mode=phaseonly \
        sol.solint=4 sol.nchan=4 sol.smoothnessconstraint='+str(smoothnessconstraint)+' sol.smoothnessreffrequency=54e6', 
        log='$nameMS_solFR.log', commandType="DP3")
        
        lib_util.run_losoto(s, 'fr_'+c, [ms+'/fr_'+c+'.h5' for ms in MSs_concat.getListStr()], 
                            [parset_dir + '/losoto-fr.parset'])

        # remove columns
        MSs_concat.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, FR_MODEL_DATA"', log='$nameMS_taql_delcol.log', commandType='general')

        # Correct FR CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr_'+c+'.h5 cor.correction=rotationmeasure000',
                log='$nameMS_corFR.log', commandType="DP3")

    ### DONE

    ######################################################
    # 4: find BP

    with w.if_todo('cal_bp_%s' % c):

        # Apply best phases (phase000 has phase solution with no rotation and no align - both already applied)
        logger.info('Ph correction...')
        MSs_concat.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.parmdb=cal-pa_'+c+'.h5 \
                cor.correction=phase000', log='$nameMS_corPH.log', commandType="DP3")

        # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
        logger.info('BL-smooth...')
        MSs_concat.run(f'BLsmooth.py -r -q -c 1 -n 8 -f {1e-2 if MSs.isLBA else .2e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3.log',
                commandType ='python', maxThreads=8)

        # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
        # Use all MSs and smoothsol
        logger.info('Calibrating BP...')
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/amp_'+c+'.h5 sol.mode=diagonal \
                sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=1e6', \
                log='concat_solAMP.log', commandType="DP3")

        lib_util.run_losoto(s, 'amp_'+c, [ms+'/amp_'+c+'.h5' for ms in MSs_concat.getListStr()],
                [parset_dir + '/losoto-flag.parset', parset_dir+'/losoto-plot-amp.parset',
                parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-bp.parset'])

    ### DONE

    ####################################################
    # Re-do correcitons in right order

    with w.if_todo('apply_all_%s' % c):
        # Pol align correction DATA -> CORRECTED_DATA
        logger.info('Polalign correction...')
        MSs_concat.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa_'+c+'.h5 cor.correction=polalign', log='$nameMS_corPA2.log', commandType="DP3")

        # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
        logger.info('AmpBP correction...')
        MSs_concat.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.parmdb=cal-amp_'+c+'.h5 \
                cor.correction=amplitudeSmooth cor.updateweights=True', log='$nameMS_corAMP.log', commandType="DP3")

        # Beam correction CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Beam correction...')
        MSs_concat.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_beam2.log', commandType="DP3")

        # Correct FR CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr_'+c+'.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DP3")

        # Correct abs FR CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Absolute Faraday rotation correction...')
        #MSs.run('createRMh5parm.py $pathMS $pathMS/rme.h5', log='$nameMS_RME.log', commandType="general", maxThreads=1)
        #MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=$pathMS/rme.h5 cor.correction=RMextract', log='$nameMS_corRME.log', commandType="DP3")

    ### DONE

    #################################################
    # 4: find iono

    with w.if_todo('cal_iono_%s' % c):
        # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
        logger.info('BL-smooth...')
        MSs_concat.run(f'BLsmooth.py -r -q -c 1 -n 8 -f {1e-2 if MSs.isLBA else .2e-3} -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth4.log',
                commandType ='python', maxThreads=8)

        # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating IONO...')
        if c == 'core': smoothnessconstraint = 1e6
        else: smoothnessconstraint = 0.1e6
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/iono_'+c+'.h5 sol.mode=scalarphase \
                sol.solint=1 sol.nchan=1 sol.smoothnessconstraint='+str(smoothnessconstraint)+' sol.smoothnessreffrequency=54e6', \
                log='$nameMS_solIONO.log', commandType="DP3")

        if iono3rd:
            lib_util.run_losoto(s, 'iono_'+c, [ms+'/iono_'+c+'.h5' for ms in MSs_concat.getListStr()],
                [parset_dir+'/losoto-plot-scalarph.parset', parset_dir+'/losoto-iono3rd.parset'])
        elif MSs_concat.isHBA:
            lib_util.run_losoto(s, 'iono_'+c, [ms + '/iono_'+c+'.h5' for ms in MSs_concat.getListStr()],
                                [parset_dir + '/losoto-plot-scalarph.parset', parset_dir + '/losoto-iono-hba.parset'])
        else:
            lib_util.run_losoto(s, 'iono_'+c, [ms+'/iono_'+c+'.h5' for ms in MSs_concat.getListStr()],
                [parset_dir+'/losoto-plot-scalarph.parset', parset_dir+'/losoto-iono.parset'])

    ### DONE

with w.if_todo('compressing_h5'):
    logger.info('Compressing caltables...')
    #os.system('cp cal-pa.h5 fullcal-pa.h5')
    #os.system('cp cal-fr.h5 fullcal-fr.h5') # no need to keep orig
    #os.system('cp cal-amp.h5 fullcal-amp.h5')
    #os.system('cp cal-iono.h5 fullcal-iono.h5')
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

    s.add('losoto -d sol000/amplitude000 cal-amp.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/amplitudeRes cal-amp.h5', log='losoto-final.log', commandType="python")
    s.run()
    s.add('losoto -d sol000/phase000 cal-amp.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-amp.h5 cal-amp-compressed.h5; mv cal-amp-compressed.h5 cal-amp.h5')
    
    s.add('losoto -d sol000/phase_offset000 cal-iono.h5', log='losoto-final.log', commandType="python")
    s.run()
    os.system('h5repack cal-iono.h5 cal-iono-compressed.h5; mv cal-iono-compressed.h5 cal-iono.h5')

### DONE

# a debug image
if imaging:

    ####################################################
    # 5: find leak

    with w.if_todo('cal_leak'):

        # Correct all CORRECTED_DATA (PA, beam, FR, BP corrected) -> CORRECTED_DATA
        logger.info('IONO correction...')
        MSs_concat.run("DP3 " + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-iono_'+c+'.h5 \
            cor.correction=phaseOrig000', log='$nameMS_corIONO.log', commandType="DP3")

        # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
        logger.info('BL-smooth...')
        MSs_concat.run('BLsmooth.py -r -q -c 1 -n 8 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth4.log',
                commandType ='python', maxThreads=8)
        
        # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Calibrating LEAK...')
        MSs_concat.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/leak_'+c+'.h5 sol.mode=fulljones sol.solint=1 sol.nchan=1 sol.minvisratio=0.8 \
                sol.solint='+str(base_solint)+' sol.nchan='+str(base_nchan), \
                log='$nameMS_solLEAK.log', commandType="DP3")
    
        lib_util.run_losoto(s, 'leak_'+c, [ms+'/leak_'+c+'.h5' for ms in MSs_concat.getListStr()], \
                [parset_dir+'/losoto-flag-leak.parset', parset_dir+'/losoto-plot-fullj.parset'])

        logger.info('LEAK correction...')
        MSs_concat.run("DP3 " + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-leak_'+c+'.h5 cor.correction=fulljones\
            cor.soltab=[amplitude000,phase000]', log='$nameMS_corIONO.log', commandType="DP3")

    ### DONE

    with w.if_todo('cal_imaging'):
        logger.info("Imaging section:")
    
        lib_util.check_rm('img')
        os.makedirs('img')
    
        #imgsizepix = int(MSs.getListObj()[0].getFWHM()*3600/3)
        imgsizepix = 1024 # debug
        if imgsizepix%2 != 0: imgsizepix += 1 # prevent odd img sizes
    
        logger.info('Cleaning normal...')
        imagename = 'img/cal'
        lib_util.run_wsclean(s, 'wscleanA.log', MSs_concat.getStrWsclean(), name=imagename, size=imgsizepix, scale='3arcsec',
                auto_mask=5, local_rms='', local_rms_method='rms-with-min',
                weight='briggs -0.3', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.5,
                baseline_averaging='', auto_threshold=2, join_channels='', fit_spectral_pol=5, channels_out=12)

    ### DONE

logger.info("Done.")
