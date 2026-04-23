#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation.
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import sys, os, glob, re
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-cal.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('uGMRT_cal','parset_dir')
skymodel = parset.get('uGMRT_cal','skymodel')
bl2flag = parset.get('flag','stations')

#############################################################
# cicle on each msID assuming telescope properties are similar in all calibrator scans of the same msID
msIDs = np.unique([ int(name[5:7]) for name in glob.glob('cals/*') ])
for msID in msIDs:

    logger.info('Working on MSid %02i' % msID)

    MSs_cals = lib_ms.AllMSs( glob.glob('cals/%02i*/*MS' % msID), s )
    calnames = [MS.getNameField() for MS in MSs_cals.getListObj()]
    for MS in MSs_cals.getListObj():
        os.system('cp -r %s %s' % (skymodel, MS.pathMS))
    
    logger.info('Add columns...')
    MSs_cals.run('addcol2ms.py -m $pathMS -c MODEL_DATA -i DATA', log="$nameMS_addcol.log", commandType="python")

    # predict to save time ms:MODEL_DATA
    for i, MS in enumerate(MSs_cals.getListObj()):
        logger.info('%s: add model to MODEL_DATA (%s)...' % (MS.pathMS, calnames[i]) )
        s.add("DP3 " + parset_dir + "/DP3-predict.parset msin="+MS.pathMS+" pre.sourcedb="+MS.pathMS+"/" + os.path.basename(skymodel) + " pre.sources=" + calnames[i], \
            log=MS.nameMS+"_pre.log", commandType="DP3")
    s.run(check=True)
    
    # Smooth data DATA -> SMOOTHED_DATA (BL-based smoothing)
    # TODO: gmrt amplitude are unstable, the smooth seems not to work fine
    #logger.info('BL-smooth...')
    #MSs_cals.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType ='python', maxThreads=10)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    MSs_cals.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/diag.h5 sol.mode=diagonal', log='$nameMS_sol.log', commandType="DP3")
    
    lib_util.run_losoto(s, 'diag', [ms+'/diag.h5' for ms in MSs_cals.getListStr()], \
            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-flag.parset', \
            parset_dir+'/losoto-pa.parset', parset_dir+'/losoto-iono.parset', parset_dir+'/losoto-bp.parset'])

    # Transfer to target
    MSs_tgts = lib_ms.AllMSs( glob.glob('tgts/*/mss/%02i*MS' % msID), s )

    logger.info('Add columns...')
    MSs_tgts.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA -i DATA', log="$nameMS_addcol.log", commandType="python")

    # Trasnfer tgt_SB.MS:DATA -> tgt.MS:CORRECTED_DATA
    logger.info('Correcting targets...')
    MSs_tgts.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS', log='$nameMS_cor.log', commandType="DP3")

    # Combine and split
    # TODO: find a way to concat in time
    #for tgt_dir in glob.glob('tgts/*'):
    #    logger.info('Splitting %s...' % tgt_dir)
    #    mss_in = glob.glob(tgt_dir+'/mss/*MS')
    #    ms_out = tgt_dir+'/mss/concat_'+tgt_dir+'.MS'
    #    log = tgt_dir.split('/')[1]+'_split.log'
    #    s.add('DP3 ' + parset_dir + '/DP3-split.parset msin="['+','.join(mss_in)+']" msin.datacolumn=CORRECTED_DATA msout='+ms_out, log=log, commandType="DP3")
    #s.run(check=True)
    
    # TODO: for now just use the un-splitted files
    logger.info('Set DATA = CORRECTED_DATA...')
    MSs_tgts.run('taql "update $pathMS set DATA = CORRECTED_DATA"', log='$nameMS_taql1.log', commandType='general')

logger.info("Done.")
