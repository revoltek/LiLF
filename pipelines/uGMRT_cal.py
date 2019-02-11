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
msIDs = np.unique([ int(name[8:10]) for name in glob.glob('cals/*') ])
for msID in msIDs:

    logger.info('Working on MSid %02i' % msID)

    MSs = lib_ms.AllMSs( glob.glob('cals/%02i*/*MS' % msID), s )
    calname = MSs.getListObj()[0].getNameField()
    for MS in MSs.getListObj():
        os.system('cp -r %s %s' % (skymodel, MS.pathMS))
    
    # predict to save time ms:MODEL_DATA
    logger.info('Add model to MODEL_DATA (%s)...' % calname)
    MSs.run("DPPP " + parset_dir + "/DPPP-predict.parset msin=$pathMS pre.sourcedb=$pathMS/" + os.path.basename(skymodel) + " pre.sources=" + calname, \
            log="$nameMS_pre.log", commandType="DPPP")
    
    # Smooth data DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType ='python', maxThreads=10)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/cal.h5 sol.mode=rotation+diagonal', log='$nameMS_sol.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'cal', [ms+'/cal.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset'])

logger.info("Done.")
