#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob, re
import numpy as np
from casacore import tables

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-init.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset('lilf.config')
parset_dir = parset.get('uGMRT','parset_dir')
data_dir = parset.get('uGMRT','data_dir')
bl2flag = parset.get('flag','antennas')

#############################
# 1. Download data from NCRA server

#############################
# 2. Convert into MS

#############################
# 3. Set-up directory structure

MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

if not os.path.exists('cals'): os.mkdir('cals')
if not os.path.exists('tgts'): os.mkdir('tgts')

for j, MS in enumerate(MSs.getListObj()):
    pathMS = MS.pathMS

    logger.info("- Split-up of original MS '" + pathMS + "' -")

    scanIDs = tables.taql("select distinct SCAN_NUMBER from $pathMS").getcol("SCAN_NUMBER")
    numberOfScans = len(scanIDs)

    # Create temporary paths for the sub-MSs.
    for i, scanID in enumerate(scanIDs):
        pathMSscan = "scanID%03i.MS" % scanID
        if os.path.exists(pathMSscan): continue
        logger.info("%i/%i - Creating MS: %s" % (i+1, numberOfScans, pathMSscan) )
        tables.taql("SELECT from $pathMS where SCAN_NUMBER = $scanID giving $pathMSscan as plain")

#############################
# 4. fix MS
MSs = lib_ms.AllMSs( glob.glob('scanID*MS'), s )
logger.info('Fixing MS...')
MSs.run('fixuGMRTms.py -v $pathMS', log='$nameMS_fixMS.log', commandType='python')

##################################
# 5. Move the sub-MSs from a temporary place to their right locus in the file tree.
for j, MS in enumerate(MSs.getListObj()):

        if (MS.isCalibrator()):
            pathDirectoryCalibrator = 'cals/%03i_%03i_%s' % (j, i, MS.getNameField())
            pathMSFinal             = pathDirectoryCalibrator + '/' + MS.getNameField() + ".MS"
            if not os.path.isdir(pathDirectoryCalibrator):
                os.mkdir(pathDirectoryCalibrator)
        else:
            # TODO: check in all tgts MS already created if there's a similar pointing,
            # in that case name this accordingly and combine
            pathDirectoryTarget = 'tgts/' + MS.getNameField()
            pathMSFinal         = pathDirectoryTarget + "/mss/" + MS.getNameField() + ".MS"
            if not os.path.isdir(pathDirectoryTarget):
                os.mkdir(pathDirectoryTarget)
                os.mkdir(pathDirectoryTarget+'/mss')
    
        logger.info(MS.nameMS + ": move to '" + pathMSFinal + "'...")
        MS.move(pathMSFinal)


##################################
# 6. Flagging
MSs = lib_ms.AllMSs( glob.glob('cals/*/*MS')+glob.glob('tgts/*/mss/*MS'), s )

logger.info('Flagging...')
MSs.run("DPPP " + parset_dir + "/DPPP-flag_user.parset msin=$pathMS ant.baseline=\"" + bl2flag+"\"", log="$nameMS_flag.log", commandType="DPPP")

logger.info('Done.')
