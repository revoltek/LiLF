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
parset_dir = parset.get('uGMRT_init','parset_dir')
data_dir = parset.get('uGMRT_init','data_dir')
bl2flag = parset.get('flag','antennas')

#############################
# 1. Download data from NCRA server

#############################
# 2. Convert into MS

#############################
# 3. Set-up directory structure

MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )
numberOfMSs = len(MSs.getListObj())

if not os.path.exists('cals'): os.mkdir('cals')
if not os.path.exists('tgts'): os.mkdir('tgts')

for msID, MS in enumerate(MSs.getListObj()):
    pathMS = MS.pathMS

    logger.info("- Split-up of original MS '%s' (msID: %02i) -" % (pathMS, msID))

    scanIDs = tables.taql("select distinct SCAN_NUMBER from $pathMS").getcol("SCAN_NUMBER")
    numberOfScans = len(scanIDs)

    # Create temporary paths for the sub-MSs.
    for i, scanID in enumerate(scanIDs):
        pathMSscan = "msID%02i_scanID%03i.MS" % (msID, scanID)
        if os.path.exists(pathMSscan): continue
        logger.info("MS %02i/%02i: %03i/%03i - Creating MS: %s" % (msID+1, numberOfMSs, i+1, numberOfScans, pathMSscan) )
        tables.taql("SELECT from $pathMS where SCAN_NUMBER = $scanID giving $pathMSscan as plain")

#############################
# 4. fix MS
MSs = lib_ms.AllMSs( glob.glob('msID*_scanID*MS'), s )
logger.info('Fixing MS...')
MSs.run('fixuGMRTms.py -v $pathMS', log='$nameMS_fixMS.log', commandType='python')

##################################
# 5. Move the sub-MSs from a temporary place to their right locus in the file tree.
for MS in MSs.getListObj():

    msID = int(MS.nameMS[4:6])
    scanID = int(MS.nameMS[14:17])

    if (MS.isCalibrator()):
        pathDirectoryCalibrator = 'cals/%02i_%03i_%s' % (msID, scanID, MS.getNameField())
        pathMSFinal             = pathDirectoryCalibrator + '/' + MS.getNameField() + ".MS"
        if not os.path.isdir(pathDirectoryCalibrator):
            os.mkdir(pathDirectoryCalibrator)
    else:
        # TODO: check in all tgts MS already created if there's a similar pointing,
        # in that case name this accordingly and combine
        pathDirectoryTarget = 'tgts/' + MS.getNameField()
        pathMSFinal         = pathDirectoryTarget + "/mss/%02i_%03i_%s.MS" % (msID, scanID, MS.getNameField())
        if not os.path.isdir(pathDirectoryTarget):
            os.mkdir(pathDirectoryTarget)
            os.mkdir(pathDirectoryTarget+'/mss')

    logger.info(MS.nameMS + ": move to '" + pathMSFinal + "'...")
    MS.move(pathMSFinal)


##################################
# 6. Flagging
MSs = lib_ms.AllMSs( glob.glob('cals/*/*MS')+glob.glob('tgts/*/mss/*MS'), s )

logger.info('Flagging...')
MSs.run("DPPP " + parset_dir + "/DPPP-flag.parset msin=$pathMS flagBaselines.baseline=\"" + bl2flag + "\" flagRFI.strategy=" + parset_dir + "/uGMRTDummy.rfis", \
        log="$nameMS_flag.log", commandType="DPPP")

logger.info('Done.')
