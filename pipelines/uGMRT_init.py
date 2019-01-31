#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob, re
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-init.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('init','parset_dir')
data_dir = parset.get('init','data_dir')
bl2flag = parset.get('flag','stations')

#############################
# 1. Download data from NCRA server

#############################
# 2. Convert into MS

#############################
# 3. Set-up directory structure

MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

os.mkdir('cals')
os.mkdir('tgts')

for MS in MS.getObjList():
	pathMS = MS.pathMS

	logger.info("- Split-up of original MS '" + pathMS + "' -")

	scanIDs = tables.taql("select distinct SCAN_NUMBER from $pathMS").getcol("SCAN_NUMBER")
	numberOfScans = len(scanIDs)

# Create temporary paths for the sub-MSs.
for i, scanID in enumerate(scanIDs):
    pathMSNew = "scanID%03i.MS" % scanID
    tables.taql("SELECT from $pathMS where SCAN_NUMBER = $scanID giving $pathMSNew as plain")
    logger.info("%i/%i - Created MS: %s" % (i, numberOfScans, pathMSNew) )

##################################
# 6. Move the sub-MSs from a temporary place to their right locus in the file tree.
for pathMSNew in pathsMSNew:

    MSObject    = lib_ms.MS(pathMSNew)

    if (MSObject.isCalibrator()):
        pathDirectoryCalibrator = pathDirectoryFieldsCalibrator + '/' + MSObject.nameMS
        pathMSFinal             = pathDirectoryCalibrator + '/' + MSObject.nameMS + ".MS"
        if (not os.path.isdir(pathDirectoryCalibrator)):
            os.mkdir(pathDirectoryCalibrator)
            os.mkdir(pathDirectoryCalibrator + "/plots")
            os.mkdir(pathDirectoryCalibrator + "/solutions")
    else:
        pathDirectoryTarget = pathDirectoryFieldsTarget + '/' + MSObject.getNameField()
        pathMSFinal         = pathDirectoryTarget + "/MSs/" + MSObject.nameMS + ".MS"
        if (not os.path.isdir(pathDirectoryTarget)):
            os.mkdir(pathDirectoryTarget)
            os.mkdir(pathDirectoryTarget + "/plots")
            os.mkdir(pathDirectoryTarget + "/MSs")
            os.mkdir(pathDirectoryTarget + "/images")

    logging.info("Moving MS at '" + pathMSNew + "' to '" + pathMSFinal + "'...")

    MSObject.move(pathMSFinal)


#############################
# 5. fix MS

MSs.run()

logger.info('Done.')
