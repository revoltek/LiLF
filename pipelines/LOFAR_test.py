#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob, re
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-test.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = '/home/baq1889/scripts/LiLF/parsets/LOFAR_cal'

#############################################################
MSs = lib_ms.AllMSs( glob.glob('*MS'), s )
calname = MSs.getListObj()[0].getNameField()
obsmode = MSs.getListObj()[0].getObsMode()

# predict to save time ms:MODEL_DATA
logger.info('Add model to MODEL_DATA...')
MSs.run("DPPP " + parset_dir + "/DPPP-predict.parset msin=$pathMS pre.sourcedb=img/calM-sources-cut.skydb", log="$nameMS_pre.log", commandType="DPPP")

logger.info('Subtract model...')
MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql2.log', commandType ='general')

logger.info("Done.")
