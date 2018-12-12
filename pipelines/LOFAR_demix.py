#!/usr/bin/python
# demix of a set of SBs from a given dir, output is in the local dir

import sys, os, glob
import numpy as np
from autocal.lib_pipeline import *

###############################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-demix.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False, max_threads = 4) # set here max number of threads here

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('demix','parset_dir')
data_dir = parset.get('demix','data_dir')

##############################################
# Demix
MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )
logger.info('Demixing...')
MSs.run('DPPP '+parset_dir+'/DPPP-demix.parset msin=$pathMS msout=$nameMS demixer.instrumentmodel=$nameMS/instrument_demix', log='$nameMS_demix.log', commandType='DPPP')

logger.info("Done.")
