#!/usr/bin/python
# demix of a set of SBs from a given dir, output is in the local dir

import sys, os, glob
import numpy as np

###############################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-demix.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_demix','parset_dir')
data_dir = parset.get('LOFAR_demix','data_dir')
skydb = parset.get('LOFAR_demix','demix_model')

##############################################
# Demix
MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

for MS in MSs.getListStr():
    lib_util.check_rm(os.path.basename(MS)+'_'+os.path.basename(skydb))
    print('cp -r '+skydb+' '+os.path.basename(MS)+'_'+os.path.basename(skydb))
    os.system('cp -r '+skydb+' '+os.path.basename(MS)+'_'+os.path.basename(skydb))

logger.info('Demixing...')
MSs.run('DPPP '+parset_dir+'/DPPP_demix.parset msin=$pathMS msout=$nameMS demixer.skymodel=$nameMS.MS_'+os.path.basename(skydb)+' demixer.instrumentmodel=$nameMS/instrument_demix', \
        log='$nameMS_demix.log', commandType='DPPP')

logger.info("Done.")
