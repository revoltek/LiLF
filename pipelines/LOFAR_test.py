#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob, re
import pyrap.tables as pt
from astropy.time import Time
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-test.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
#parset = lib_util.getParset()
parset_dir = '/home/baq1889/scripts/LiLF/parsets/LOFAR_dd'

#############################################################
MSs = lib_ms.AllMSs( glob.glob('*MS'), s, check_flags=False )

skymodel_voro_skydb = 'skymodel00_voro.skydb'
c=0

# Calibration - ms:SMOOTHED_DATA
logger.info('Calibrating...')
MSs.run('DPPP '+parset_dir+'/DPPP-solDD.parset msin=$pathMS ddecal.h5parm=$pathMS/cal-core.h5 \
        ddecal.sourcedb='+skymodel_cl_skydb+' ddecal.solint=15 ddecal.nchan=30', \
        log='$nameMS_solDD-core.log', commandType='DPPP')

lib_util.run_losoto(s, 'core', [MS+'/cal-core.h5' for MS in MSs.getListStr()], \
        [parset_dir+'/losoto-core.parset'])

logger.info('re-calibration...')
# predict and corrupt each facet
logger.info('Reset MODEL_DATA...')
MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

for i, d in enumerate(directions):
    # predict - ms:MODEL_DATA
    logger.info('Patch '+d.name+': predict...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_DIR pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
        log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

    # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
    logger.info('Patch '+d.name+': corrupt...')
    MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS msin.datacolumn=MODEL_DATA_DIR msout.datacolumn=MODEL_DATA_DIR cor.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor.direction=['+d.name+']', \
        log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

    logger.info('Patch '+d.name+': subtract...')
    MSs.run('taql "update $pathMS set MODEL_DATA = MODEL_DATA + MODEL_DATA_DIR"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')



logger.info("Done.")
