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
#parset_dir = '/home/baq1889/scripts/LiLF/parsets/LOFAR_cal'

#############################################################
MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

with pt.table(MSs.getListStr()[0]+'/OBSERVATION', readonly=True, ack=False) as obs:
    t = Time(obs.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
    time = np.int(t.iso.replace('-','')[0:8])

# Rescale visibilities by 1e3 if before 2014-03-19 (old correlator), and by 1e-2 otherwise
if time < 20140500:
    MSs.run('taql "update $pathMS set DATA = 1e10*DATA"', log='$nameMS_taql.log', commandType='general')
else:
    MSs.run('taql "update $pathMS set DATA = 1e-4*DATA"', log='$nameMS_taql.log', commandType='general')

logger.info("Done.")
