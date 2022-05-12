#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob, re
import pyrap.tables as pt
from astropy.time import Time
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-quality.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_quality','parset_dir')
self_dir = parset.get('LOFAR_quality','self_dir')
dd_dir = parset.get('LOFAR_quality','dd_dir')

#############################################################
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_flag=False)

# MS flags

# self images [noise per cycle]

# dd images [noise per cycle, astrometry, fluxscale]


logger.info("Done.")
