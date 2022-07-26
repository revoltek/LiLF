#!/usr/bin/env python3
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
ddcal_dir = parset.get('LOFAR_quality','ddcal_dir')

#############################################################
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_flag=False)

# MS flags, count all flags and print %

# self images [noise per cycle]
if os.path.exists('self'):
    img_self_c0 = lib_img.Image(self_dir+'/images/wideM-0-MFS-residual.fits')
    logging.info('Self residual rms noise (cycle 0): %f mJy/b' % (imag_self_c0.getNoise()*1e-3))
    img_self_c1 = lib_img.Image(self_dir+'/images/wideM-1-MFS-residual.fits')
    logging.info('Self residual rms noise (cycle 1): %f mJy/b' % (imag_self_c1.getNoise()*1e-3))
else:
    logging.warning('Skip "self" tests, missing dir.')

# ddcal images [noise per cycle, astrometry, fluxscale]
if os.path.exists('ddcal'):
    img_ddcal_c0 = lib_img.Image(sorted(glob.glob(ddcal_dir+'/c00/images/wideDD-c00.residual*.fits'))[-1])
    logging.info('DDcal residual rms noise (cycle 0): %f mJy/b' % (imag_ddcal_c0.getNoise()*1e-3))
    img_ddcal_c1 = lib_img.Image(sorted(glob.glob(ddcal_dir+'/c01/images/wideDD-c01.residual*.fits'))[-1])
    logging.info('DDcal residual rms noise (cycle 1): %f mJy/b' % (imag_ddcal_c1.getNoise()*1e-3))
else:
    logging.warning('Skip "ddcal" tests, missing dir.')

logger.info("Done.")
