#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob, re
import pyrap.tables as pt
from astropy.time import Time
from astropy.io.fits import getdata
import numpy as np
from astropy.stats import median_absolute_deviation

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

def get_noise(fitsfile):
    eps = 1e-3; niter = 100; sigma = 5
    img_data = getdata(fitsfile)
    data = img_data[ ~np.isnan(img_data) & (img_data != 0) ] # remove nans and 0s
    initial_len = len(data)
    if initial_len == 0: return 0
    mad_old = 0.
    for i in range(niter):
         mad = median_absolute_deviation(data)
         logger.debug('%s: MAD noise: %.1f mJy on %f%% data' % (fitsfile, mad*1e3, 100*len(data)/initial_len))
         if np.isnan(mad): return 0
         if np.abs(mad_old-mad)/mad < eps:
             rms = np.nanstd( data )
             logger.debug('%s: Noise: %.1f mJy/b (data len: %i -> %i - %.2f%%)' % (fitsfile, rms*1e3, initial_len, len(data), 100*len(data)/initial_len))
             return rms

         data = data[np.abs(data) < (sigma*mad)]
         mad_old = mad


#############################################################
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_flags=False)

# MS flags, count all flags and print %

# self images [noise per cycle]
if os.path.exists('self'):
    img_self_c0 = self_dir+'/images/wideM-0-MFS-residual.fits'
    logger.info('Self residual rms noise (cycle 0): %.1f mJy/b' % (get_noise(img_self_c0)*1e3))
    img_self_c1 = self_dir+'/images/wideM-1-MFS-residual.fits'
    logger.info('Self residual rms noise (cycle 1): %.1f mJy/b' % (get_noise(img_self_c1)*1e3))
else:
    logger.warning('Skip "self" tests, missing dir.')

# ddcal images [noise per cycle, astrometry, fluxscale]
if os.path.exists('ddcal'):
    img_ddcal_c0 = sorted(glob.glob(ddcal_dir+'/c00/images/wideDD-c00.residual*.fits'))[-1]
    logger.info('DDcal residual rms noise (cycle 0): %.1f mJy/b' % (get_noise(img_ddcal_c0)*1e3))
    img_ddcal_c1 = sorted(glob.glob(ddcal_dir+'/c01/images/wideDD-c01.residual*.fits'))[-1]
    logger.info('DDcal residual rms noise (cycle 1): %.1f mJy/b' % (get_noise(img_ddcal_c1)*1e3))
else:
    logger.warning('Skip "ddcal" tests, missing dir.')

logger.info("Done.")
