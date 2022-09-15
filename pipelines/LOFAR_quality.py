#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob, re, pickle
from astropy.io.fits import getdata
import numpy as np
from astropy.table import Table
from astropy.stats import median_absolute_deviation
import bdsf

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_cat
logger_obj = lib_log.Logger('pipeline-quality')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-quality.walker')

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
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('quality')
    os.makedirs('quality')

### DONE

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_flags=False)
ra, dec = MSs.getListObj()[0].getPhaseCentre()
fwhm = MSs.getListObj()[0].getFWHM()
qdict = {'self_c0_rms': None, 'self_c1_rms': None, 'ddcal_c0_rms': None,
                'ddcal_c1_rms': None, 'nvss_ratio': None}
# MS flags, count all flags and print %

# self images [noise per cycle]
if os.path.exists('self'):
    img_self_c0 = self_dir+'/images/wideM-0-MFS-residual.fits'
    qdict['self_c0_rms'] = get_noise(img_self_c0)
    logger.info(f'Self residual rms noise (cycle 0): %.1f mJy/b' % (qdict["self_c0_rms"]*1e3))
    img_self_c1 = self_dir+'/images/wideM-1-MFS-residual.fits'
    qdict['self_c1_rms'] = get_noise(img_self_c1)
    logger.info(f'Self residual rms noise (cycle 1): %.1f mJy/b' % (qdict["self_c1_rms"]*1e3))
else:
    logger.warning('Skip "self" tests, missing dir.')

# ddcal images [noise per cycle, astrometry, fluxscale]
if os.path.exists('ddcal'):
    img_ddcal_c0 = sorted(glob.glob(ddcal_dir+'/c00/images/wideDD-c00.residual*.fits'))[-1]
    qdict['ddcal_c0_rms'] = get_noise(img_ddcal_c0)
    logger.info('DDcal residual rms noise (cycle 0): %.1f mJy/b' % (qdict['ddcal_c0_rms']*1e3))
    img_ddcal_c1 = sorted(glob.glob(ddcal_dir+'/c01/images/wideDD-c01.residual*.fits'))[-1]
    qdict['ddcal_c1_rms'] = get_noise(img_ddcal_c1)
    logger.info('DD-cal residual rms noise (cycle 1): %.1f mJy/b' % (qdict['ddcal_c1_rms']*1e3))

    with w.if_todo('process_ddimage'):
        os.chdir(f'{ddcal_dir}/c01/images/') # bdsf raises error if image not in wdir?
        img = bdsf.process_image(f'wideDD-c01.int.restored.fits', detection_image=f'wideDD-c01.app.restored.fits',
                                 thresh_isl=4.0, thresh_pix=5.0, rms_box=(150, 15),
                                 rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True,
                                 adaptive_thresh=50, rms_box_bright=(60, 15),group_by_isl=False, group_tol=10.0,
                                 output_opts=True, output_all=True, atrous_do=False, atrous_jmax=4, flagging_opts=True,
                                 flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None)
        os.chdir('../../../')
        img.write_catalog(outfile='quality/wideDD-c01.int.cat.fits', catalog_type='srl', format='fits', correct_proj='True',
                          clobber=True)

    nvss = lib_cat.RadioCat(f'{parset_dir}/NVSS_small.fits', 'NVSS', log=logger, col_pflux=None, col_maj=None)
    nvss.filter(sigma=5, circle=[ra, dec, fwhm/2], isolation=90)
    lofar = lib_cat.RadioCat('quality/wideDD-c01.int.cat.fits', 'LOFAR', log=logger)
    lofar.filter(sigma=5, circle=[ra, dec, fwhm/2], isolation=30, minflux=0.06, size=25)
    lofar.match(nvss, 10)
    lofar.write('quality/wideDD-c01.int.cat_match_nvss.fits', overwrite=True)
    median_nvss_ratio = lofar.flux_ratio('NVSS')
    qdict['nvss_ratio'] = median_nvss_ratio
else:
    logger.warning('Skip "ddcal" tests, missing dir.')

with open('quality/quality.pickle', 'wb') as f:
    pickle.dump(qdict, f)

logger.info("Done.")
