#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, glob, pickle
from astropy.io.fits import getdata
import numpy as np
from casacore.tables import table
from astropy.stats import median_absolute_deviation
import bdsf

########################################################
from LiLF import lib_ms, lib_util, lib_log, lib_cat, lib_img
logger_obj = lib_log.Logger('pipeline-quality')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-quality.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_quality','parset_dir')
ddparallel_dir = parset.get('LOFAR_quality','ddparallel_dir')
ddserial_dir = parset.get('LOFAR_quality','ddserial_dir')

#############################################################
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('quality')
    os.makedirs('quality')

### DONE

# get the target name if available
if os.path.exists('../target.txt'):
    with open("../target.txt", "r") as file:
        target = file.readline()[:-1]
else:
    target = None

MSs = lib_ms.AllMSs( glob.glob('mss-avg/TC*[0-9].MS'), s, check_flags=False)
ra, dec = MSs.getListObj()[0].getPhaseCentre()
fwhm = MSs.getListObj()[0].getFWHM(elliptical=True)
qdict = {'ddparallel_c0_rms': None, 'ddparallel_c1_rms': None, 'ddserial_c0_rms': None,
                'ddserial_c1_rms': None, 'nvss_ratio': None, 'nvss_match': None, 'flag_frac':None}
# MS flags, count all flags and print %

# ddparallel images [noise per cycle]
if os.path.exists(ddparallel_dir):
    img = lib_img.Image(ddparallel_dir+'/images/wideDDP-c0-MFS-residual.fits')
    qdict['ddparallel_c0_rms'] = img.getNoise(useMask=False)
    logger.info('ddparallel residual rms noise (cycle 0): %.1f mJy/b' % (qdict["ddparallel_c0_rms"]*1e3))
    img = lib_img.Image(ddparallel_dir+'/images/wideDDP-c1-MFS-residual.fits')
    qdict['ddparallel_c1_rms'] = img.getNoise(useMask=False)
    logger.info('ddparallel residual rms noise (cycle 1): %.1f mJy/b' % (qdict["ddparallel_c1_rms"]*1e3))
else:
    logger.warning('Skip "ddparallel" tests, missing dir.')

# flag fraction in blocks of 10 min
flags = []
weights = []
for MS in MSs.getListObj():
    t_start, t_end = MS.getTimeRange()
    t_diff = (t_end - t_start)/6
    t = table(MS.pathMS, ack=False)
    for i in range(6): # ASSUME 60 MIN time-split MSs and query in each 10 min interval
        f = t.query(query=f'(TIME>{t_start+i*t_diff}) AND ({t_start+(i+1)*t_diff}>TIME)', columns='FLAG').getcol('FLAG')
        this_flag_frac = np.mean(f)
        if this_flag_frac < 0.2:
            logger.info(f'{MS.nameMS}: {this_flag_frac:.1%} flagged in interval {i}')
        else:
            logger.warning(f'{MS.nameMS}: {this_flag_frac:.1%} flagged in interval {i}')
        flags.append(this_flag_frac)
        weights.append(t_diff)

flag_frac = np.average(flags, weights=weights) # weighting just in case the MSs are not the same length
if flag_frac < 0.2:
    logger.info(f'Total flagged: {flag_frac:.1%}')
else:
    logger.warning(f'Total flagged: {flag_frac:.1%}')
qdict['flag_frac'] = flag_frac

# ddserial images [noise per cycle, astrometry, fluxscale]
if os.path.exists(ddserial_dir):
    img = lib_img.Image(ddserial_dir+'/c00/images/wideDDS-c0-MFS-residual.fits')
    qdict['ddserial_c0_rms'] = img.getNoise(useMask=False)
    logger.info('ddserial residual rms noise (cycle 0): %.1f mJy/b' % (qdict['ddserial_c0_rms']*1e3))
    img = lib_img.Image(ddserial_dir+'/c00/images/wideDDS-c0-MFS-image.fits')
    img.plotimage('quality/wideDDS-c0-MFS-image.png', regionfile=ddserial_dir+'/peelingRegion.reg', minmax=(-5, 100), title=target)
    img.plotimage('quality/wideDDS-c0-MFS-image-reg.png', regionfile=ddserial_dir+'/c00/solutions/facetsS-c0.reg', minmax=(-5, 100), title=target)
    img = lib_img.Image(ddserial_dir+'/c00/images/wideDDS-lres-c0-MFS-image.fits')
    img.plotimage('quality/wideDDS-lres-c0-MFS-image.png', regionfile=ddserial_dir+'/peelingRegion.reg', minmax=(-5, 100), title=target+ ' (lres)')

    with w.if_todo('process_ddimage'):
        os.chdir(f'{ddserial_dir}/c00/images/') # bdsf raises error if image not in wdir?
        img = bdsf.process_image('wideDDS-c0-MFS-image-pb.fits', detection_image='wideDDS-c0-MFS-image.fits',
                                 thresh_isl=3.0, thresh_pix=5.0, rms_box=(150, 15),
                                 rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True,
                                 adaptive_thresh=50, rms_box_bright=(60, 15),group_by_isl=False, group_tol=10.0,
                                 output_opts=True, output_all=False, atrous_do=False, atrous_jmax=4, flagging_opts=True,
                                 flag_maxsize_fwhm=0.4, advanced_opts=True, blank_limit=None, quiet=True, debug=False)
        os.chdir('../../../')
        img.write_catalog(outfile='quality/wideDDS-c0-MFS-image-pb.cat.fits', catalog_type='srl', format='fits', correct_proj='True',
                          clobber=True)
        
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(f'{ddserial_dir}/c00/images/wideDDS-c0-MFS-image-pb.fits') as hdul:
        wcs = WCS(hdul[0].header)
        wcs = wcs.dropaxis(2); wcs = wcs.dropaxis(2) # remove 3rd and 4th axis

    nvss = lib_cat.RadioCat(f'{parset_dir}/NVSS_small.fits', 'NVSS', log=logger, col_pflux=None, col_maj=None, wcs=wcs)
    nvss.filter(sigma=5, ellipse=[ra, dec, fwhm[0]/2, fwhm[1]/2, 0], isolation=120)
    nvss.write('quality/debug_nvss.fits', overwrite=True, format='fits')
    lofar = lib_cat.RadioCat('quality/wideDDS-c0-MFS-image-pb.cat.fits', 'LOFAR', log=logger, wcs=wcs)
    lofar.filter(sigma=5, ellipse=[ra, dec, fwhm[0]/2, fwhm[1]/2, 0], isolation=45, minflux=0.06, size=25)
    lofar.write('quality/debug_lofar.fits', overwrite=True, format='fits')
    lofar.match(nvss, 10)
    n_match = len(lofar.get_matches('NVSS'))
    lofar.write('quality/wideDDS-c0-MFS-image-pb.cat_match_nvss.fits', overwrite=True, format='fits')
    median_nvss_ratio = lofar.flux_ratio('NVSS')
    qdict['nvss_ratio'] = median_nvss_ratio
    qdict['nvss_match'] = n_match
else:
    logger.warning('Skip "ddserial" tests, missing dir.')


with open('quality/quality.pickle', 'wb') as f:
    pickle.dump(qdict, f)

w.alldone()