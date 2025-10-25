#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import os, glob
import numpy as np
import shutil

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, lib_h5
logger_obj = lib_log.Logger('pipeline-dynspec')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-dynspec.walker')

with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('mss-dynspec')
### DONE
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_consistency=True)
pixscale = MSs.getListObj()[0].getPixelScale()
imgsizepix = int(1.85*max(MSs.getListObj()[0].getFWHM(freq='max', elliptical=True)) * 3600 / pixscale) # roughly to smallest null
if imgsizepix > 10000: imgsizepix = 10000 # keep SPARSE doable
if imgsizepix % 2 != 0: imgsizepix += 1  # prevent odd img sizes

MSs_avg = lib_ms.AllMSs( glob.glob('mss-avg/TC*[0-9].MS'), s)
fwhm = max(MSs_avg.getListObj()[0].getFWHM(freq='mid', elliptical=True))
workingReg = 'ddserial/workingRegion.reg' # sources outside of this region will be ignored (and not peeled)
MSs_avg.getListObj()[0].makeBeamReg(workingReg, freq='min', to_pbval=0)
peelReg = 'ddserial/peelingRegion.reg' # sources outside of this region will be peeled
MSs_avg.getListObj()[0].makeBeamReg(peelReg, freq='max', to_pbval=0.12) # this is slighly smaller than the null
freq_min = np.min(MSs_avg.getFreqs())
freq_mid = np.mean(MSs_avg.getFreqs())
phase_center = MSs_avg.getListObj()[0].getPhaseCentre()
timeint = MSs_avg.getListObj()[0].getTimeInt()
ch_out = MSs_avg.getChout(4e6)  # for full band (48e6 MHz) is 12
mode = MSs.getListObj()[0].getAntennaSet()

MSs_lres = lib_ms.AllMSs( glob.glob('mss-lres/TC*[0-9].MS'), s)

parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_ddserial'])))
parset_dir = parset.get('LOFAR_ddserial','parset_dir')
maxIter = parset.getint('LOFAR_ddserial','maxIter')
min_cal_flux60 = parset.getfloat('LOFAR_ddserial','minCalFlux60') # default: 0.8
solve_amp = parset.getboolean('LOFAR_ddserial','solve_amp')
manual_dd_cal = parset.get('LOFAR_ddserial','manual_dd_cal') # ds9 circle region file containing a manual dd-calibrator
develop = parset.getboolean('LOFAR_ddserial', 'develop') # for development, make more output/images
use_shm = parset.getboolean('LOFAR_ddserial', 'use_shm') # use shared memory for wsclean
userReg = parset.get('model','userReg')

with w.if_todo('predict'):
        # wsclean predict - from ddparallel in cycle 0, otherwise from previous iteration
        logger.info('Predict full model...')
        s.add(f'wsclean -predict -padding 1.8 -name ddserial/c00/images/wideDD-c00 -j {s.max_cpucores} -channels-out {ch_out} \
                -facet-regions ddserial/c00/solutions/facets-c00.reg -apply-facet-solutions ddserial/c00/solutions/interp.h5 phase000,amplitude000 \
                -apply-facet-beam -use-differential-lofar-beam -facet-beam-update 120 -no-solution-directions-check \
                -reorder -parallel-reordering 4 {MSs_avg.getStrWsclean()}',
              log='wscleanPRE-c' + str('00') + '.log', commandType='wsclean')
        s.run(check=True)

#Leakage calibration
with w.if_todo('corr-leakage'):
       logger.info('Correct amp-di (fulljones)...')
       MSs_avg.run(f'DP3 {parset_dir}/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
               cor.parmdb=ddserial/c00/solutions/cal-leak.h5 cor.correction=fulljones cor.soltab=[amplitude000,phase000] \
               cor.updateweights=False', log='$nameMS_leakcorr.log', commandType='DP3')

##############################################################################################################
### Calibration finished - additional images with scientific value

# Low res as this is relevant only for transient detection
### TODO dynspec part

#with w.if_todo('output-lressub'):
#    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
#    MSs_avg.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
#        log='$nameMS_taql.log', commandType='general')
### DONE
with w.if_todo('make_catalog'):
    os.makedirs('dynspec/', exist_ok=True)
    logger.info('Making catalog...')
    s.add(f'LOFAR_FoV_dynspec.py \
  --catalog /project/lspc/Data/edler/Long_project/mss/LiLF/models/Dynspec_catalog.fits \
  --beam ddserial/primarybeam.fits \
  --quality-catalog quality/wideDD-c00-MFS-image-pb.cat.fits \
  --mask-threshold 0.3 \
  --n-offbeam 100 \
  --min-sep-arcmin 1.0 \
  --outdir dynspec \
  -v', log='target.log', commandType='general')
  	
    s.run(check=True)

# HE: that should not be required anymore in newest DP3
# TDODO HE also: add h5parm correction and beam correction once DP3 feature is ready.
with w.if_todo('make_dynspec'):
    MSs_avg.run(f'DP3 {parset_dir}/DP3-dynspec.parset msin=$pathMS dynspec.sourcelist=dynspec/TARGET_LIST.skymodel dynspec.fitsprefix=dynspec/$nameMS_', log='dynspec.log', commandType='DP3')
    s.run(check=True)    
    shutil.rmtree('dynspec-out.ms')
w.alldone()
