#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_ddfacet
logger_obj = lib_log.Logger('pipeline-ddfacet.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)#, maxThreads = 4)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_ddfacet','parset_dir')
maxniter = parset.getint('LOFAR_ddfacet','maxniter')
calFlux = parset.getfloat('LOFAR_ddfacet','calFlux')
userReg = parset.get('model','userReg')

####################################################
MSs_self = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs_self.getListObj()[0].getPhaseCentre()
fwhm = MSs_self.getListObj()[0].getFWHM(freq='mid')

############################
logger.info('Cleaning...')
lib_util.check_rm('ddfcal')
os.makedirs('ddfcal/masks')
os.makedirs('ddfcal/plots')
os.makedirs('ddfcal/images')
os.makedirs('ddfcal/solutions')
os.makedirs('ddfcal/skymodels')

############################################################
# use SUBTRACTED_DATA (no pre-correction - subtraction would not work) or CORRECTED_DATA (DIE iono correction)?
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs_self.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
MSs.getListObj()[0].makeBeamReg('ddcal/beam.reg', freq='mid')
beamReg = 'ddcal/beam.reg'
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg)
mosaic_image.selectCC()

rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    # calibrate
    lib_ddfacet.killms_data()
    
    # image
    lib_ddfacet.ddf_image()
