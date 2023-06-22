#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re
import numpy as np
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_ddfacet
from LiLF.make_mask import make_mask
logger_obj = lib_log.Logger('pipeline-ddfacet.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)#, maxThreads = 4)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_ddfacet','parset_dir')
maxniter = parset.getint('LOFAR_ddfacet','maxniter')
calFlux = parset.getfloat('LOFAR_ddfacet','calFlux')
userReg = parset.get('model','userReg')
uvrange=[0.03,1e6]
image_robust=-0.15
ncluster = 5

####################################################
MSs_self = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs_self.getListObj()[0].getPhaseCentre()
fwhm = MSs_self.getListObj()[0].getFWHM(freq='mid')
nchan = MSs_self.getListObj()[0].getNchan()

############################
logger.info('Cleaning...')
lib_util.check_rm('ddfcal')
os.makedirs('ddfcal/masks')
os.makedirs('ddfcal/plots')
os.makedirs('ddfcal/images')
os.makedirs('ddfcal/solutions')
os.makedirs('ddfcal/cache')
lib_util.check_rm('img')
os.makedirs('img')

# Clear the shared memory
s.add('ClearSHM.py', log='clearSHM.log', commandType='singularity', processors='max')
s.run(check=True)

############################################################
# use SUBTRACTED_DATA (no pre-correction - subtraction would not work) or CORRECTED_DATA (DIE iono correction)?
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    for n in range(1,nchan%122):
        MSs_self.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msout=mss-dd/$nameMS-n%i.MS msin.datacolumn=CORRECTED_DATA msin.nchan=122 msin.startchan=%i \
                avg.freqstep=1 avg.timestep=1' % (n, n*122), log='$nameMS_avg.log', commandType='DP3')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9]-n[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
MSs.getListObj()[0].makeBeamReg('ddfcal/beam.reg', freq='mid')
beamReg = 'ddfcal/beam.reg'
rms_noise_pre = np.inf

# DIE clean
logger.info('DIE cleaning...')
rootimg1 = 'img/image_dirin_SSD'
lib_ddfacet.ddf_image(s, 'ddf1.log', MSs, rootimg1, cleanmode='SSD', \
                      majorcycles=1, robust=image_robust, reuse_psf=False, reuse_dirty=False,
                      peakfactor=0.02, rms_factor=3, colname='DATA', clusterfile=None, \
                      automask=True, automask_threshold=15, \
                      apply_weights=False, use_weightspectrum=True, uvrange=uvrange)

sys.exit()

# make mask
logger.info('Making mask...')
maskname = rootimg1+'.app.mask.fits'
make_mask(rootimg1+'.app.restored.fits', mask_name=maskname, threshisl=5, atrous_do=False, rmsbox=(50,20))


# DIE cleaning (continue)
logger.info('DIE cleaning (2)...')
rootimg2 = 'img/image_dirin_SSD'
lib_ddfacet.ddf_image(s, 'ddf2.log', MSs, rootimg2, cleanmode='SSD', \
                      majorcycles=1, robust=image_robust, reuse_psf=True, reuse_dirty=True, \
                      peakfactor=0.01, rms_factor=3, colname='DATA', clusterfile=None,
                      automask=True, automask_threshold=15, cleanmask=maskname, \
                      apply_weights=False, use_weightspectrum=True, uvrange=uvrange,
                    RMSFactorInitHMP=1., MaxMinorIterInitHMP=10000, PredictSettings=("Clean","DD_PREDICT"))

# make source catalogue
logger.info('Making source catalogue...')
maskname = rootimg2+'.app.mask.fits'
make_mask(rootimg2+'.app.restored.fits', mask_name=maskname, threshisl=5, atrous_do=False, rmsbox=(30,10), write_srl=True)

logger.info('Grouping...')
s.add("ClusterCat.py --SourceCat %s.app.restored.pybdsm.srl.fits --DoPlot 0 --NGen 100 --NCluster %i --NCPU %i" % (rootimg2,ncluster,64), \
        log='clustercat.log', commandType='singularity', processors='max')
s.run(check=True)

clusterfile = 'img/image_dirin_SSD.app.restored.pybdsm.srl.fits.ClusterCat.npy'

## Incorporate direction with DDF
logger.info('DIE cleaning (3)...')
rootimg3 = 'img/image_dirin_SSD_c'
lib_ddfacet.ddf_image(s, 'ddf3.log', MSs, rootimg3, cleanmode='SSD', \
                                           majorcycles=1, robust=image_robust, \
                                           peakfactor=0.001, rms_factor=0, colname='DATA', clusterfile=clusterfile, \
                                           automask=True, automask_threshold=15, cleanmask=maskname, \
                                           apply_weights=False, use_weightspectrum=True, uvrange=uvrange,
                                           use_dicomodel=True, dicomodel_base=rootimg2,
                                           RMSFactorInitHMP=1., MaxMinorIterInitHMP=10000, PredictSettings=("Clean","DD_PREDICT"))

rootimg_old = rootimg3

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    solname = 'DDS%02i' % c
    rootimg = 'img/image-c%02i_dde_SSD' % c

    # calibrate
    if c == 0: dt = 2
    else: dt = 1
    logger.info('Calibrating...')
    CurrentDDkMSSolName = lib_ddfacet.killms_data(s, 'killms%02i.log' % c, MSs, rootimg_old, solname, dicomodel='%s.DicoModel'%rootimg_old,
                                        cache_dir = 'ddfcal/cache', sols_dir = 'ddfcal/solutions',
                                        clusterfile=clusterfile, CovQ=0.02, niterkf=6,
                                        uvrange=uvrange, wtuv=None, robust=image_robust, dt=dt, NChanSols=12, MergeSmooth=True)

    # run losoto
    lib_util.run_losoto(s, 'g-c'+str(c), glob.glob('ddfcal/solutions/TC*.MS/*npz'), \
                        [parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset'])
    os.system('mv plots-g-c'+str(c)+' ddfcal/plots')
    os.system('mv cal-g-c'+str(c)+'.h5 ddfcal/solutions')

    # image
    if c <= 1: applysols = 'P'
    else: applysols = 'P'

    logger.info('Imaging...')
    lib_ddfacet.ddf_image(s, 'ddf%02i.log' % c, MSs, rootimg,cleanmode='SSD',\
                                           ddsols=solname,applysols=applysols,\
                                           majorcycles=2,robust=image_robust,\
                                           peakfactor=0.001,colname='DATA',\
                                           automask=True,automask_threshold=10, cleanmask=maskname, \
                                           normalization='BLBased',apply_weights=True,use_weightspectrum=False,uvrange=uvrange,
                                           use_dicomodel=True,dicomodel_base=rootimg_old,
                                           RMSFactorInitHMP=1.,MaxMinorIterInitHMP=10000,PredictSettings=("Clean","DD_PREDICT"))

    # make source catalogue
    logger.info('Making source catalogue...')
    maskname = rootimg+'.app.mask.fits'
    make_mask(rootimg+'.app.restored.fits', mask_name=maskname, threshisl=5, atrous_do=False, rmsbox=(30,10), write_srl=True)

    rootimg_old = rootimg

logging.info('Done.')
