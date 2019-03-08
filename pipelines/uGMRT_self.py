#!/usr/bin/env python
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

parset = lib_util.getParset()
parset_dir = parset.get('uGMRT_self','parset_dir')
sourcedb = parset.get('model','sourcedb') # relative to tgts dir "./tgts/xxx/model.txt
userReg = parset.get('model','userReg') # relative to tgts dir "./tgts/xxx/region.ref"

niter = 10

tgts = glob.glob('tgts/*')

#############################################################################
# Clear
logger.info('Cleaning...')

lib_util.check_rm('img')
os.makedirs('img')

# here images, models, solutions for each group will be saved
lib_util.check_rm('self')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

MSs = lib_ms.AllMSs( glob.glob('mss/*.MS'), s )

# make beam
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg')
beamReg = 'self/beam.reg'

# set image size
imgsizepix =  1.2*MSs.getListObj()[0].getFWHM()*3600/2.

#################################################################
# Get online model
if sourcedb is None:
    if not os.path.exists('tgts.skydb'):
        fwhm = 3 # deg
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "http://172.104.228.177/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f"' % (radeg, decdeg, fwhm/2.))
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<0.1')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')
        apparent = False

    sourcedb = 'tgts.skydb'

##################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/'+sourcedb_basename)
    logger.debug('Copy: '+sourcedb+' -> '+MS)
    os.system('cp -r '+sourcedb+' '+MS)

logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA,MODEL_DATA_LOWRES,SUBTRACTED_DATA -i DATA', log="$nameMS_addcol.log", commandType="python")

logger.info('Add model to MODEL_DATA...')
MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=$pathMS/'+sourcedb_basename, log='$nameMS_pre.log', commandType='DPPP')

#####################################################################################################
# Self-cal cycle
for c in range(0, niter):

    if c > 1:
        incol = 'SUBTRACTED_DATA'
    else:
        incol = 'DATA'

    logger.info('Start selfcal cycle: '+str(c))

    # Smooth DATA -> SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -r -f 0.2 -i '+incol+' -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1-c'+str(c)+'.log', commandType='python')

    # solve - concat*.MS:SMOOTHED_DATA
    logger.info('Solving G...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solGdd.parset msin=$pathMS sol.h5parm=$pathMS/gs.h5 sol.solint=1 sol.nchan=1 sol.mode=complexgain sol.smoothnessconstraint=1e6', \
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
    MSs.run('DPPP '+parset_dir+'/DPPP-solGdd.parset msin=$pathMS sol.h5parm=$pathMS/tecs.h5 sol.solint=1 sol.nchan=1 sol.mode=tec sol.smoothnessconstraint=1e6', \
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')

    # LoSoTo plot
    lib_util.run_losoto(s, 'gs'+str(c), [MS+'/gs.h5' for MS in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-flag.parset'])
    os.system('mv plots-gs'+str(c)+'* self/solutions/')
    os.system('mv cal-gs'+str(c)+'*.h5 self/solutions/')
    lib_util.run_losoto(s, 'tecs'+str(c), [MS+'/tecs.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-tec.parset'])
    os.system('mv plots-tecs'+str(c)+'* self/solutions/')
    os.system('mv cal-tecs'+str(c)+'*.h5 self/solutions/')

    # correct phases - MS:DATA -> MS:CORRECTED_DATA
    logger.info('Correcting Gp...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.correction=phase000', \
                log='$nameMS_corGp-c'+str(c)+'.log', commandType='DPPP')
    if c>2:
        # correct amplitudes - MS:DATA -> MS:CORRECTED_DATA
        logger.info('Correcting Ga...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.correction=amplitude000', \
                log='$nameMS_corGa-c'+str(c)+'.log', commandType='DPPP')

    # clean mask clean
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='2arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=20, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl=4, atrous_do=False)
    
    # baseline averaging possible as we cut longest baselines (also it is in time, where smearing is less problematic)
    # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imgsizepix, scale='2arcsec', \
            weight='briggs 0.', niter=300000, no_update_model_required='', mgain=0.85, \
            #multiscale='', multiscale_scales='0,5,10,20,40', \
            baseline_averaging=5, auto_threshold=1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')

    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl=5, atrous_do=False)
    im.selectCC()

    if c != niter-1:
        # predict - ms: MODEL_DATA
        logger.info('Predict model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    if c == 1:
        # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')
    
        # reclean low-resolution
        # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, save_source_list='', temp_dir='./', size=imgsizepix, scale='10arcsec', \
                weight='briggs 0.', niter=100000, no_update_model_required='', maxuv_l=5000, mgain=0.85, \
                baseline_averaging=5, auto_threshold=0.5, \
                join_channels='', fit_spectral_pol=2, channels_out=8)
        
        im = lib_img.Image(imagename_lr+'-MFS-image.fits', beamReg=beamReg)
        im.makeMask(threshisl=5, atrous_do=False)
        im.selectCC(keepInBeam=False)

        # predict - ms: MODEL_DATA_LOWRES
        # must be done with DPPP to remove sources in beam
        logger.info('Predict low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_LOWRES pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-lr.log', commandType='DPPP')

        ##############################################
        # Flag on empty dataset

        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA_LOWRES -> concat.MS:CORRECTED_DATA (empty)
        #logger.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_LOWRES)...')
        #MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_LOWRES"', log='$nameMS_taql2-c'+str(c)+'.log', commandType='general')

        # Flag on residuals (CORRECTED_DATA)
        #logger.info('Flagging residuals...')
        #MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS', log='$nameMS_flag-c'+str(c)+'.log', commandType='DPPP')

        ##############################################
        # Prepare SUBTRACTED_DATA

        # corrupt model with TEC solutions - ms:MODEL_DATA_LOWRES -> ms:MODEL_DATA_LOWRES
        logger.info('Corrupt low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corG.parset msin=$pathMS msin.datacolumn=MODEL_DATA_LOWRES msout.datacolumn=MODEL_DATA_LOWRES  \
                cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.invert=false', \
                log='$nameMS_corrupt.log', commandType='DPPP')
    
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA_LOWRES -> concat.MS:CORRECTED_DATA (empty)
        logger.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES)...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES"', log='$nameMS_taql3-c'+str(c)+'.log', commandType='general')


# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in xrange(niter) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv logs self')

logger.info("Done.")
