#!/usr/bin/env python
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

# Temporary
if 'LBAsurvey' in os.getcwd():
    obs = os.getcwd().split('/')[-1]
    if not os.path.exists('mss'):
        os.makedirs('mss')
        for i, tc in enumerate(glob.glob('../../c*-o*/%s/mss/*' % obs)):
            tc_ren = 'TC%02i.MS' % i
            print('cp -r %s mss/%s' % (tc,tc_ren))
            os.system('cp -r %s mss/%s' % (tc,tc_ren))

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_self','parset_dir')
sourcedb = parset.get('model','sourcedb')
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')

#############################################################################
# Clear
logger.info('Cleaning...')

lib_util.check_rm('img')
os.makedirs('img')

# here images, models, solutions for each group will be saved
lib_util.check_rm('self')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg') # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null
beamReg = 'self/beam.reg'

# set image size
imgsizepix =  1.2*MSs.getListObj()[0].getFWHM()*3600/10

#################################################################
# Get online model
if sourcedb is None:
    if not os.path.exists('tgts.skydb'):
        fwhm = MSs.getListObj()[0].getFWHM()
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "http://172.104.228.177/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f"' % (radeg, decdeg, fwhm/2.))
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<1')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')
        apparent = False

    sourcedb = 'tgts.skydb'

#################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/'+sourcedb_basename)
    logger.debug('Copy: '+sourcedb+' -> '+MS)
    os.system('cp -r '+sourcedb+' '+MS)

# Create columns (non compressed)
logger.info('Creating MODEL_DATA_LOWRES and SUBTRACTED_DATA...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA_LOWRES,SUBTRACTED_DATA', log='$nameMS_addcol.log', commandType='python')

logger.info('Add model to MODEL_DATA...')
if apparent:
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.usebeammodel=false pre.sourcedb=$pathMS/'+sourcedb_basename, log='$nameMS_pre.log', commandType='DPPP')
else:
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.usebeammodel=true pre.sourcedb=$pathMS/'+sourcedb_basename, log='$nameMS_pre.log', commandType='DPPP')

#####################################################################################################
# Self-cal cycle
for c in range(2):

    logger.info('Start selfcal cycle: '+str(c))

    if c != 0:
        incol = 'SUBTRACTED_DATA'
    else:
        incol = 'DATA'

    # Smooth DATA -> SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -r -f 0.2 -i '+incol+' -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1-c'+str(c)+'.log', commandType='python')

    #logger.info('Solving G...')
    #MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS ddecal.h5parm=$pathMS/tecmax70km.h5 ddecal.solint=3 ddecal.nchan=1', \
    #            log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
    #MSs.run('DPPP '+parset_dir+'/DPPP-solGdd.parset msin=$pathMS sol.h5parm=$pathMS/gspmax70kmavg3x12.h5 sol.solint=3 sol.nchan=12 sol.mode=scalarphase sol.uvmmax=70e3', \
    #            log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
    #sys.exit()

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving TEC...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS ddecal.h5parm=$pathMS/tec.h5', \
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DPPP')

    # LoSoTo plot
    lib_util.run_losoto(s, 'tec'+str(c), [MS+'/tec.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-tec.parset'])
    os.system('mv plots-tec'+str(c)+'* self/solutions/')
    os.system('mv cal-tec'+str(c)+'*.h5 self/solutions/')

    # correct TEC - group*_TC.MS:(SUBTRACTED_)DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correcting TEC...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn='+incol+' cor.parmdb=self/solutions/cal-tec'+str(c)+'.h5 cor.correction=tec000', \
                log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DPPP')

    #####################################################################################################
    # Faraday rotation correction
    #if c >= 0:
     
    #    # To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (circular)
    #    logger.info('Convert to circular...')
    #    MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', \
    #            log='$nameMS_circ2lin-c'+str(c)+'.log', commandType='python', maxThreads=4)
 
    #    # Smooth CORRECTED_DATA -> SMOOTHED_DATA
    #    logger.info('BL-based smoothing...')
    #    MSs.run('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2-c'+str(c)+'.log', commandType='python')

    #    # Solve G SB.MS:SMOOTHED_DATA (only solve)
    #    logger.info('Solving G...')
    #    #MSs.run('DPPP ' + parset_dir + '/DPPP-solGdd.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.mode=rotation+diagonal sol.solint=30 sol.nchan=8 \
    #    MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.parmdb=$pathMS/fr.h5 sol.solint=30 sol.nchan=8', \
    #                 log='$nameMS_solFR.log', commandType="DPPP")

    #    #lib_util.run_losoto(s, 'fr'+str(c), [MS+'/fr.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-fr.parset'])
    #    lib_util.run_losoto(s, 'fr'+str(c), [MS+'/fr.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-fr.parset'])
    #    os.system('mv plots-fr'+str(c)+'* self/solutions/')
    #    os.system('mv cal-fr'+str(c)+'*.h5 self/solutions/')
    #   
    #    # Correct FR SB.MS:(SUBTRACTED_)DATA -> CORRECTED_DATA
    #    logger.info('Faraday rotation correction...')
    #    h5 = 'self/solutions/cal-fr'+str(c)+'.h5'
    #    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn='+incol+' cor.parmdb='+h5+' cor.correction=rotationmeasure000', \
    #                log='$nameMS_corFR-c'+str(c)+'.log', commandType='DPPP')

    #    #################################
    #    # Smooth CORRECTED_DATA -> SMOOTHED_DATA
    #    logger.info('BL-based smoothing...')
    #    MSs.run('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3-c'+str(c)+'.log', commandType='python')

    #    # Solve G SB.MS:SMOOTHED_DATA (only solve)
    #    logger.info('Solving G...')
    #    MSs.run('DPPP '+parset_dir+'/DPPP-solGdd.parset msin=$pathMS sol.h5parm=$pathMS/amp.h5 sol.solint=60 sol.nchan=8', \
    #                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')

    #    lib_util.run_losoto(s, 'amp'+str(c), [MS+'/amp.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-amp.parset'])
    #    os.system('mv plots-amp'+str(c)+'* self/solutions/')
    #    os.system('mv cal-amp'+(str(c))+'*.h5 self/solutions/')

    #    # Correct beam amp SB.MS:SUBTRACTED_DATA->CORRECTED_DATA
    #    #logger.info('Beam amp correction...')
    #    #h5 = 'self/solutions/cal-amp'+str(c)+'.h5'
    #    #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA cor.parmdb='+h5+' cor.correction=amplitude000', \
    #    #        log='$nameMS_corAMP-c'+str(c)+'.log', commandType='DPPP')

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

    # do beam-corrected+deeper image at last cycle
    # TODO: find a way to save the beam image
    if c == niter-1:
        logger.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        lib_util.run_wsclean(s, 'wscleanBeam-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, temp_dir='./', size=imgsizepix, scale='10arcsec', \
                weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.85, \
                multiscale='', multiscale_scale_bias=0.5, multiscale_scales='0,10,20', \
                use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
                parallel_deconvolution=256, auto_mask=10, auto_threshold=1, \
                join_channels='', fit_spectral_pol=2, channels_out=8)
        os.system('cat logs/wscleanBeam-c'+str(c)+'.log | grep "background noise"')

        #logger.info('Cleaning beam high-res (cycle: '+str(c)+')...')
        #imagename = 'img/wideBeamHR'
        #lib_util.run_wsclean(s, 'wscleanBeamHR-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, temp_dir='./', size=int(imgsizepix*2), scale='2.5arcsec', \
        #        weight='uniform', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
        #        use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
        #        parallel_deconvolution=256, auto_mask=10, auto_threshold=1, \
        #        join_channels='', fit_spectral_pol=2, channels_out=8)
        #
        #logger.info('Cleaning beam low-res (cycle: '+str(c)+')...')
        #imagename = 'img/wideBeamLR'
        #lib_util.run_wsclean(s, 'wscleanBeamLR-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, temp_dir='./', size=imgsizepix/5, scale='60arcsec', \
        #        weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=1000, mgain=0.85, \
        #        use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
        #        parallel_deconvolution=256, auto_mask=10, auto_threshold=1, pol='IQUV', \
        #        join_channels='', fit_spectral_pol=2, channels_out=8)

    # clean mask clean (cut at 5k lambda)
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='10arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=20, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg, beamReg=beamReg)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 4)
    
    # baseline averaging possible as we cut longest baselines (also it is in time, where smearing is less problematic)
    # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imgsizepix, scale='10arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.85, \
            multiscale='', multiscale_scales='0,10,20', \
            baseline_averaging=5, auto_threshold=1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')

    if c != niter-1:

        im = lib_img.Image(imagename+'-MFS-image.fits', beamReg=beamReg)
        im.selectCC(keepInBeam=True)

        # predict only sources in beam - ms: MODEL_DATA
        logger.info('Predict model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.usebeammodel=false pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    if c == 0:
        # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    
        # reclean low-resolution
        # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, save_source_list='', temp_dir='./', size=imgsizepix, scale='20arcsec', \
                weight='briggs 0.', niter=50000, no_update_model_required='', minuv_l=30, maxuv_m=5000, mgain=0.85, \
                baseline_averaging=5, auto_threshold=1, \
                join_channels='', fit_spectral_pol=2, channels_out=8)
        
        im = lib_img.Image(imagename_lr+'-MFS-image.fits', beamReg=beamReg)
        im.selectCC(keepInBeam=False)

        # predict low-res model for CS-CS baselines only (sols are too inaccurate and data smeared for the RS) - ms: MODEL_DATA_LOWRES
        logger.info('Set MODEL_DATA_LOWRES to 0')
        MSs.run('taql "update $pathMS set MODEL_DATA_LOWRES = 0"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        logger.info('Predict low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msin.baselines=CS*& msout.datacolumn=MODEL_DATA_LOWRES pre.usebeammodel=false pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-lr.log', commandType='DPPP')

        ##############################################
        # Flag on empty dataset

        # Subtract low-res model - CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_LOWRES
        logger.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_LOWRES)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_LOWRES"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # Flag on residuals (CORRECTED_DATA)
        logger.info('Flagging residuals...')
        MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS', log='$nameMS_flag-c'+str(c)+'.log', commandType='DPPP')

        ##############################################
        # Prepare SUBTRACTED_DATA

        # corrupt model with TEC solutions - ms:MODEL_DATA_LOWRES -> ms:MODEL_DATA_LOWRES
        logger.info('Corrupt low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA_LOWRES msout.datacolumn=MODEL_DATA_LOWRES  \
                cor.parmdb=$pathMS/tec.h5 cor.correction=tec000 cor.invert=false', \
                log='$nameMS_corrupt.log', commandType='DPPP')
    
        # Subtract low-res model - SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES
        logger.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES)...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # TODO: remove only from short baselines

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in range(niter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in range(niter) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wideBeam-MFS-image.fits  img/wideBeam-MFS-image-pb.fits self/images')
os.system('mv logs self')

logger.info("Done.")
