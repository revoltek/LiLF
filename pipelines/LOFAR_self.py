#!/usr/bin/env python
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal corrected source subtracted data in CORRECTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions
# last high/low resolution models are copied in the "self/models" dir
# last high/low resolution images + masks + empty images (CORRECTED_DATA) are copied in the "self/images" dir
# h5parm solutions and plots are copied in the "self/solutions" dir

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt

# Temporary
#if 'bootes' in os.getcwd():
#    sourcedb = '/home/fdg/scripts/model/Bootes_HBA.corr.skydb'
#    apparent = False
#    userReg = None
#elif '3Csurvey' in os.getcwd():
#    sourcedb = '/home/fdg/scripts/model/3c31_LBA.skydb'
#    apparent = False
#    userReg = None
if 'LBAsurvey' in os.getcwd():
    obs = os.getcwd().split('/')[-1]
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/skymodels/%s.skydb' % obs
    apparent = False
    userReg = None
    if not os.path.exists('mss'):
        os.makedirs('mss')
        for i, tc in enumerate(glob.glob('../../c*-o*/%s/mss/*' % obs)):
            tc_ren = 'TC%02i.MS' % i
            print 'cp -r %s mss/%s' % (tc,tc_ren)
            os.system('cp -r %s mss/%s' % (tc,tc_ren))

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
lib_log.set_logger('pipeline-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(dry = False)

parset = lib_util.getParset()
parset_dir = parset.get('self','parset_dir')
sourcedb = parset.get('model','sourcedb')
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')
assert os.path.exists(sourcedb)

niter = 3

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

# se timage size
obsmode = MSs.getListObj()[0].getObsMode()
if 'INNER' in obsmode: size = 6000
elif 'SPARSE' in obsmode: size = 4500
elif 'OUTER' in obsmode: size = 3000
else: 
    logging.error('Observing mode not recognised.')
    sys.exit(1)

##################################################################################################
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
for c in xrange(0, niter):

#    logger.info('Start selfcal cycle: '+str(c))
#
#    if c >= 2:
#        incol = 'SUBTRACTED_DATA'
#    else:
#        incol = 'DATA'
#
#    # Smooth DATA -> SMOOTHED_DATA
#    logger.info('BL-based smoothing...')
#    MSs.run('BLsmooth.py -r -f 0.2 -i '+incol+' -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1-c'+str(c)+'.log', commandType='python', maxThreads=6)

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving TEC...')
    for MS in MSs.getListStr():
        lib_util.check_rm(MS+'/tec.h5')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTECdd.parset msin=$pathMS ddecal.h5parm=$pathMS/tec.h5', \
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DPPP')

    # LoSoTo plot
    lib_util.run_losoto(s, 'tec'+str(c), [MS+'/tec.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-tec.parset'])
    os.system('mv plots-tec'+str(c)+'* self/solutions/')
    os.system('mv cal-tec'+str(c)+'*.h5 self/solutions/')

    # correct TEC - group*_TC.MS:(SUBTRACTED_)DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correcting TEC...')
    MSs.run('DPPP '+parset_dir+'/DPPP-corTEC.parset msin=$pathMS msin.datacolumn='+incol+' cor1.parmdb=$pathMS/tec.h5 cor2.parmdb=$pathMS/tec.h5', \
                log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DPPP')

    #####################################################################################################
    # Faraday rotation correction
    if c >= 1:

#        # To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (circular)
#        logger.info('Convert to circular...')
#        MSs.run('/home/fdg/scripts/mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin-c'+str(c)+'.log', commandType='python', maxThreads=4)
 
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2-c'+str(c)+'.log', commandType='python', maxThreads=6)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Solving G...')
        for MS in MSs.getListStr():
            lib_util.check_rm(MS+'/fr.h5')
        MSs.run('DPPP ' + parset_dir + '/DPPP-solGdd.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.mode=rotation+diagonal \
                     sol.solint=30 sol.nchan=8', log='$nameMS_solFR.log', commandType="DPPP")

        lib_util.run_losoto(s, 'fr'+str(c), [MS+'/fr.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-fr.parset'])
        os.system('mv plots-fr'+str(c)+'* self/solutions/')
        os.system('mv cal-fr'+str(c)+'*.h5 self/solutions/')
       
#        # To linear - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (linear)
#        logger.info('Convert to linear...')
#        MSs.run('/home/fdg/scripts/mslin2circ.py -r -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', \
#                log='$nameMS_circ2lin-c'+str(c)+'.log', commandType='python', maxThreads=4)
        
        # Correct FR SB.MS:(SUBTRACTED_)DATA -> CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        h5 = 'self/solutions/cal-fr'+str(c)+'.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn='+incol+' cor.parmdb='+h5+' cor.correction=rotationmeasure000', \
                    log='$nameMS_corFR-c'+str(c)+'.log', commandType='DPPP')

        #################################
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3-c'+str(c)+'.log', commandType='python', maxThreads=6)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Solving G...')
        for MS in MSs.getListStr():
            lib_util.check_rm(MS+'/amp.h5')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.parmdb=$pathMS/amp.h5 sol.solint=30 sol.nchan=8', \
                    log='$nameMS_sol-g2-c'+str(c)+'.log', commandType='DPPP')

        lib_util.run_losoto(s, 'amp'+str(c), [MS+'/amp.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-amp.parset'])
        os.system('mv plots-amp'+str(c)+'* self/solutions/')
        os.system('mv cal-amp'+(str(c))+'*.h5 self/solutions/')

        # Correct beam amp SB.MS:SUBTRACTED_DATA->CORRECTED_DATA
        #logger.info('Beam amp correction...')
        #h5 = 'self/solutions/cal-amp'+str(c)+'.h5'
        #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA cor.parmdb='+h5+' cor.correction=amplitude000', \
        #        log='$nameMS_corAMP-c'+str(c)+'.log', commandType='DPPP')
        ## Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        #logger.info('Faraday rotation correction...')
        #h5 = 'self/solutions/cal-fr'+str(c)+'.h5'
        #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb='+h5+' cor.correction=rotationmeasure000', \
        #            log='$nameMS_corFR2-c'+str(c)+'.log', commandType='DPPP')
        #
        ## Finally re-calculate TEC
        #logger.info('BL-based smoothing...')
        #MSs.run('BLsmooth.py -r -f 0.2 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3-c'+str(c)+'.log', commandType='python', maxThreads=6)

        # solve TEC - group*_TC.MS:SMOOTHED_DATA
        logger.info('Solving TEC...')
        for MS in MSs.getListStr():
            lib_util.check_rm(MS+'/tec.h5')
        MSs.run('DPPP '+parset_dir+'/DPPP-solTECdd.parset msin=$pathMS ddecal.h5parm=$pathMS/tec.h5', \
                    log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DPPP')

        # LoSoTo plot
        lib_util.run_losoto(s, 'tec'+str(c)+'b', [MS+'/tec.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-tec.parset'])
        os.system('mv plots-tec'+str(c)+'b* self/solutions')
        os.system('mv cal-tec'+str(c)+'b*.h5 self/solutions')

        # correct TEC - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting TEC...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corTEC.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor1.parmdb=$pathMS/tec.h5 cor2.parmdb=$pathMS/tec.h5', \
                    log='$nameMS_corTECb-c'+str(c)+'.log', commandType='DPPP')

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

    # do beam-corrected+deeper image at last cycle
    if c == niter-1:
        # beam corrected: -use-differential-lofar-beam' - no baseline avg!
        logger.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(int(size*1.5)) + ' ' + str(int(size*1.5)) + ' -mem 90 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
                -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -apply-primary-beam -use-differential-lofar-beam -minuv-l 30 '+MSs.getStrWsclean(), \
                log='wscleanBeam-c'+str(c)+'.log', commandType='wsclean', processors='max')
        s.run(check=True)
        os.system('cat logs/wscleanBeam-c'+str(c)+'.log | grep "background noise"')

        logger.info('Cleaning beam high-res (cycle: '+str(c)+')...')
        imagename = 'img/wideBeamHR'
        s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size*2) + ' ' + str(size*2) + ' -mem 90 -j '+str(s.max_processors)+' \
                -scale 2.5arcsec -weight briggs -1.5 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -apply-primary-beam -use-differential-lofar-beam -minuv-l 30 '+MSs.getStrWsclean(), \
                log='wscleanBeamHR-c'+str(c)+'.log', commandType='wsclean', processors='max')
        s.run(check=True)

        logger.info('Cleaning beam low-res (cycle: '+str(c)+')...')
        imagename = 'img/wideBeamLR'
        s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size/5) + ' ' + str(size/5) + ' -mem 90 -j '+str(s.max_processors)+' \
                -scale 60arcsec -weight briggs 1.5 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -pol IUQV -join-channels -fit-spectral-pol 2 -channels-out 10 -apply-primary-beam -use-differential-lofar-beam -maxuv-l 1000 '+MSs.getStrWsclean(), \
                log='wscleanBeamLR-c'+str(c)+'.log', commandType='wsclean', processors='max')
        s.run(check=True)


    # clean mask clean (cut at 5k lambda)
    # no MODEL_DATA update with -baseline-averaging
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size) + ' ' + str(size) + ' -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 10arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.9 \
            -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 20 -minuv-l 30 '+MSs.getStrWsclean(), \
            log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 3)

    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size) + ' ' + str(size) + ' -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 10arcsec -weight briggs 0.0 -niter 1000000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
            -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 0.1 -minuv-l 30 -save-source-list -fits-mask '+im.maskname+' '+MSs.getStrWsclean(), \
            log='wscleanM-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-c'+str(c)+'.log | grep "background noise"')

    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg, beamReg=beamReg)
    im.selectCC(keepInBeam=True)

    # predict
    logger.info('Predict (ft)...')
    if c != niter:
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.usebeammodel=false pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    if c == 1:
        # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')
    
        # reclean low-resolution
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        s.add('wsclean -reorder -name ' + imagename_lr + ' -size ' + str(size) + ' ' + str(size) + ' -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
                -scale 20arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 2000 -mgain 0.8 \
                -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 1 -minuv-l 100 -save-source-list '+MSs.getStrWsclean(), \
                log='wsclean-lr.log', commandType='wsclean', processors='max')
        s.run(check=True)
        
        im = lib_img.Image(imagename_lr+'-MFS-image.fits', beamReg=beamReg)
        im.selectCC(keepInBeam=False)

        # predict - ms: MODEL_DATA_LOWRES
        logger.info('Predict low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_LOWRES pre.usebeammodel=false pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-lr.log', commandType='DPPP')

        # corrupt model with TEC solutions - ms:MODEL_DATA_LOWRES -> ms:MODEL_DATA_LOWRES
        logger.info('Corrupt low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corTEC.parset msin=$pathMS msin.datacolumn=MODEL_DATA_LOWRES msout.datacolumn=MODEL_DATA_LOWRES  \
                cor1.parmdb=$pathMS/tec.h5 cor1.invert=false cor2.parmdb=$pathMS/tec.h5 cor2.invert=false', \
                log='$nameMS_corrupt.log', commandType='DPPP')
    
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA_LOWRES -> concat.MS:CORRECTED_DATA (empty)
        logger.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES)...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA_LOWRES"', log='$nameMS_taql2-c'+str(c)+'.log', commandType='general')

    ###############################################################################################################
    # TODO
    # Flag on residuals (CORRECTED_DATA)
    #logger.info('Flagging residuals...')
    #MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS', log='$nameMS_flag-c'+str(c)+'.log', commandType='DPPP')
    
# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in xrange(niter) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wideBeam-MFS-image.fits  img/wideBeam-MFS-image-pb.fits self/images')
os.system('mv img/wideBeamHR-MFS-image.fits  img/wideBeamHR-MFS-image-pb.fits self/images')
os.system('mv logs self')

logger.info("Done.")
