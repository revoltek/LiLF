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

# here images, models, solutions for each group will be saved
lib_util.check_rm('self')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')
if not os.path.exists('self/plots'): os.makedirs('self/plots')
lib_util.check_rm('img')
os.makedirs('img')

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )
try:
    MSs.plot_HAcov('HAcov.png')
except:
    logger.error('Problem with HAcov, continue anyway.')

# make beam to the first mid null
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg', freq='mid', to_null=True)
beamReg = 'self/beam.reg'

# set image size
imgsizepix = int(2.1*MSs.getListObj()[0].getFWHM(freq='mid')*3600/10.)
if imgsizepix%2 != 0: imgsizepix += 1 # prevent odd img sizes

#################################################################
# Get online model
if sourcedb is None:
    if not os.path.exists('tgts.skydb'):
        fwhm = MSs.getListObj()[0].getFWHM(freq='min')
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
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

# Create columns
logger.info('Creating SUBTRACTED_DATA...')
MSs.run('addcol2ms.py -m $pathMS -c SUBTRACTED_DATA,CORRECTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

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
 
    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving TEC...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS sol.h5parm=$pathMS/tec.h5', \
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DPPP')
 
    # LoSoTo plot dejump
    for MS in MSs.getListObj():
        lib_util.run_losoto(s, 'tec-c'+str(c)+'-'+MS.nameMS, MS.pathMS+'/tec.h5',[parset_dir+'/losoto-tec.parset'])
    os.system('mv plots-tec-c'+str(c)+'* self/plots/')
    s.add('H5parm_collector.py -V -s sol000 -o self/solutions/cal-tec-c'+str(c)+'.h5 '+' '.join(glob.glob('cal-tec-c'+str(c)+'*.h5')),\
            log='losotoTEC-c'+str(c)+'.log', commandType="python", processors='max')
    s.run(check = True)
    lib_util.check_rm('cal-tec'+str(c)+'*.h5')

    # correct TEC - group*_TC.MS:(SUBTRACTED_)DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correcting TEC...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn='+incol+' cor.parmdb=self/solutions/cal-tec-c'+str(c)+'.h5 cor.correction=tec000', \
               log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DPPP')

    # AMP+LEAK DIE correction
    if c >= 0:

        # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Converting to circular...')
        MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)

        # DIE Calibration - ms:CORRECTED_DATA
        logger.info('Solving slow G...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.h5parm=$pathMS/g.h5', \
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
        lib_util.run_losoto(s, 'g-c'+str(c), [MS+'/g.h5' for MS in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-fr.parset', parset_dir+'/losoto-amp.parset'])
        os.system('mv plots-g-c'+str(c)+' self/plots/')
        os.system('mv cal-g-c'+str(c)+'.h5 self/solutions/')

        # Convert back to linear CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Converting to linear...')
        MSs.run('mslin2circ.py -r -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)

        # TEST: correct G - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        #logger.info('Correcting G...')
        #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=self/solutions/cal-g-c'+str(c)+'.h5 cor.correction=amplitudeG', \
        #        log='$nameMS_corG-c'+str(c)+'.log', commandType='DPPP')

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA

    # clean mask clean (cut at 5k lambda)
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='10arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=3, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg, beamReg=beamReg)
    im.makeMask(threshisl = 4, only_beam=True)
   
    # baseline averaging possible as we cut longest baselines (also it is in time, where smearing is less problematic)
    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), do_predict=True, name=imagename, save_source_list='', size=imgsizepix, scale='10arcsec', \
            weight='briggs 0.', niter=1000000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.85, \
            multiscale='', \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=1, auto_mask=3., fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')
       
    # add model and remove first sidelobe
    if c == 0:

        # TEST: reclean low-resolution
        #logger.info('TEST: Cleaning low resolution...')
        #imagename_lr = 'img/TESTpre-wide-lr'
        #lib_util.run_wsclean(s, 'wscleanLR-pre.log', MSs.getStrWsclean(), name=imagename_lr, temp_dir='./', size=imgsizepix, scale='30arcsec', \
        #        weight='briggs 0.', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=5000, mgain=0.8, \
        #        parallel_deconvolution=256, baseline_averaging=5, auto_mask=3, auto_threshold=0.5, \
        #        join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

        # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # Making beam mask
        lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name='img/tmp', size=imgsizepix, scale='30arcsec')
        os.system('mv img/tmp-image.fits img/wide-lr-mask.fits')
        lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 0.)
        lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 1., inverse=True)

        # not all vis are overwritten by wsclean
        logger.info('Reset MODEL_DATA...')
        MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # reclean low-resolution
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), do_predict=True, name=imagename_lr, temp_dir='./', size=imgsizepix, scale='30arcsec', \
                weight='briggs 0.', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=5000, mgain=0.85, \
                parallel_deconvolution=256, baseline_averaging=5, auto_mask=3, auto_threshold=1, fits_mask='img/wide-lr-mask.fits', \
                join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)
        
        ##############################################
        # Flag on empty dataset

        # Subtract low-res model - CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA
        logger.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # Flag on residuals (CORRECTED_DATA)
        logger.info('Flagging residuals...')
        MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS', log='$nameMS_flag-c'+str(c)+'.log', commandType='DPPP')

        ##############################################
        # Prepare SUBTRACTED_DATA

        # corrupt model with TEC solutions - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Corrupt low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                cor.parmdb=self/solutions/cal-tec-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False', \
                log='$nameMS_corrupt.log', commandType='DPPP')
    
        # Subtract low-res model - SUBTRACTED_DATA = DATA - MODEL_DATA
        logger.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # Recreate MODEL_DATA
        logger.info('Predict model...')
        s.add('wsclean -predict -name img/wideM-0 -j '+str(s.max_processors)+' -channels-out 9 '+MSs.getStrWsclean(), \
               log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
        s.run(check=True)

        #MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.usebeammodel=false pre.sourcedb='+im.skydb, \
        #        log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

        # TEST: reclean low-resolution
        #logger.info('TEST: Cleaning low resolution...')
        #imagename_lr = 'img/TESTpost-wide-lr'
        #lib_util.run_wsclean(s, 'wscleanLR-after.log', MSs.getStrWsclean(), name=imagename_lr, temp_dir='./', size=imgsizepix, scale='30arcsec', \
        #        weight='briggs 0.', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=5000, mgain=0.8, \
        #        parallel_deconvolution=256, baseline_averaging=5, auto_mask=3, auto_threshold=0.5, \
        #        join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)


    # do beam-corrected+fullstokes image at last cycle
    if c == 1:

        logger.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        lib_util.run_wsclean(s, 'wscleanBeam-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, temp_dir='./', size=imgsizepix, scale='10arcsec', \
                weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.85, \
                pol='IQUV', join_polarizations='', \
                multiscale='', \
                use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=600, \
                parallel_deconvolution=256, auto_threshold=1, auto_mask=3., fits_mask=im.maskname, \
                join_channels='', channels_out=9)
        os.system('cat logs/wscleanBeam-c'+str(c)+'.log | grep "background noise"')
        os.system('makepb.py -o img/avgbeam.fits -i '+imagename)
 

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in range(2) ]
[ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in range(2) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wideBeam-MFS-*-image.fits  img/wideBeam-MFS-*-image-pb.fits img/avgbeam.fits self/images')
os.system('mv logs self')

logger.info("Done.")
