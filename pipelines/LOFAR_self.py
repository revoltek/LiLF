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

## here images, models, solutions for each group will be saved
#lib_util.check_rm('self')
#if not os.path.exists('self/images'): os.makedirs('self/images')
#if not os.path.exists('self/solutions'): os.makedirs('self/solutions')
#if not os.path.exists('self/plots'): os.makedirs('self/plots')

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
imgsizepix = int(1.2*MSs.getListObj()[0].getFWHM(freq='min')*3600/10.)
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
logger.info('Creating MODEL_DATA_LOWRES and SUBTRACTED_DATA...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA_LOWRES,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

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
        # DIE Calibration - ms:CORRECTED_DATA
        logger.info('Solving slow G...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.h5parm=$pathMS/g.h5', \
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
        lib_util.run_losoto(s, 'g-c'+str(c), [MS+'/g.h5' for MS in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-leak.parset'])
        os.system('mv plots-g-c'+str(c)+' self/plots/')
        os.system('mv cal-g-c'+str(c)+'.h5 self/solutions/')

        # TEST: correct G - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting G...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=self/solutions/cal-g-c'+str(c)+'.h5 cor.correction=amplitude000', \
                log='$nameMS_corG-c'+str(c)+'.log', commandType='DPPP')

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA

    # clean mask clean (cut at 5k lambda)
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='10arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.8, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=3, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 4)
   
    # baseline averaging possible as we cut longest baselines (also it is in time, where smearing is less problematic)
    # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imgsizepix, scale='10arcsec', \
            weight='briggs 0.', niter=1000000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.8, \
            multiscale='', \
            baseline_averaging=5, auto_threshold=0.5, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')

    print('test apply amplitudeG')
    sys.exit()

    # do beam-corrected+deeper image at last cycle
    if c == 1:

        logger.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        lib_util.run_wsclean(s, 'wscleanBeam-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, temp_dir='./', size=imgsizepix, scale='10arcsec', \
                weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=5000, mgain=0.8, \
                multiscale='', \
                use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=600, \
                parallel_deconvolution=256, auto_threshold=0.5, fits_mask=im.maskname, \
                join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)
        os.system('cat logs/wscleanBeam-c'+str(c)+'.log | grep "background noise"')
        os.system('makepb.py -o img/avgbeam.fits -i '+imagename)
        
        logger.info('Cleaning V (cycle: '+str(c)+')...')
        imagename = 'img/wideV'
        lib_util.run_wsclean(s, 'wscleanV-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='10arcsec', pol='V', \
            weight='briggs 0.', niter=1000, no_update_model_required='', minuv_l=30, maxuv_l=5000, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=600)

    # add model and remove first sidelobe
    if c == 0:

        im = lib_img.Image(imagename+'-MFS-image.fits', beamReg=beamReg)
        im.selectCC(keepInBeam=True)

        # predict only sources in beam - ms: MODEL_DATA
        logger.info('Predict model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.usebeammodel=false pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

        # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    
        # reclean low-resolution
        # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, save_source_list='', temp_dir='./', size=imgsizepix, scale='30arcsec', \
                weight='briggs 0.', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=5000, mgain=0.8, \
                baseline_averaging=5, auto_mask=3, auto_threshold=0.5, \
                join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)
        
        im = lib_img.Image(imagename_lr+'-MFS-image.fits', beamReg=beamReg)
        im.selectCC(keepInBeam=False)

        # predict low-res model for CS-CS baselines only (sols are too inaccurate and data smeared for the RS) - ms: MODEL_DATA_LOWRES
        logger.info('Set MODEL_DATA_LOWRES to 0')
        MSs.run('taql "update $pathMS set MODEL_DATA_LOWRES = 0"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        logger.info('Predict low-res model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msin.baseline="CS*&" msout.datacolumn=MODEL_DATA_LOWRES pre.usebeammodel=false pre.sourcedb='+im.skydb, \
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

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in range(2) ]
[ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in range(2) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wideV-image.fits self/images')
os.system('mv img/wideBeam-MFS-image.fits  img/wideBeam-MFS-image-pb.fits img/avgbeam.fits self/images')
os.system('mv logs self')

logger.info("Done.")
