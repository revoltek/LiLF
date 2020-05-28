#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation.
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import sys, os, glob, re
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-cal.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-cal.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_cal','parset_dir')
data_dir = parset.get('LOFAR_cal','data_dir')
skymodel = parset.get('LOFAR_cal','skymodel')
imaging = parset.getboolean('LOFAR_cal','imaging')
bl2flag = parset.get('flag','stations')

#############################################################
MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s, check_flags = False )
# copy data
logger.info('Copy data...')
for MS in MSs.getListObj():
    if min(MS.getFreqs()) > 30.e6:
        MS.move(MS.nameMS+'.MS', keepOrig=True, overwrite=False)

MSs = lib_ms.AllMSs( glob.glob('*MS'), s, check_flags = False )
calname = MSs.getListObj()[0].getNameField()
for MS in MSs.getListObj():
    os.system('cp -r %s %s' % (skymodel, MS.pathMS))

######################################################
# flag bad stations, flags will propagate
if w.todo('flag'):
    logger.info("Flagging...")
    MSs.run("DPPP " + parset_dir + "/DPPP-flag.parset msin=$pathMS ant.baseline=\"" + bl2flag+"\"", log="$nameMS_flag.log", commandType="DPPP")
    
    # extend flags
    logger.info('Remove bad time/freq stamps...')
    MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    w.done('flag')
### DONE

if w.todo('predict'):
    # predict to save time ms:MODEL_DATA
    logger.info('Add model to MODEL_DATA (%s)...' % calname)
    MSs.run("DPPP " + parset_dir + "/DPPP-predict.parset msin=$pathMS pre.sourcedb=$pathMS/" + os.path.basename(skymodel) + " pre.sources=" + calname, \
            log="$nameMS_pre.log", commandType="DPPP", maxThreads=30)
    
    w.done('predict')
### DONE

###################################################
# 1: find PA

if w.todo('cal_pa'):
    # Smooth data DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -c 1 -n 8 -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType
    ='python', maxThreads=8)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating PA...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal', log='$nameMS_solPA.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'pa', [ms+'/pa.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])
    
    # Pol align correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DPPP")

    w.done('cal_pa')
### DONE

########################################################
# 2: find FR

if w.todo('cal_fr'):
    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType="DPPP")
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -c 1 -n 8 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2.log',
            commandType ='python', maxThreads=8)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating FR...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.mode=rotation+diagonal', log='$nameMS_solFR.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'fr', [ms+'/fr.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-fr.parset'])
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DPPP")

    w.done('cal_fr')
### DONE

#######################################################
# 4: find BP

if w.todo('cal_bp'):
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -c 1 -n 8 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3.log',
            commandType ='python', maxThreads=8)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating BP...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/amp.h5 sol.mode=diagonal', log='$nameMS_solAMP.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'amp', [ms+'/amp.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-flag.parset',parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-plot-ph.parset',parset_dir+'/losoto-bp.parset'])

    w.done('cal_bp')
### DONE

####################################################
# Re-do corrections in right order

if w.todo('apply_all'):
    # Pol align correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA2.log', commandType="DPPP")
    
    # Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('AmpBP correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-amp.h5 \
            cor.correction=amplitudeSmooth cor.updateweights=True', log='$nameMS_corAMP.log', commandType="DPPP")
    
    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam2.log', commandType="DPPP")
    
    # # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    # logger.info('Faraday rotation correction...')
    # MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DPPP")
    
    w.done('apply_all')
### DONE

#################################################
# 4: find iono

# if w.todo('cal_iono'):
#     # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
#     logger.info('BL-smooth...')
#     MSs.run('BLsmooth.py -r -c 1 -n 8 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth4.log',
#             commandType ='python', maxThreads=8)
#
#     # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#     logger.info('Calibrating IONO...')
#     MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/iono.h5 sol.mode=diagonal', log='$nameMS_solIONO.log', commandType="DPPP")
#
#     lib_util.run_losoto(s, 'iono', [ms+'/iono.h5' for ms in MSs.getListStr()], \
#             [parset_dir+'/losoto-flag.parset', parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-iono.parset'])
#
#     w.done('cal_iono')
# ### DONE

# a debug image
if imaging:
    logger.info("Imaging section:")

    # Correct all CORRECTED_DATA (PA, beam, FR, BP corrected) -> CORRECTED_DATA
    logger.info('IONO correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.steps=[ph,amp] cor.ph.parmdb=cal-iono.h5 cor.amp.parmdb=cal-iono.h5 \
        cor.ph.correction=phaseOrig000 cor.amp.correction=amplitude000 cor.amp.updateweights=False', log='$nameMS_corIONO.log', commandType="DPPP")

    logger.info('Subtract model...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql2.log', commandType ='general')

    logger.info('Flag...')
    MSs.run('DPPP '+parset_dir+'/DPPP-flag2.parset msin=$pathMS', log='$nameMS_flag2.log', commandType='DPPP')

    lib_util.check_rm('img')
    os.makedirs('img')

    imgsizepix =  MSs.getListObj()[0].getFWHM()*3600/5

    logger.info('Cleaning low-res...')
    imagename = 'img/calLR'
    lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix/5, scale='60arcsec', \
            weight='briggs 0.', taper_gaussian='240arcsec', niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=2000, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            parallel_deconvolution=512, \
            auto_mask=10, auto_threshold=1, pol='IQUV', join_polarizations='', join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    logger.info('Cleaning normal...')
    imagename = 'img/cal'
    lib_util.run_wsclean(s, 'wscleanA.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='5arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', minuv_l=30, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=512, \
            auto_threshold=20, join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits')
    im.makeMask(threshisl = 3)

    logger.info('Cleaning w/ mask...')
    lib_util.run_wsclean(s, 'wscleanB.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='5arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=512, \
            auto_threshold=0.1, fits_mask=im.maskname, join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)
    os.system('cat logs/wscleanA.log logs/wscleanB.log | grep "background noise"')

    # make new mask
    im.makeMask(threshisl = 5)

logger.info("Done.")
