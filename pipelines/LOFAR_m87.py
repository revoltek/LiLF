#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob, re
import numpy as np

patch = 'VirA'
nouseblrange = ''
f = lambda nu: 10**(2.4466 - 0.8116 * ((np.log10(nu/1.e9))**1) - 0.0483 * ((np.log10(nu/1.e9))**2) ) # PB17

skymodel = '/home/fdg/scripts/model/A-team_4_CC.skydb'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-m87.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_m87'])))
parset_dir = parset.get('LOFAR_m87','parset_dir')
data_dir = parset.get('LOFAR_m87','data_dir')
bl2flag = parset.get('flag','stations')

##########################################################
logger.info('Cleaning...')
lib_util.check_rm('cal*h5')
lib_util.check_rm('plots*')
lib_util.check_rm('img')
lib_util.check_rm('mss')
os.makedirs('img')
os.makedirs('mss')

MSs = lib_ms.AllMSs( sorted(glob.glob(data_dir+'/*MS')), s )
# copy data (avg to 1ch/sb and 10 sec)
nchan = 1 #int(MSs.getListObj()[0].getNchan()) # no avg in freq
timeint = MSs.getListObj()[0].getTimeInt()
avg_time = int(np.rint(24./timeint)) # average to 24 s

logger.info('Copy data...')
mss_bkp = sorted(glob.glob(data_dir+'/*MS'))
MSs.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msout=mss/$nameMS avg.freqstep=%i avg.timestep=%i' % (nchan, avg_time),\
            log='$nameMS_avg.log', commandType='DP3')

################################################################
MSs = lib_ms.AllMSs( glob.glob('mss/*MS'), s )

# HBA/LBA
if MSs.isLBA():
    flag_steps = "[ant, uvmin, elev, count]"
else: 
    flag_steps = "[ears, ant, uvmin, elev, count]"

########################################################   
# flag bad stations, and low-elev
logger.info('Flagging...')
MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. steps=\"'+flag_steps+'\" ant.baseline=\"'+bl2flag+'\"', \
            log='$nameMS_flag.log', commandType='DP3')

if MSs.isHBA(): model_dir = '/home/fdg/scripts/model/AteamHBA/'+patch
else: model_dir = '/home/fdg/scripts/model/AteamLBA/'+patch

if os.path.exists(model_dir+'/img-MFS-model.fits'):
    im = lib_img.Image(model_dir+'/img-MFS-image.fits')
    im.rescaleModel(f)
    n = len(glob.glob(model_dir+'/img-[0-9]*-model.fits'))
    logger.info('Predict (wsclean: %s - chan: %i)...' % (model_dir, n))
    s.add('wsclean -predict -name '+model_dir+'/img -j '+str(s.max_processors)+' -channels-out '+str(n)+' '+MSs.getStrWsclean(), \
          log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    s.run(check=True)
else:
    logger.info('Predict (DP3)...')
    MSs.run('DP3 '+parset_dir+'/DP3-predict.parset msin=$pathMS pre.sourcedb='+skymodel+' pre.sources='+patch, log='$nameMS_pre.log', commandType='DP3')

# TESTTESTTEST
# BL Smooth DATA -> DATA
#if lofar_system == 'hba':
#    logger.info('BL-based smoothing...')
#    MSs.run('BLsmooth.py -r -s 0.7 -i DATA -o DATA $pathMS', log='$nameMS_smooth.log', commandType='python')

for c in range(100):

    logger.info('== Start cycle: %s ==' % c)

    #logger.info('Remove bad timestamps...')
    MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    ####################################################
    # 1: find PA and remove it

    # Solve cal_SB.MS:DATA (only solve)
    if not os.path.exists('cal-pa-c0.h5'):
        logger.info('Solving PA...')
        MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal \
            sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solPA.log', commandType="DP3")

        lib_util.run_losoto(s, 'pa-c'+str(c), [ms+'/pa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

    #################################################
    # 2: find the FR and remve it
    
    # Beam correction DATA -> CORRECTED_DATA
    logger.info('PA correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c0.h5 cor.correction=polalign', \
            log='$nameMS_corPA2.log', commandType="DP3")

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType='DP3')
    
    # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=5)
    
    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Solving FR...')
    MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/fr.h5 sol.mode=diagonal \
             sol.smoothnessconstraint=1e6 sol.solint=3 sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solFR.log', commandType="DP3")
    
    lib_util.run_losoto(s, 'fr-c'+str(c), [ms+'/fr.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-fr.parset'])

   #####################################################
   # 3: find BANDPASS/IONO

    # Beam correction DATA -> CORRECTED_DATA
    logger.info('PA correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c0.h5 cor.correction=polalign', \
            log='$nameMS_corPA3.log', commandType="DP3")

    # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    if c == 0 and lofar_system == 'lba':
        MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_corBEAM3.log', commandType='DP3')
    else:
        MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_corBEAM3.log', commandType='DP3')
 
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR3.log', commandType="DP3")

    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Solving IONO...')
    MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/iono.h5 sol.mode=scalarcomplexgain \
                                    sol.smoothnessconstraint=1e6 sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solIONO3.log', commandType="DP3")

    lib_util.run_losoto(s, 'iono-c'+str(c), [ms+'/iono.h5' for ms in MSs.getListStr()], \
                        [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-plot-ph.parset'])

    # Correct all CORRECTED_DATA -> CORRECTED_DATA
    logger.info('IONO correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-iono-c'+str(c)+'.h5 cor.correction=phase000', \
                                log='$nameMS_corIONO3.log', commandType='DP3')

    # Solve MS:CORRECTED_DATA (only solve)
    logger.info('Solving BP...')
    MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/bp.h5 sol.mode=fulljones \
            sol.uvlambdarange='+str(nouseblrange)+' sol.smoothnessconstraint=2e6 sol.nchan=1 sol.solint=50', log='$nameMS_solBP3.log', commandType="DP3")
    
    lib_util.run_losoto(s, 'bp-c'+str(c), [ms+'/bp.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-plot-ph.parset'])


    # Correct BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('BP correction...')
    if c == 0 and lofar_system == 'lba':
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=True cor.parmdb=cal-bp-c'+str(c)+'.h5 cor.correction=fulljones \
                cor.soltab=\[amplitude000,phase000\]', \
                log='$nameMS_corBP3.log', commandType='DP3')
    else:
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-bp-c'+str(c)+'.h5 cor.correction=fulljones \
               cor.soltab=\[amplitude000,phase000\]', \
               log='$nameMS_corBP3.log', commandType='DP3')

    logger.info('Cleaning (cycle %02i)...' % c)
    imagename = 'img/img-c%02i' % c
        
    if MSs.isLBA():
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=1500, scale='2arcsec', \
                weight='briggs -1.0', niter=50000, no_update_model_required='', nmiter=50, mgain=0.4, \
                multiscale='', multiscale_scale_bias=0.6, multiscale_scales='0,5,10,20,40,80,160', \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirAlba.fits', \
                baseline_averaging=10, auto_threshold=1, \
                join_channels='', deconvolution_channels=8, fit_spectral_pol=4, channels_out=61)
   #     for modelfile in glob.glob(imagename+'*model*'):
   #         rev_reg(modelfile,'/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/virgohole.reg')

    if MSs.isHBA():
        lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1000, scale='2arcsec', \
                weight='briggs -0.2', niter=350, update_model_required='', mgain=0.5, \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirAphba.fits', \
                join_channels='', deconvolution_channels=5, fit_spectral_pol=5, channels_out=10) # use cont=True
        lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), cont=True, name=imagename, size=1000, scale='2arcsec', \
                weight='briggs -0.2', niter=5000000, no_update_model_required='', nmiter=100, mgain=0.4, \
                multiscale='', multiscale_scale_bias=0.7, \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirAhba.fits', \
                auto_threshold=1, \
                join_channels='', deconvolution_channels=5, fit_spectral_pol=5, channels_out=10)

    # TEST: rescale model
    #im = lib_img.Image(imagename+'-MFS-image.fits')
    #im.rescaleModel(f)

    logger.info('Predict (wsclean: %s)...' % imagename)
    s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channels-out 10 '+MSs.getStrWsclean(), \
          log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

#    # every 5 cycles: sub model and rescale model
#    if c%5 == 0 and c != 0:
#
#        logger.info('Sub model...')
#        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql1.log', commandType='general')
#
#        logger.info('Cleaning wide (cycle %i)...' % c)
#        imagename = 'img/imgsub-c'+str(c)
#        lib_util.run_wsclean(s, 'wscleanSUB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1000, scale='15arcsec', \
#                weight='briggs 0.4', taper_gaussian='100arcsec', niter=10000, no_update_model_required='', mgain=0.85, \
#                baseline_averaging=5, deconvolution_channels=4, \
#                auto_threshold=1, join_channels='', fit_spectral_pol=2, channels_out=16)
# 
#        #logger.info('Predict wide (wsclean)...')
#        #s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channelsout 32 '+MSs.getStrWsclean(), \
#        #      log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
#        #s.run(check = True)
#
#        #logger.info('Sub low-res model...')
#        #MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql2.log', commandType='general')

logger.info("Done.")
