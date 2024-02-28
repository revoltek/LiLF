#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob, re
import numpy as np

patch = 'VirA'
nouseblrange = ''
f = lambda nu: 10**(2.4466 - 0.8116 * ((np.log10(nu/1.e9))**1) - 0.0483 * ((np.log10(nu/1.e9))**2) ) # PB17

skymodel = '/home/baq1889/scripts/model/A-team_4_CC.skydb'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-m87')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-m87.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_m87'])))
parset_dir = parset.get('LOFAR_m87','parset_dir')
data_dir = parset.get('LOFAR_m87','data_dir')
updateweights = parset.getboolean('LOFAR_m87','updateweights')
skipmodel = parset.getboolean('LOFAR_m87','skipmodel')
model_dir = parset.get('LOFAR_m87','model_dir')
bl2flag = parset.get('flag','stations')

##########################################################
with w.if_todo('clean'):
    logger.info('Cleaning...')
    lib_util.check_rm('cal*h5')
    lib_util.check_rm('plots*')
    lib_util.check_rm('img')
    lib_util.check_rm('mss')
    lib_util.check_rm('mssIS')
    os.makedirs('img')
    os.makedirs('mss')
    os.makedirs('mssIS')

MSs = lib_ms.AllMSs( sorted(glob.glob(data_dir+'/*MS')), s )
# copy data (avg to 1ch/sb and 10 sec)
nchan = 1 #int(MSs.getListObj()[0].getNchan()) # no avg in freq
timeint = MSs.getListObj()[0].getTimeInt()
avg_time = int(np.rint(2./timeint)) # average to 2 s

with w.if_todo('AveragingDutch...'):
    logger.info('Copy data dutch...')
    mss_bkp = sorted(glob.glob(data_dir+'/*MS'))
    MSs.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msin.baselines=[CR]S*& \
            msout=mss/$nameMS.MS avg.freqstep=1 avg.timestep=2' % (nchan, avg_time),\
                log='$nameMS_avg.log', commandType='DP3')

with w.if_todo('AveragingIS'):
    logger.info('Copy data IS...')
    mss_bkp = sorted(glob.glob(data_dir+'/*MS'))
    MSs.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msin.baselines=*& msin.startchan=1200 msin.nchan=2000 \
             msout=mssIS/$nameMS.MS avg.freqstep=1 avg.timestep=1' % (nchan, avg_time),\
            log='$nameMS_avg.log', commandType='DP3')


MSs_dutch = lib_ms.AllMSs( glob.glob('mss/*MS'), s )
MSs_IS = lib_ms.AllMSs( glob.glob('mssIS/*MS'), s )

if MSs.isLBA and not MSs.hasIS:
    channels_out = 20
else:
    channels_out = 12

########################################################   
# flag bad stations, and low-elev
with w.if_todo('flag'):
    logger.info('Flagging...')
    if MSs.isLBA:
        MSs_dutch.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. steps=[ant,uvmin,elev,count] ant.baseline=\"'+bl2flag+'\"', \
                log='$nameMS_flag.log', commandType='DP3')
        MSs_IS.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. steps=[ant,uvmin,elev,count] ant.baseline=\"'+bl2flag+'\"', \
                log='$nameMS_flag.log', commandType='DP3')
    else:
        MSs_dutch.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. steps=[ears,ant,uvmin,elev,count] \
                ant.baseline=\"'+bl2flag+'\" ears.type=preflagger ears.baseline=\"/(.*)HBA0&\\1HBA1/\"', \
                log='$nameMS_flag.log', commandType='DP3')
        MSs_IS.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. steps=[ears,ant,uvmin,elev,count] \
                ant.baseline=\"'+bl2flag+'\" ears.type=preflagger ears.baseline=\"/(.*)HBA0&\\1HBA1/\"', \
                log='$nameMS_flag.log', commandType='DP3')
        
    MSs_dutch.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
    MSs_IS.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

if model_dir == '' and MSs.isHBA: model_dir = '/home/fdg/scripts/model/AteamHBA/'+patch
if model_dir == '' and MSs.isLBA: model_dir = '/home/fdg/scripts/model/AteamLBA/'+patch

####################################################################
with w.if_todo('model'):
    if not skipmodel and os.path.exists(model_dir+'/img-MFS-model.fits'):
        im = lib_img.Image(model_dir+'/img-MFS-image.fits')
        im.rescaleModel(f)
        n = len(glob.glob(model_dir+'/img-[0-9]*-model.fits'))
        logger.info('Predict Dutch (wsclean: %s - chan: %i)...' % (model_dir, n))
        s.add('wsclean -predict -name '+model_dir+'/img -j '+str(s.max_processors)+' -channels-out '+str(n)+' '+MSs_dutch.getStrWsclean(), \
              log='wscleanPRE-init.log', commandType='wsclean', processors='max')
        s.run(check=True)
        logger.info('Predict IS (wsclean: %s - chan: %i)...' % (model_dir, n))
        s.add('wsclean -predict -name '+model_dir+'/img -j '+str(s.max_processors)+' -channels-out '+str(n)+' '+MSs_IS.getStrWsclean(), \
              log='wscleanPRE-init.log', commandType='wsclean', processors='max')
        s.run(check=True)
    elif not skipmodel:
        logger.info('Predict (DP3)...')
        MSs.run('DP3 '+parset_dir+'/DP3-predict.parset msin=$pathMS pre.sourcedb='+skymodel+' pre.sources='+patch, log='$nameMS_pre.log', commandType='DP3')

########################################################################
with w.if_todo('phaseupcore'):
    MSs_dutch.run(
        "DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS msin.datacolumn=DATA corrbeam.updateweights=True',
        log='$nameMS_beam.log', commandType="DP3")

    # Smooth data concat_core.MS:DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs_dutch.run(
        f'BLsmooth.py -r -q -c 1 -n 8 -f {1e-3 if MSs.isLBA else .2e-3} -i DATA -o SMOOTHED_DATA $pathMS',
        log='$nameMS_smooth3.log', commandType='python', maxThreads=8)

    # Solve concat_core.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating fullj cores...')
    MSs_dutch.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
            sol.h5parm=$pathMS/fullj.h5 sol.mode=fulljones \
            sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=1e6', \
                        log='$nameMS_solBP.log', commandType="DP3")
    lib_util.run_losoto(s, 'core_fullj', [ms + '/fullj.h5' for ms in MSs_dutch.getListStr()],
                        [parset_dir + '/losoto-flag.parset', parset_dir + '/losoto-plot-fullj.parset'])

    # Correct cores concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA (for cores)
    logger.info('FullJ correction...')
    MSs_IS.run("DP3 " + parset_dir + '/DP3-beam.parset msin=$pathMS msin.datacolumn=DATA corrbeam.updateweights=True \
                       corrbeam.noapplystations=[' + ','.join(RSISlist) + ']', log='$nameMS_beam.log',
                       commandType="DP3")
    MSs_IS.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-core_fullj.h5 \
                cor.correction=fulljones cor.soltab=[amplitude000,phase000] cor.updateweights=True',
                       log='$nameMS_corFullJ.log', commandType="DP3")

    # Phasing up the cose stations
    logger.info('Phasing up Core Stations...')
    lib_util.check_rm('concat_all-phaseup.MS')
    MSs_IS.run('DP3 ' + parset_dir + '/DP3-phaseup.parset msin=$pathMS msout=concat_all-phaseup.MS',
                       log='$nameMS_phaseup.log', commandType="DP3")

### DONE

for c in range(100):

    logger.info('== Start cycle: %s ==' % c)
    
    with w.if_todo('cal_c%02i' % c):
        logger.info('Calibrating...')
        MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS\
                sol.h5parm=$pathMS/cal.h5 sol.mode=fulljones \
                sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=.5e6', \
                log='$nameMS_sol.log', commandType="DP3")

        lib_util.run_losoto(s, '%i' % c, [ms+'/cal.h5' for ms in MSs.getListStr()],
                [parset_dir+'/losoto-plot-fullj.parset'])
                
        logger.info('Correcting...')
        if c==1:
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.parmdb=cal-'+str(c)+'.h5 \
                cor.correction=fulljones cor.soltab=[amplitude000,phase000] cor.updateweights=True', log='$nameMS_cor.log', commandType="DP3")
        else:
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.parmdb=cal-'+str(c)+'.h5 \
                cor.correction=fulljones cor.soltab=[amplitude000,phase000] cor.updateweights=False', log='$nameMS_cor.log', commandType="DP3")

        
    """
    ####################################################
    # 1: find PA and remove it

    with w.if_todo('pa_c%02i' % c):
        # Solve cal_SB.MS:DATA (only solve)
        if not os.path.exists('cal-pa-c0.h5'):
            logger.info('Solving PA...')
            MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal \
                sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solPA.log', commandType="DP3")

            lib_util.run_losoto(s, 'pa-c'+str(c), [ms+'/pa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

    #################################################
    # 2: find the FR and remve it
    
    if MSs.isLBA:
        with w.if_todo('fr_c%02i' % c):
            # PA correction DATA -> CORRECTED_DATA
            logger.info('PA correction...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c0.h5 cor.correction=polalign', \
                    log='$nameMS_corPA2.log', commandType="DP3")
    
            # Beam correction CORRECTED_DATA -> CORRECTED_DATA
            #logger.info('Beam correction...')
            #MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType='DP3')
            
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

    with w.if_todo('sol3_c%02i' % c):
        # PA correction DATA -> CORRECTED_DATA
        logger.info('PA correction...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c0.h5 cor.correction=polalign', \
                log='$nameMS_corPA3.log', commandType="DP3")
    
        # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Beam correction...')
        #if c == 0 and MSs.isLBA:
        #    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights='+str(updateweights), log='$nameMS_corBEAM3.log', commandType='DP3')
        #else:
        #    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_corBEAM3.log', commandType='DP3')
     
        if MSs.isLBA:
            # Correct FR CORRECTED_DATA -> CORRECTED_DATA
            logger.info('Faraday rotation correction...')
            MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
                    log='$nameMS_corFR3.log', commandType="DP3")
    
        # Solve cal_SB.MS:CORRECTED_DATA (only solve)
        logger.info('Solving IONO...')
        MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/iono.h5 sol.mode=scalarcomplexgain \
                                        sol.smoothnessconstraint=1e6 sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solIONO3.log', commandType="DP3")
    
        lib_util.run_losoto(s, 'iono-c'+str(c), [ms+'/iono.h5' for ms in MSs.getListStr()], \
                            [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-plot-amp-nopol.parset',parset_dir+'/losoto-plot-ph-nopol.parset'])
    
        # Correct all CORRECTED_DATA -> CORRECTED_DATA
        logger.info('IONO correction...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-iono-c'+str(c)+'.h5 cor.correction=phase000', \
                                    log='$nameMS_corIONO3.log', commandType='DP3')
    
        # Solve MS:CORRECTED_DATA (only solve)
        logger.info('Solving BP...')
        if MSs.isLBA:
            MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/bp.h5 sol.mode=fulljones \
                sol.uvlambdarange='+str(nouseblrange)+' sol.smoothnessconstraint=2e6 sol.nchan=1 sol.solint=50', log='$nameMS_solBP3.log', commandType="DP3")
        else:
            MSs.run('DP3 ' + parset_dir + '/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/bp.h5 sol.mode=fulljones \
                sol.uvlambdarange='+str(nouseblrange)+' sol.nchan='+str(channels_out)+' sol.solint=10', log='$nameMS_solBP3.log', commandType="DP3")
        
        lib_util.run_losoto(s, 'bp-c'+str(c), [ms+'/bp.h5' for ms in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-amp-nopol.parset',parset_dir+'/losoto-plot-ph-nopol.parset'])
    
    
        # Correct BP CORRECTED_DATA -> CORRECTED_DATA
        logger.info('BP correction...')
        if c == 0 and MSs.isLBA:
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights='+str(updateweights)+' cor.parmdb=cal-bp-c'+str(c)+'.h5 cor.correction=fulljones \
                    cor.soltab=\[amplitude000,phase000\]', \
                    log='$nameMS_corBP3.log', commandType='DP3')
        else:
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-bp-c'+str(c)+'.h5 cor.correction=fulljones \
                   cor.soltab=\[amplitude000,phase000\]', \
                   log='$nameMS_corBP3.log', commandType='DP3')
    """

    with w.if_todo('image_c%02i' % c):
        logger.info('Cleaning (cycle %02i)...' % c)
        imagename = 'img/img-c%02i' % c
        
        if MSs.isLBA and not MSs.hasIS:
            lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, no_update_model_required='', baseline_averaging=8, parallel_gridding=4,\
                    reorder='', parallel_reordering=4, gridder='wgridder', size=1500, scale='2arcsec', padding=1.2, \
                    weight='briggs -1.0', niter=50000, nmiter=50, mgain=0.4, \
                    multiscale='', multiscale_scale_bias=0.6, \
                    fits_mask='/home/baq1889/LiLF/parsets/LOFAR_ateam/masks/VirAlba.fits', auto_threshold=1, \
                    join_channels='', channels_out=channels_out)

        if MSs.isLBA and  MSs.hasIS:
            lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, no_update_model_required='', baseline_averaging=8, parallel_gridding=4,\
                    reorder='', parallel_reordering=4, gridder='wgridder', size=2000, scale='1arcsec', padding=1.2, \
                    weight='briggs 0', taper_gaussian='0.75arcsec', niter=15000, nmiter=50, mgain=0.4, \
                    multiscale='', multiscale_scale_bias=0.6, \
                    fits_mask='/home/baq1889/LiLF/parsets/LOFAR_ateam/masks/VirAlbaIS.fits', auto_threshold=1, \
                    join_channels='', channels_out=channels_out)
    
        if MSs.isHBA:
            #lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1000, scale='2arcsec', \
            #        weight='briggs -0.2', niter=350, update_model_required='', mgain=0.5, \
            #        fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirAphba.fits', \
            #        join_channels='', deconvolution_channels=5, fit_spectral_pol=5, channels_out=channels_out) # use cont=True
            lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, no_update_model_required='', baseline_averaging=6, minuv_l=175, \
                    reorder='', parallel_reordering=4, wgridder='wgridder', size=1600, scale='1arcsec', padding=1.6, \
                    weight='briggs -1', niter=1000000, nmiter=100, mgain=0.85, \
                    multiscale='', multiscale_scales='0,20,40,80,160,320', \
                    fits_mask='/home/baq1889/LiLF/parsets/LOFAR_ateam/masks/VirAhba.fits', auto_threshold=1, \
                    join_channels='',  channels_out=channels_out)

    with w.if_todo('predict_c%02i' % c):
        # TEST: rescale model
        #im = lib_img.Image(imagename+'-MFS-image.fits')
        #im.rescaleModel(f)
        logger.info('Predict (wsclean: %s)...' % imagename)
        s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channels-out '+str(channels_out)+' '+MSs.getStrWsclean(), \
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
