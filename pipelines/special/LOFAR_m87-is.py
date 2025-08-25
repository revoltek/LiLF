#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob, re
import numpy as np

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
    os.makedirs('img')
    os.makedirs('mss')

MSs = lib_ms.AllMSs( sorted(glob.glob(data_dir+'/*MS')), s )
# copy data (avg to 1ch/sb and 10 sec)
nchan = 1 #int(MSs.getListObj()[0].getNchan()) # no avg in freq
timeint = MSs.getListObj()[0].getTimeInt()
avg_time = int(np.rint(2./timeint)) # average to 2 s
channels_out = 12

with w.if_todo('Averaging'):
    logger.info('Copy data...')
    MSs.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msout=mss/$nameMS.MS avg.freqstep=%i avg.timestep=%i' % (nchan, avg_time),\
                log='$nameMS_avg.log', commandType='DP3')

MSs = lib_ms.AllMSs( glob.glob('mss/*MS'), s )

########################################################   
# flag bad stations, and low-elev
with w.if_todo('flag'):
    logger.info('Flagging...')
    MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. steps=[ears,ant,uvmin,elev,count] \
                ant.baseline=\"'+bl2flag+'\" ears.type=preflagger ears.baseline=\"/(.*)HBA0&\\1HBA1/\"', \
                log='$nameMS_flag.log', commandType='DP3')

with w.if_todo('model'):
    im = lib_img.Image(model_dir+'/img-MFS-image.fits')
    #im.rescaleModel(f)
    n = len(glob.glob(model_dir+'/img-[0-9]*-model.fits'))
    logger.info('Predict (wsclean: %s - chan: %i)...' % (model_dir, n))
    s.add('wsclean -predict -name '+model_dir+'/img -j '+str(s.max_processors)+' -channels-out '+str(n)+' '+MSs.getStrWsclean(), \
              log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    s.run(check=True)

# init CS+RS correction DATA -> CORRECTED_DATA_DUTCH
with w.if_todo('initial_correction'):
    logger.info('PA correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=init-cal-pa.h5 cor.correction=polalign', \
        log='$nameMS_init_corPA.log', commandType="DP3")
    logger.info('IONO correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.parmdb=init-cal-iono.h5 cor.correction=phase000', \
        log='$nameMS_init_corIONO.log', commandType='DP3')
    logger.info('BP correction...')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msout.datacolumn=CORRECTED_DATA_DUTCH cor.parmdb=init-cal-bp.h5 cor.correction=fulljones \
        cor.soltab=\[amplitude000,phase000\]', log='$nameMS_init_corBP.log', commandType='DP3')

# TESTTESTTEST
# BL Smooth DATA -> DATA
#if MSs.isHBA:
#    logger.info('BL-based smoothing...')
#    MSs.run('BLsmooth.py -r -s 0.7 -i DATA -o DATA $pathMS', log='$nameMS_smooth.log', commandType='python')

for c in range(100):

    logger.info('== Start cycle: %s ==' % c)

    with w.if_todo('flag_c%02i' % c):
        logger.info('Remove bad timestamps...')
        MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    ####################################################
    # 1: find PA and remove it

    with w.if_todo('pa_c%02i' % c):
        # solve cal_SB.MS:CORRECTED_DATA_DUTCH (only solve)
        if not os.path.exists('cal-pa-c0.h5'):
            logger.info('solving PA...')
            MSs.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA_DUTCH sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal  sol.rotationdiagonalmode=diagonalphase sol.uvlambdarange='+str(nouseblrange)+ \
                    ' sol.antennaconstraint=[[CS001HBA0,CS001HBA1,CS002HBA0,CS002HBA1,CS003HBA0,CS003HBA1,CS004HBA0,CS004HBA1,CS005HBA0,CS005HBA1,CS006HBA0,CS006HBA1,CS007HBA0,CS007HBA1,CS011HBA0,CS011HBA1,CS013HBA0,CS013HBA1,CS017HBA0,CS017HBA1,CS021HBA0,CS021HBA1,CS024HBA0,CS024HBA1,CS026HBA0,CS026HBA1,CS028HBA0,CS028HBA1,CS030HBA0,CS030HBA1,CS031HBA0,CS031HBA1,CS032HBA0,CS032HBA1,CS101HBA0,CS101HBA1,CS103HBA0,CS103HBA1,CS201HBA0,CS201HBA1,CS301HBA0,CS301HBA1,CS302HBA0,CS302HBA1,CS401HBA0,CS401HBA1,CS501HBA0,CS501HBA1,RS106HBA,RS205HBA,RS208HBA,RS210HBA,RS305HBA,RS306HBA,RS307HBA,RS310HBA,RS406HBA,RS407HBA,RS409HBA,RS503HBA,RS508HBA,RS509HBA]]', \
                    log='$nameMS_solPA.log', commandType="DP3")

            lib_util.run_losoto(s, 'pa-c'+str(c), [ms+'/pa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

    #################################################
    # 2: find the FR and remve it
    
    with w.if_todo('fr_c%02i' % c):
        # PA correction CORRECTED_DATA_DUTCH -> CORRECTED_DATA
        logger.info('PA correction...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA_DUTCH cor.parmdb=cal-pa-c0.h5 cor.correction=polalign', \
                    log='$nameMS_corPA2.log', commandType="DP3")
    
        # Beam correction CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Beam correction...')
        #MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType='DP3')
            
        # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Converting to circular...')
        MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=5)
            
        # solve cal_SB.MS:CORRECTED_DATA (only solve)
        logger.info('solving FR...')
        MSs.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/fr.h5 sol.mode=diagonal sol.datause=dual \
                    sol.smoothnessconstraint=1e6 sol.solint=3 sol.uvlambdarange='+str(nouseblrange)+ \
                    ' sol.antennaconstraint=[[CS001HBA0,CS001HBA1,CS002HBA0,CS002HBA1,CS003HBA0,CS003HBA1,CS004HBA0,CS004HBA1,CS005HBA0,CS005HBA1,CS006HBA0,CS006HBA1,CS007HBA0,CS007HBA1,CS011HBA0,CS011HBA1,CS013HBA0,CS013HBA1,CS017HBA0,CS017HBA1,CS021HBA0,CS021HBA1,CS024HBA0,CS024HBA1,CS026HBA0,CS026HBA1,CS028HBA0,CS028HBA1,CS030HBA0,CS030HBA1,CS031HBA0,CS031HBA1,CS032HBA0,CS032HBA1,CS101HBA0,CS101HBA1,CS103HBA0,CS103HBA1,CS201HBA0,CS201HBA1,CS301HBA0,CS301HBA1,CS302HBA0,CS302HBA1,CS401HBA0,CS401HBA1,CS501HBA0,CS501HBA1,RS106HBA,RS205HBA,RS208HBA,RS210HBA,RS305HBA,RS306HBA,RS307HBA,RS310HBA,RS406HBA,RS407HBA,RS409HBA,RS503HBA,RS508HBA,RS509HBA]]', \
                    log='$nameMS_solFR.log', commandType="DP3")
            
        lib_util.run_losoto(s, 'fr-c'+str(c), [ms+'/fr.h5' for ms in MSs.getListStr()], \
                    [parset_dir + '/losoto-fr.parset'])

    #####################################################
    # 3: find BANDPASS/IONO

    with w.if_todo('iono_c%02i' % c):
        # PA correction CORRECTED_DATA_DUTCH -> CORRECTED_DATA
        logger.info('PA correction...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA_DUTCH cor.parmdb=cal-pa-c0.h5 cor.correction=polalign', \
                log='$nameMS_corPA3.log', commandType="DP3")
    
        # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
        #logger.info('Beam correction...')
        #if c == 0 and MSs.isLBA:
        #    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights='+str(updateweights), log='$nameMS_corBEAM3.log', commandType='DP3')
        #else:
        #    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_corBEAM3.log', commandType='DP3')
     
        # Correct FR CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
                log='$nameMS_corFR3.log', commandType="DP3")
    
        # solve cal_SB.MS:CORRECTED_DATA (only solve)
        logger.info('solving IONO...')
        MSs.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/iono.h5 sol.mode=scalarphase sol.datause=single sol.smoothnessconstraint=1e6 sol.uvlambdarange='+str(nouseblrange)+ \
                ' sol.antennaconstraint=[[CS001HBA0,CS001HBA1,CS002HBA0,CS002HBA1,CS003HBA0,CS003HBA1,CS004HBA0,CS004HBA1,CS005HBA0,CS005HBA1,CS006HBA0,CS006HBA1,CS007HBA0,CS007HBA1,CS011HBA0,CS011HBA1,CS013HBA0,CS013HBA1,CS017HBA0,CS017HBA1,CS021HBA0,CS021HBA1,CS024HBA0,CS024HBA1,CS026HBA0,CS026HBA1,CS028HBA0,CS028HBA1,CS030HBA0,CS030HBA1,CS031HBA0,CS031HBA1,CS032HBA0,CS032HBA1,CS101HBA0,CS101HBA1,CS103HBA0,CS103HBA1,CS201HBA0,CS201HBA1,CS301HBA0,CS301HBA1,CS302HBA0,CS302HBA1,CS401HBA0,CS401HBA1,CS501HBA0,CS501HBA1,RS106HBA,RS205HBA,RS208HBA,RS210HBA,RS305HBA,RS306HBA,RS307HBA,RS310HBA,RS406HBA,RS407HBA,RS409HBA,RS503HBA,RS508HBA,RS509HBA]]', \
                log='$nameMS_solIONO.log', commandType="DP3")
    
        lib_util.run_losoto(s, 'iono-c'+str(c), [ms+'/iono.h5' for ms in MSs.getListStr()], [parset_dir+'/losoto-plot-ph-nopol.parset'])
    
        # Correct all CORRECTED_DATA -> CORRECTED_DATA
        logger.info('IONO correction...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-iono-c'+str(c)+'.h5 cor.correction=phase000', \
                                    log='$nameMS_corIONO3.log', commandType='DP3')
    
        # solve MS:CORRECTED_DATA (only solve)
        logger.info('solving BP...')
        MSs.run('DP3 ' + parset_dir + '/DP3-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/bp.h5 sol.mode=fulljones \
                sol.uvlambdarange='+str(nouseblrange)+' sol.nchan='+str(channels_out)+' sol.solint=10 \
                sol.antennaconstraint=[[CS001HBA0,CS001HBA1,CS002HBA0,CS002HBA1,CS003HBA0,CS003HBA1,CS004HBA0,CS004HBA1,CS005HBA0,CS005HBA1,CS006HBA0,CS006HBA1,CS007HBA0,CS007HBA1,CS011HBA0,CS011HBA1,CS013HBA0,CS013HBA1,CS017HBA0,CS017HBA1,CS021HBA0,CS021HBA1,CS024HBA0,CS024HBA1,CS026HBA0,CS026HBA1,CS028HBA0,CS028HBA1,CS030HBA0,CS030HBA1,CS031HBA0,CS031HBA1,CS032HBA0,CS032HBA1,CS101HBA0,CS101HBA1,CS103HBA0,CS103HBA1,CS201HBA0,CS201HBA1,CS301HBA0,CS301HBA1,CS302HBA0,CS302HBA1,CS401HBA0,CS401HBA1,CS501HBA0,CS501HBA1,RS106HBA,RS205HBA,RS208HBA,RS210HBA,RS305HBA,RS306HBA,RS307HBA,RS310HBA,RS406HBA,RS407HBA,RS409HBA,RS503HBA,RS508HBA,RS509HBA]]', \
                log='$nameMS_solBP3.log', commandType="DP3")
        
        lib_util.run_losoto(s, 'bp-c'+str(c), [ms+'/bp.h5' for ms in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-amp-nopol.parset',parset_dir+'/losoto-plot-ph-nopol.parset'])
    
        # Correct BP CORRECTED_DATA -> CORRECTED_DATA
        logger.info('BP correction...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-bp-c'+str(c)+'.h5 cor.correction=fulljones \
                   cor.soltab=\[amplitude000,phase000\]', log='$nameMS_corBP3.log', commandType='DP3')

    with w.if_todo('image_c%02i' % c):
        logger.info('Cleaning (cycle %02i)...' % c)
        imagename = 'img/img-c%02i' % c
        
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, no_update_model_required='', baseline_averaging=6, minuv_l=175, \
                    reorder='', parallel_reordering=4, gridder='wgridder', size=1600, scale='0.1arcsec', padding=1.6, \
                    weight='briggs -1', niter=1000000, nmiter=100, mgain=0.85, \
                    multiscale='', multiscale_scales='0,5,10,20,40,80,160,320', \
                    auto_threshold=1, \
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
