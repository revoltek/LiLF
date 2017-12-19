#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt

parset_dir = '/home/fdg/scripts/autocal/AteamLBA/parset_ateam/'

# Temporary!
if 'VirA2013' in os.getcwd():
    patch = 'VirA'
    datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts2013-bkp'
    bl2flag = 'CS013LBA\;CS031LBA'
    blrange = '[0,1e30]'
elif 'VirA2015' in os.getcwd():
    patch = 'VirA'
    datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts2015-bkp'
    bl2flag = 'CS017LBA\;RS407LBA'
    blrange = '[0,1e30]'
elif 'VirA2017' in os.getcwd():
    patch = 'VirA'
    datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts2017-bkp'
    bl2flag = ''
    blrange = '[0,1e30]'
elif 'TauA' in os.getcwd():
    patch = 'TauA'
    datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/TauA/tgts-bkp'
    bl2flag = 'RS310LBA\;RS210LBA\;RS407LBA\;RS409LBA'
    blrange = '[0,1000,5000,1e30]'
elif 'CasA' in os.getcwd():
    patch = 'CasA'
    datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CasA/tgts1-bkp'
    #datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CasA/tgts2-bkp'
    bl2flag = 'CS031LBA'
    blrange = '[0,30000]'
elif 'CygA' in os.getcwd():
    patch = 'CygA'
    datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CygA/tgts1-bkp'
    #datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CygA/tgts2-bkp'
    bl2flag = 'CS031LBA'
    blrange = '[0,30000]'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
lib_log.set_logger('pipeline-ateam.logger')
logger = lib_log.logger
lib_util.check_rm('logs')
s = lib_util.Scheduler(dry = False)
lib_util.check_rm('img')
os.makedirs('img')
mss = sorted(glob.glob(datadir+'/*MS'))
MSs = lib_ms.AllMSs( mss[len(mss)*2/5:len(mss)*4/5], s ) # Use only 1/2 of the SBs

# copy data
logger.info('Copy data...')
for MS in MSs.getListObj():
    MS.move(MS.nameMS+'.MS', keepOrig=True)

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

#############################################################   
## flag bad stations, and low-elev
logger.info('Flagging...')
MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS msout=. flag1.baseline='+bl2flag+' msin.datacolumn=DATA', \
            log='$nameMS_flag.log', commandType='DPPP')
 
# predict to save time ms:MODEL_DATA
if os.path.exists('/home/fdg/scripts/model/AteamLBA/'+patch+'/wideM-MFS-model.fits'):
    logger.info('Predict (wsclean)...')
    s.add('wsclean -predict -name /home/fdg/scripts/model/AteamLBA/'+patch+'/wideM -mem 90 -j '+str(s.max_processors)+' -channelsout 15 '+MSs.getStrWsclean(), \
          log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    s.run(check=True)
else:
    logger.info('Predict (DPPP)...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=/home/fdg/scripts/model/A-team_4_CC.skydb pre.sources='+patch, log='$nameMS_pre.log', commandType='DPPP')

for c in xrange(10):

    #################################################
    # 1: find the FR and remve it
    
    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS', log='$nameMS_beam.log', commandType='DPPP')
    
    # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python')
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType='python')
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for MS in MSs.getListStr():
        lib_util.check_rm(MS+'/fr.h5')
    MSs.run('DPPP '+parset_dir+'/DPPP-sol.parset msin=$pathMS filter.blrange='+blrange+' sol.parmdb=$pathMS/fr.h5', log='$nameMS_sol1.log', commandType='DPPP')
    
    lib_util.run_losoto(s, 'fr-c'+str(c), [ms+'/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
    
    #####################################################
    # 2: find BANDPASS

    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType="DPPP")
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DPPP")
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2.log', commandType ='python', maxThreads=20)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for MS in MSs.getListStr():
        lib_util.check_rm(MS+'/amp.h5')
    MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS filter.blrange='+blrange+' sol.parmdb=$pathMS/amp.h5', log='$nameMS_sol2.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'amp', [ms+'/amp.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-align.parset'])

    #################################################
    # 3: recalibrate without FR

    # Correct ampBP DATA (beam corrected) -> CORRECTED_DATA
    logger.info('Cross delay+ampBP correction...')
    if c == 0:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin='+ms+' msin.datacolumn=DATA cor.updateweights=True cor.parmdb=$pathMS/amp.h5 cor.correction=amplitudeSmooth000', \
                log=ms+'$nameMS_corAMP.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin='+ms+' msin.datacolumn=DATA cor.updateweights=False cor.parmdb=$pathMS/amp.h5 cor.correction=amplitudeSmooth000', \
                log=ms+'$nameMS_corAMP.log', commandType='DPPP')
 
    # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    if c == 0:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA corrbeam.updateweights=True', log='$nameMS_beam2.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA corrbeam.updateweights=False', log='$nameMS_beam2.log', commandType='DPPP')
       
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR.log', commandType='DPPP')
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3.log', commandType='python')
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for MS in MSs.getListStr():
        lib_util.check_rm(MS+'/iono.h5')
    MSs.run('DPPP '+parset_dir+'/DPPP-sol.parset msin=$pathMS filter.blrange='+blrange+' sol.parmdb=$pathMS/iono.h5', log='$nameMS_sol3.log', commandType='DPPP')
    
    run_losoto(s, 'final-c'+str(c), mss, [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-ph.parset'])
    
    # Correct all CORRECTED_DATA (beam, CD, FR corrected) -> CORRECTED_DATA
    logger.info('Amp/ph correction...')
    for ms in mss:
        if c == 0:
            s.add('DPPP '+parset_dir+'/DPPP-cor.parset msin='+ms+' cor.updateweights=True cor.parmdb='+ms+'/instrument cor.correction=gain', log=ms+'_corG.log', commandType='DPPP')
        else:
            # update weight only first time, it should be at first order correct
            s.add('DPPP '+parset_dir+'/DPPP-cor.parset msin='+ms+' cor.updateweights=False cor.parmdb='+ms+'/instrument cor.correction=gain', log=ms+'_corG.log', commandType='DPPP')
    s.run(check=True)
    
    # briggs: -1.2 for virgo
    logger.info('Cleaning (cycle %i)...' % c)
    imagename = 'img/wideM-c'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 1700 1700 -trim 1500 1500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 1.5 \
            -scale 2arcsec -weight briggs -1.5 -niter 100000 -no-update-model-required -mgain 0.7 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,4,8,16,32 -auto-mask 5\
            -pol I -joinchannels -fit-spectral-pol 3 -channelsout 15 -threshold 0.005 '+' '.join(mss), \
            log='wscleanB-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

    logger.info('Predict (ft)...')
    s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 15 '+' '.join(mss), \
            log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

#    logger.info('Sub model...')
#    for ms in mss:
#        s.add('taql "update '+ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log=ms+'_taql.log', commandType='general')
#    s.run(check=True)
#
#    logger.info('Cleaning sub (cycle %i)...' % c)
#    imagename = 'img/wideMsub-c'+str(c)
#    s.add('wsclean -reorder -name ' + imagename + ' -size 500 500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 1.5 \
#            -scale 10arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -mgain 0.7 -taper-gaussian 45arcsec \
#            -pol I -joinchannels -fit-spectral-pol 3 -channelsout 15 '+' '.join(mss), \
#            log='wscleanB-c'+str(c)+'.log', commandType='wsclean', processors='max')
#    s.run(check=True)

logger.info("Done.")
