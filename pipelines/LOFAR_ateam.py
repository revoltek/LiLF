#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt

# Temporary!
if 'Vir' in os.getcwd():
    patch = 'VirA'
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

skymodel = '/home/fdg/scripts/model/A-team_4_CC.skydb'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
lib_log.set_logger('pipeline-ateam.logger')
logger = lib_log.logger
s = lib_util.Scheduler(dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('ateam','parset_dir')
bl2flag = parset.get('flag','stations')
data_dir = '../tgts-bkp/'

##########################################################
lib_util.check_rm('img')
os.makedirs('img')
mss = sorted(glob.glob(data_dir+'/*MS'))
MSs = lib_ms.AllMSs( mss[len(mss)*2/5:len(mss)*4/5], s ) # Use only 1/2 of the SBs

# copy data
logger.info('Copy data...')
for MS in MSs.getListObj():
    MS.move(MS.nameMS+'.MS', keepOrig=True)

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

###########################################################   
# flag bad stations, and low-elev
logger.info('Flagging...')
MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS msout=. ant.baseline='+bl2flag, \
            log='$nameMS_flag.log', commandType='DPPP')

# predict to save time MODEL_DATA
#if os.path.exists('/home/fdg/scripts/model/AteamLBA/'+patch+'/wideM-MFS-model.fits'):
#    logger.info('Predict (wsclean)...')
#    s.add('wsclean -predict -name /home/fdg/scripts/model/AteamLBA/'+patch+'/wideM -mem 90 -j '+str(s.max_processors)+' -channelsout 15 '+MSs.getStrWsclean(), \
#          log='wscleanPRE-init.log', commandType='wsclean', processors='max')
#    s.run(check=True)
#else:

logger.info('Predict (DPPP)...')
MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel+' pre.sources='+patch, log='$nameMS_pre.log', commandType='DPPP')

for c in xrange(10):

    ####################################################
    # 1: find PA and remove it

    # Smooth data DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType='python', maxThreads=6)

    logger.info('Calibrating...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal', log='$nameMS_solPA.log', commandType="DPPP")

    lib_util.run_losoto(s, 'pa-c'+str(c), [ms+'/pa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

    #################################################
    # 1: find the FR and remve it
    
    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c'+str(c)+'.h5 cor.correction=polalign', \
            log='$nameMS_corPA2.log', commandType="DPPP")

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType='DPPP')
    
    # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2.log', commandType='python', maxThreads=6)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-soldd.parset msin=$pathMS sol.parmdb=$pathMS/fr.h5 sol.mode=diagonal', log='$nameMS_fr.log', commandType='DPPP')
    
    lib_util.run_losoto(s, 'fr-c'+str(c), [ms+'/fr.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-fr.parset'])
    
    #####################################################
    # 2: find BANDPASS

    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c'+str(c)+'.h5 cor.correction=polalign', \
            log='$nameMS_corPA3.log', commandType="DPPP")

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-beam.parset msin=$pathMS', log='$nameMS_beam3.log', commandType="DPPP")
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR3.log', commandType="DPPP")
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3.log', commandType ='python', maxThreads=6)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.parmdb=$pathMS/amp.h5 sol.mode=diagonal', log='$nameMS_amp.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'amp', [ms+'/amp.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-bp.parset'])

    #################################################
    # 3: apply all

    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c'+str(c)+'.h5 cor.correction=polalign', \
            log='$nameMS_corPA4.log', commandType="DPPP")

    # Correct BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Cross bp correction...')
    if c == 0:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin='+ms+' cor.updateweights=True cor.parmdb=cal-amp-c'+str(c)+'.h5 cor.correction=amplitudeSmooth000', \
                log=ms+'$nameMS_corAMP4.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin='+ms+' cor.updateweights=False cor.parmdb=cal-amp-c'+str(c)+'.h5 cor.correction=amplitudeSmooth000', \
                log=ms+'$nameMS_corAMP4.log', commandType='DPPP')
 
    # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    if c == 0:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam4.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_beam4.log', commandType='DPPP')
       
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR4.log', commandType='DPPP')
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth4.log', commandType='python', maxThreads=6)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-soldd.parset msin=$pathMS sol.parmdb=$pathMS/iono.h5 sol.mode=diagonal', log='$nameMS_iono.log', commandType='DPPP')
    
    run_losoto(s, 'iono-c'+str(c), mss, [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-ph.parset'])
    
    # Correct all CORRECTED_DATA (PA, beam, BP, FR corrected) -> CORRECTED_DATA
    logger.info('Amp/ph correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.steps=[ph,amp] cor.ph.parmdb=cal-iono-c'+str(c)+'.h5 cor.amp.parmdb=cal-iono.h5 \
                    cor.ph.correction=phaseOrig000 cor.amp.correction=amplitude000 cor.amp.updateweights=False', log='$nameMS_corG4.log', commandType="DPPP")
    
    # briggs: -1.2 for virgo
    logger.info('Cleaning (cycle %i)...' % c)
    imagename = 'img/wideM-c'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 1700 1700 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 1.5 \
            -scale 2arcsec -weight briggs -1.0 -niter 100000 -mgain 0.7 -minuv-l 100 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,4,8,16,32 -auto-mask 5\
            -pol I -join-channels -fit-spectral-pol 3 -channels-out 15 -auto-threshold 0.005 '+MSs.getStrWsclean(), \
            log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors = 'max')

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
