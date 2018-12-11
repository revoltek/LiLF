#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt

if 'Vir' in os.getcwd():
    patch = 'VirA'
    nouseblrange = '[0..30]'
    f = lambda nu: 1226. * 10**(-0.79 * (np.log10(nu/150.e6))**1)
elif 'Tau' in os.getcwd():
    patch = 'TauA'
    nouseblrange = '[500..5000]'
    f = lambda nu: 1838. * 10**(-0.299 * (np.log10(nu/150.e6))**1)
elif 'Cas' in os.getcwd():
    patch = 'CasA'
    nouseblrange = '[15000..1e30]'
    f = lambda nu: 11733. * 10**(-0.77 * (np.log10(nu/150.e6))**1)
elif 'Cyg' in os.getcwd():
    patch = 'CygA'
    nouseblrange = '[15000..1e30]'
    f = lambda nu: 10690. * 10**(-0.67 * (np.log10(nu/150.e6))**1) * 10**(-0.204 * (np.log10(nu/150.e6))**2) * 10**(-0.021 * (np.log10(nu/150.e6))**3)

skymodel = '/home/fdg/scripts/model/A-team_4_CC.skydb'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
lib_log.Logger('pipeline-ateam.logger')
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
MSs = lib_ms.AllMSs( mss, s )

# HBA/LBA
if min(MSs.getFreqs()) < 80.e6: hba = False
else: hba = True

# wsclean temp-dir
temp_dir = '.'
if s.get_cluster() == 'Hamburg_fat': temp_dir = '/localwork.ssd'

# copy data (avg to 1ch/sb and 10 sec)
nchan = MSs.getListObj()[0].getNchan()
timeint = MSs.getListObj()[0].getTimeInt()
avg_time = int(np.rint(10./timeint))

logger.info('Copy data...')
for MS in MSs.getListObj():
    if os.path.exists(MS.nameMS+'.MS'): continue
    s.add('DPPP '+parset_dir+'/DPPP-avg.parset msin='+MS.pathMS+' msout='+MS.nameMS+'.MS msin.datacolumn=DATA avg.freqstep=%i avg.timestep=%i' % (nchan, avg_time), \
            log=MS.nameMS+'_avg.log', commandType='DPPP')
s.run(check=True, maxThreads=20) # limit threads to prevent I/O isssues

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

# TEST
#logger.info("Put data to Jy...")
#MSs.run('taql "update $pathMS set DATA = 1e2*DATA"', log='$nameMS_taql.log', commandType='general')

#######################################################   
# Create columns (non compressed)
logger.info('Creating MODEL_DATA_LOWRES and SUBTRACTED_DATA...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA_HIGHRES', log='$nameMS_addcol.log', commandType='python')

########################################################   
# flag bad stations, and low-elev
logger.info('Flagging...')
MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS msout=. ant.baseline=\"'+bl2flag+'\"', \
            log='$nameMS_flag.log', commandType='DPPP')

# predict to save time MODEL_DATA
if hba: model_dir = '/home/fdg/scripts/model/AteamHBA/'+patch
else: model_dir = '/home/fdg/scripts/model/AteamLBA/'+patch

if os.path.exists(model_dir+'/img-MFS-model.fits'):
    logger.info('Predict (wsclean)...')
    s.add('wsclean -predict -name '+model_dir+'/img -j '+str(s.max_processors)+' -channelsout 15 '+MSs.getStrWsclean(), \
          log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    s.run(check=True)
else:
    logger.info('Predict (DPPP)...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel+' pre.sources='+patch, log='$nameMS_pre.log', commandType='DPPP')

for c in xrange(100):

    ####################################################
    # 1: find PA and remove it

    # Solve cal_SB.MS:DATA (only solve)
    logger.info('Calibrating PA...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal \
            sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solPA.log', commandType="DPPP")

    lib_util.run_losoto(s, 'pa-c'+str(c), [ms+'/pa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

    #################################################
    # 1: find the FR and remve it
    
    # Beam correction DATA -> CORRECTED_DATA
    logger.info('PA correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c'+str(c)+'.h5 cor.correction=polalign', \
            log='$nameMS_corPA2.log', commandType="DPPP")

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType='DPPP')
    
    # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)
    
    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Calibrating FR...')
    MSs.run('DPPP '+parset_dir+'/DPPP-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.parmdb=$pathMS/fr.h5 sol.caltype=diagonal', log='$nameMS_fr.log', commandType='DPPP')
    
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
    
    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Calibrating BP...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.parmdb=$pathMS/amp.h5 sol.caltype=diagonal', log='$nameMS_amp.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'amp-c'+str(c), [ms+'/amp.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-flag.parset',parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-bp.parset'])

    #################################################
    # 3: apply all

    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c'+str(c)+'.h5 cor.correction=polalign', \
            log='$nameMS_corPA4.log', commandType="DPPP")

    # Correct BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('BP correction...')
    if c == 0 and not hba:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.updateweights=True cor.parmdb=cal-amp-c'+str(c)+'.h5 cor.correction=amplitudeSmooth000', \
                log='$nameMS_corAMP4.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-amp-c'+str(c)+'.h5 cor.correction=amplitudeSmooth000', \
                log='$nameMS_corAMP4.log', commandType='DPPP')
 
    # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    if c == 0 and not hba:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam4.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_beam4.log', commandType='DPPP')
       
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR4.log', commandType='DPPP')
    
    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Calibrating IONO...')
    MSs.run('DPPP '+parset_dir+'/DPPP-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/iono.h5 sol.mode=diagonal', log='$nameMS_iono.log', commandType='DPPP')
    
    lib_util.run_losoto(s, 'iono-c'+str(c), [ms+'/iono.h5' for ms in MSs.getListStr()], \
        [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-plot-ph.parset'])
    
    # Correct all CORRECTED_DATA (PA, beam, BP, FR corrected) -> CORRECTED_DATA
    logger.info('IONO correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.steps=[ph,amp] cor.ph.parmdb=cal-iono-c'+str(c)+'.h5 cor.amp.parmdb=cal-iono-c'+str(c)+'.h5 \
                    cor.ph.correction=phase000 cor.amp.correction=amplitude000 cor.amp.updateweights=False', log='$nameMS_corG4.log', commandType="DPPP")

    # briggs: -1.2 for virgo; -1.0 for subtraction to get good minihalo?
    logger.info('Cleaning (cycle %i)...' % c)
    imagename = 'img/img-c'+str(c)
    if patch == 'CygA':
        s.add('wsclean -reorder -temp-dir '+ temp_dir +' -name ' + imagename + ' -size 1000 1000 -j '+str(s.max_processors)+' \
            -scale 1arcsec -weight uniform -niter 50000 -update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -multiscale -multiscale-scales 0,4,8,16,32 \
            -auto-threshold 1 \
            -use-idg \
            -join-channels -fit-spectral-pol 3 -channels-out 15 '+MSs.getStrWsclean(), \
            log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors = 'max')
        s.run(check = True)
    elif patch == 'VirA' and hba:
        s.add('wsclean -reorder -temp-dir '+ temp_dir +' -name ' + imagename + ' -size 2500 2500 -j '+str(s.max_processors)+' \
            -scale 1arcsec -weight briggs -1. -niter 1000 -update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -use-idg \
            -join-channels -fit-spectral-pol 3 -channels-out 15 '+MSs.getStrWsclean(), \
            log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors = 'max')
        s.run(check = True)
        s.add('wsclean -continue -reorder -temp-dir '+ temp_dir +' -name ' + imagename + ' -size 2500 2500 -j '+str(s.max_processors)+' \
            -scale 1arcsec -weight briggs -1. -niter 50000 -update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -multiscale -multiscale-scales 0,4,8,16,32,64 \
            -auto-threshold 1 \
            -use-idg \
            -join-channels -fit-spectral-pol 3 -channels-out 15 '+MSs.getStrWsclean(), \
            log='wscleanB-c'+str(c)+'.log', commandType='wsclean', processors = 'max')
        s.run(check = True)
    else:
        s.add('wsclean -reorder -temp-dir '+ temp_dir +' -name ' + imagename + ' -size 1500 1500 -j '+str(s.max_processors)+' \
            -scale 2arcsec -weight briggs -1.2 -niter 50000 -update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -multiscale -multiscale-scales 0,4,8,16,32 \
            -auto-threshold 1 \
            -use-idg \
            -join-channels -fit-spectral-pol 3 -channels-out 15 '+MSs.getStrWsclean(), \
            log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors = 'max')
        s.run(check = True)

    logger.info('Sub model...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql1.log', commandType='general')
    #logger.info('Copy MODEL_DATA...')
    #MSs.run('taql "update $pathMS set MODEL_DATA_HIGHRES = MODEL_DATA"', log='$nameMS_taql2.log', commandType='general')

    logger.info('Cleaning sub (cycle %i)...' % c)
    imagename = 'img/imgsub-c'+str(c)
    s.add('wsclean -reorder -temp-dir '+ temp_dir +' -name ' + imagename + ' -size 1000 1000 -j '+str(s.max_processors)+' \
            -scale 15arcsec -weight briggs -1.0 -taper-gaussian 80arcsec -niter 10000 -no-update-model-required -minuv-l 30 -mgain 0.85 -clean-border 1 \
            -multiscale -multiscale-scales 0,4,8,16 \
            -use-idg \
            -join-channels -fit-spectral-pol 2 -channels-out 15 '+MSs.getStrWsclean(), \
            log='wscleanSUB-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

    #logger.info('Combining MODEL_DATA_HIGHRES and MODEL_DATA...')
    #MSs.run('taql "update $pathMS set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"', log='$nameMS_taql3.log', commandType='general')

    imagename = 'img/img-c'+str(c)
    im = lib_img.Image(imagename)
    im.rescaleModel(f)
    logger.info('Predict (wsclean)...')
    s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channelsout 15 '+MSs.getStrWsclean(), \
          log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check = True)

logger.info("Done.")
