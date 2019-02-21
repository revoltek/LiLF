#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt

if 'Vir' in os.getcwd():
    patch = 'VirA'
    nouseblrange = ''
    f = lambda nu: 1226. * 10**(-0.79 * (np.log10(nu/150.e6))**1)
elif 'Tau' in os.getcwd():
    patch = 'TauA'
    nouseblrange = '' #'[500..5000]' # below is a point, above 10 times is hopefully resolved out
    f = lambda nu: 1838. * 10**(-0.299 * (np.log10(nu/150.e6))**1)
elif 'Cas' in os.getcwd():
    patch = 'CasA'
    nouseblrange = '' #'[15000..1e30]'
    f = lambda nu: 11733. * 10**(-0.77 * (np.log10(nu/150.e6))**1)
elif 'Cyg' in os.getcwd():
    patch = 'CygA'
    nouseblrange = ''
    nouseblrange = '[15000..1e30]'
    f = lambda nu: 10690. * 10**(-0.67 * (np.log10(nu/150.e6))**1) * 10**(-0.204 * (np.log10(nu/150.e6))**2) * 10**(-0.021 * (np.log10(nu/150.e6))**3)

skymodel = '/home/fdg/scripts/model/A-team_4_CC.skydb'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-ateam.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_ateam','parset_dir')
bl2flag = parset.get('flag','stations')
data_dir = '../tgts-bkp/'

##########################################################
logger.info('Cleaning...')
lib_util.check_rm('cal*h5')
lib_util.check_rm('plots*')
lib_util.check_rm('img')
os.makedirs('img')
MSs = lib_ms.AllMSs( sorted(glob.glob(data_dir+'/*MS')), s )

# copy data (avg to 1ch/sb and 10 sec)
nchan = MSs.getListObj()[0].getNchan()
timeint = MSs.getListObj()[0].getTimeInt()
avg_time = int(np.rint(10./timeint))

logger.info('Copy data...')
for obs in set([ os.path.basename(ms).split('_')[0] for ms in MSs.getListStr() ]):
    mss_toconcat = glob.glob(data_dir+'/'+obs+'*MS')
    MS_concat = obs+'_concat.MS'
    MS_concat_bkp = obs+'_concat.MS-bkp'
    if os.path.exists(MS_concat_bkp): 
        os.system('rm -r %s' % MS_concat)
        os.system('cp -r %s %s' % (MS_concat_bkp, MS_concat) )
    else:
        s.add('DPPP '+parset_dir+'/DPPP-avg.parset msin=\"'+str(mss_toconcat)+'\" msout='+MS_concat+' avg.freqstep=%i avg.timestep=%i' % (nchan, avg_time),\
            log=obs+'_avg.log', commandType='DPPP')
s.run(check=True, maxThreads=2)

################################################################
MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

# bkp
for MS in MSs.getListStr():
    MS_bkp = MS+'-bkp'
    if not os.path.exists(MS_bkp):
        logger.info('Making backup...')
        os.system('cp -r %s %s' % (MS, MS_bkp) ) # do not use MS.move here as it resets the MS path to the moved one

# HBA/LBA
if min(MSs.getFreqs()) < 80.e6:
    lofar_system = 'lba'
    flag_steps = "[ant, uvmin, elev, count]"
else: 
    lofar_system = 'hba'
    flag_steps = "[ears, ant, uvmin, elev, count]"

########################################################   
# flag bad stations, and low-elev
logger.info('Flagging...')
MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS msout=. steps=\"'+flag_steps+'\" ant.baseline=\"'+bl2flag+'\"', \
            log='$nameMS_flag.log', commandType='DPPP')

# predict to save time MODEL_DATA
if lofar_system == 'hba': model_dir = '/home/fdg/scripts/model/AteamHBA/'+patch
else: model_dir = '/home/fdg/scripts/model/AteamLBA/'+patch

if os.path.exists(model_dir+'/img-MFS-model.fits'):
    im = lib_img.Image(model_dir+'/img')
    im.rescaleModel(f)
    n = len(glob.glob(model_dir+'/img-[0-9]*-model.fits'))
    logger.info('Predict (wsclean: %s - chan: %i)...' % (model_dir, n))
    s.add('wsclean -predict -name '+model_dir+'/img -j '+str(s.max_processors)+' -channels-out '+str(n)+' '+MSs.getStrWsclean(), \
          log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    s.run(check=True)
else:
    logger.info('Predict (DPPP)...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel+' pre.sources='+patch, log='$nameMS_pre.log', commandType='DPPP')

for c in xrange(100):

    logger.info('== Start cycle: %s ==' % c)

    #logger.info('Remove bad timestamps...')
    #MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    ####################################################
    # 1: find PA and remove it

    # Solve cal_SB.MS:DATA (only solve)
    logger.info('Solving PA...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal \
            sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solPA.log', commandType="DPPP")

    lib_util.run_losoto(s, 'pa-c'+str(c), [ms+'/pa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

    #################################################
    # 2: find the FR and remve it
    
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
    logger.info('Solving FR...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/fr.h5 sol.mode=diagonal \
            sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solFR.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'fr-c'+str(c), [ms+'/fr.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-fr.parset'])

    #####################################################
    # 3: find BANDPASS/IONO

    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Polalign correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa-c'+str(c)+'.h5 cor.correction=polalign', \
            log='$nameMS_corPA3.log', commandType="DPPP")

    # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    if c == 0:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_corBEAM3.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_corBEAM3.log', commandType='DPPP')
 
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR3.log', commandType="DPPP")

    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Solving IONO...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/iono.h5 sol.mode=diagonal \
            sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solIONO3.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'iono-c'+str(c), [ms+'/iono.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-plot-ph.parset'])
    
    # Correct all CORRECTED_DATA -> CORRECTED_DATA
    logger.info('IONO correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-iono-c'+str(c)+'.h5 cor.correction=phase000', \
                log='$nameMS_corIONO3.log', commandType='DPPP')

    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Solving BP...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/amp.h5 sol.mode=diagonal \
            sol.uvlambdarange='+str(nouseblrange)+' sol.nchan=2 sol.solint=10', log='$nameMS_solAMP3.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'amp-c'+str(c), [ms+'/amp.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-plot-amp.parset'])

    # Correct BP CORRECTED_DATA -> CORRECTED_DATA
    logger.info('BP correction...')
    if c == 0:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.updateweights=True cor.parmdb=cal-amp-c'+str(c)+'.h5 cor.correction=amplitude000', \
                log='$nameMS_corAMP3.log', commandType='DPPP')
    else:
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.updateweights=False cor.parmdb=cal-amp-c'+str(c)+'.h5 cor.correction=amplitude000', \
                log='$nameMS_corAMP3.log', commandType='DPPP')
       
    # briggs: -1.2 for virgo; -1.0 for subtraction to get good minihalo?
    logger.info('Cleaning (cycle %i)...' % c)
    imagename = 'img/img-c'+str(c)
    if patch == 'CygA':
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1000, scale='1.5arcsec', \
                weight='briggs -1', niter=50000, no_update_model_required='', mgain=0.5, \
                iuwt='', gain=0.2, \
                #multiscale='', \
                #multiscale_scale_bias=0.6, \
                #multiscale_scales='0,10,20,40', \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/CygA.fits', \
                baseline_averaging=5, deconvolution_channels=12, \
                auto_threshold=1, join_channels='', fit_spectral_pol=5, channels_out=61)

    elif patch == 'CasA':
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1300, scale='2arcsec', \
                weight='briggs -1', niter=75000, no_update_model_required='', mgain=0.5, \
                multiscale='', multiscale_scale_bias=0.7, \
                # multiscale_scales='0,5,10,20,40,80', \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/CasA.fits', \
                baseline_averaging=5, deconvolution_channels=12, \
                auto_threshold=1, join_channels='', fit_spectral_pol=5, channels_out=61)

    elif patch == 'TauA':
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1200, scale='2arcsec', \
                weight='briggs -1', niter=75000, no_update_model_required='', mgain=0.5, \
                multiscale='', multiscale_scale_bias=0.7, \
                multiscale_scales='0,5,10,20,40,80', \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/TauA.fits', \
                baseline_averaging=5, deconvolution_channels=12, \
                auto_threshold=1, join_channels='', fit_spectral_pol=5, channels_out=61)

    elif patch == 'VirA' and lofar_system == 'lba':
        #lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1500, scale='2arcsec', \
        #        weight='briggs -1.', niter=1000, update_model_required='', mgain=0.85, \
        #        deconvolution_channels=12, \
        #        join_channels='', fit_spectral_pol=2, channels_out=61) # use cont=True
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1500, scale='2arcsec', \
                weight='briggs -1.', niter=50000, no_update_model_required='', mgain=0.5, \
                multiscale='', multiscale_scale_bias=0.7, \
                # multiscale_scales='0,5,10,20,40,80', \
                #casa_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirA.crtf', \
                baseline_averaging=5, deconvolution_channels=12, \
                auto_threshold=1, join_channels='', fit_spectral_pol=5, channels_out=61)

    elif patch == 'VirA' and lofar_system == 'hba':
        lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=2500, scale='1arcsec', \
                weight='briggs -1.', niter=1000, update_model_required='', mgain=0.5, \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirAphba.fits', \
                deconvolution_channels=12, \
                join_channels='', fit_spectral_pol=5, channels_out=61) # use cont=True
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), cont=True, name=imagename, size=2500, scale='1arcsec', \
                weight='briggs -1.', niter=50000, no_update_model_required='', mgain=0.5, \
                multiscale='', multiscale_scale_bias=0.7, \
                # multiscale_scales='0,5,10,20,40,80', \
                fits_mask='/home/fdg/scripts/LiLF/parsets/LOFAR_ateam/masks/VirAhba.fits', \
                baseline_averaging=5, deconvolution_channels=12, \
                auto_threshold=0.1, join_channels='', fit_spectral_pol=5, channels_out=61)

    logger.info('Predict (wsclean: %s)...' % imagename)
    s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channels-out 61 '+MSs.getStrWsclean(), \
          log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

    # every 5 cycles: sub model and rescale model
    # TODO: add low-res model to global model
    if c%5 == 0 and c != 0:
    
        logger.info('Sub model...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql1.log', commandType='general')

        #logger.info('Reweight...')
        #MSs.run('reweight.py $pathMS -v -m residual -d CORRECTED_DATA', log='$nameMS_weights.log', commandType='python')

        logger.info('Cleaning sub (cycle %i)...' % c)
        imagename = 'img/imgsub-c'+str(c)
        lib_util.run_wsclean(s, 'wscleanSUB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=1000, scale='15arcsec', \
                weight='briggs 0.', taper_gaussian='120arcsec', niter=10000, no_update_model_required='', baseline_averaging=5, minuv_l=30, mgain=0.85, \
                multiscale='', multiscale_scales='0,4,8,16', \
                auto_threshold=1, join_channels='', fit_spectral_pol=3, channels_out=10)
 
        #imagename = 'img/img-c'+str(c)
        #im = lib_img.Image(imagename)
        #im.rescaleModel(f)
        #logger.info('Predict (wsclean)...')
        #s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channelsout 15 '+MSs.getStrWsclean(), \
        #      log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
        #s.run(check = True)

logger.info("Done.")
