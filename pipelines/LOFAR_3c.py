#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation

import sys, os, glob, re
import numpy as np

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-3c.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_3c','parset_dir')
bl2flag = parset.get('flag','stations')
target = os.getcwd().split('/')[-1]
data_dir = '/home/fdg/lofar5/3Csurvey/%s' % target

def get_cal_dir(timestamp):
    """
    Get the proper cal directory from a timestamp
    """
    for cal_dir in glob.glob('../cals/*/*'):
        cal_timestamps = set( [ms.split('_')[1][1:] for ms in glob.glob(cal_dir+'/cals/*MS')] )
        if timestamp in cal_timestamps:
            logger.info('Calibrator found: %s (t=%s)' % (cal_dir, timestamp))
            return cal_dir

    logger.error('Missing calibrator.')
    sys.exit()

##########################################################
logger.info('Cleaning...')
lib_util.check_rm('cal*h5')
lib_util.check_rm('plots*')
lib_util.check_rm('img')
os.makedirs('img')

logger.info('Copy data...')
MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s, check_flags=False)
for timestamp in set([ os.path.basename(ms).split('_')[1][1:] for ms in MSs.getListStr() ]):
    mss_toconcat = glob.glob(data_dir+'/'+target+'_t'+timestamp+'_SB*.MS')
    MS_concat = target+'_t'+timestamp+'_concat.MS'
    MS_concat_bkp = target+'_t'+timestamp+'_concat.MS-bkp'
    if os.path.exists(MS_concat_bkp): 
        os.system('rm -r %s' % MS_concat)
        os.system('cp -r %s %s' % (MS_concat_bkp, MS_concat) )

    else:
        logger.info('Making %s...' % MS_concat)
        s.add('DPPP '+parset_dir+'/DPPP-avg.parset msin=\"'+str(mss_toconcat)+'\" msout='+MS_concat,\
            log=MS_concat+'_avg.log', commandType='DPPP')
	s.run(check=True, maxThreads=1)

        MSs = lib_ms.AllMSs( [MS_concat], s )
        
        # flag bad stations, and low-elev
        logger.info('Flagging...')
        MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS msout=. ant.baseline=\"'+bl2flag+'\"', \
                    log='$nameMS_flag.log', commandType='DPPP')
        
        cal_dir = get_cal_dir(timestamp)
        h5_pa = cal_dir+'/cal-pa.h5'
        h5_amp = cal_dir+'/cal-amp.h5'
        h5_iono = cal_dir+'/cal-iono.h5'
        assert os.path.exists(h5_pa)
        assert os.path.exists(h5_amp)
        assert os.path.exists(h5_iono)
        
        # Correct fist for BP(diag)+TEC+Clock and then for beam
        
        # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
        logger.info('Apply solutions (pa)...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.steps=[pa] \
                cor.pa.parmdb='+h5_pa+' cor.pa.correction=polalign', log='$nameMS_cor1.log', commandType='DPPP')
        
        # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
        logger.info('Apply solutions (amp/ph)...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
                cor.amp.parmdb='+h5_amp+' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
                cor.ph.parmdb='+h5_iono+' cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log', commandType='DPPP')
        
        # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
        logger.info('Beam correction...')
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DPPP')

        # Move CORRECTED_DATA -> DATA
        logger.info('Move CORRECTED_DATA -> DATA...')
        MSs.run('taql "update $pathMS set DATA = CORRECTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        # bkp
        logger.info('Making backup...')
        os.system('cp -r %s %s' % (MS_concat, MS_concat_bkp) ) # do not use MS.move here as it resets the MS path to the moved one

MSs = lib_ms.AllMSs( glob.glob('*concat.MS'), s, check_flags=False )

#####################################################
# Model
logger.info('Preparing model...')
if sourcedb is None:
    sourcedb = 'tgts.skydb'
    if not os.path.exists(sourcedb):
	phasecentre = MSs.getListObj()[0].getPhaseCentre()
        fwhm = MSs.getListObj()[0].getFWHM(freq='min')
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<1')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')

logger.info('Predict (DPPP)...')
MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+sourcedb, log='$nameMS_pre.log', commandType='DPPP')

# Smooth DATA -> SMOOTHED_DATA
logger.info('BL-based smoothing...')
MSs.run('BLsmooth.py -r -f 0.2 -i '+incol+' -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1-c'+str(c)+'.log', commandType='python')


###############################################################
# Selfcal
for c in range(100):

    logger.info('== Start cycle: %s ==' % c)

    logger.info('Remove bad timestamps...')
    MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    ####################################################
    # 1: solving

    # solve G - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving diag...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS sol.h5parm=$pathMS/cal.h5 sol.mode=rotation+diagonal', \
            log='$nameMS_solG-c'+str(c)+'.log', commandType="DPPP")

    lib_util.run_losoto(s, 'cal-c'+str(c), [ms+'/cal.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset'])

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving TEC...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS ddecal.h5parm=$pathMS/cal.h5', \
            log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DPPP')

    lib_util.run_losoto(s, 'cal-c'+str(c), [ms+'/cal.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-tec.parset'])

    sys.exit()


    #################################################
    # 2: find the FR and remve it
    

    
    # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)
    
    # Solve cal_SB.MS:CORRECTED_DATA (only solve)
    logger.info('Solving FR...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/fr.h5 sol.mode=diagonal \
            sol.uvlambdarange='+str(nouseblrange), log='$nameMS_solFR.log', commandType="DPPP")
    
    lib_util.run_losoto(s, 'fr-c'+str(c), [ms+'/fr.h5' for ms in MSs.getListStr()], \
            [parset_dir + '/losoto-fr.parset'])

    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
            log='$nameMS_corFR3.log', commandType="DPPP")


    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/img-'+str(c)
    lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=2000, scale='5arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
            multiscale='', multiscale_scales='0,10,20', \
            baseline_averaging=5, auto_threshold=1, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    logger.info('Predict (wsclean: %s)...' % imagename)
    s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_processors)+' -channels-out 61 '+MSs.getStrWsclean(), \
          log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

logger.info("Done.")
