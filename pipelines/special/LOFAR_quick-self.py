#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# they need to be in "./mss/"

import os, glob
import lsmtool

########################################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-quick-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-quick-self.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_quick-self','parset_dir')
sourcedb = parset.get('model','sourcedb')
userReg = parset.get('model','userReg')
data_dir = parset.get('LOFAR_quick-self','data_dir')

#############################################################################
# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')
    lib_util.check_rm('cal-quick.h5')

### DONE
    
################################################################
# Copy dataset
with w.if_todo('copy'):
    MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )
    s.add('DP3 ' + parset_dir + '/DP3-concat.parset msin="['+MSs.getStrDDF()+']" msout=quick.MS',
            log='concat.log', commandType='DP3')
    s.run(check=True)
    
### DONE

MSs = lib_ms.AllMSs( ['quick.MS'], s )

# set image size
imgsizepix = int(2.1*MSs.getListObj()[0].getFWHM(freq='mid')*3600/60.)
if imgsizepix%2 != 0: imgsizepix += 1 # prevent odd img sizes

# set clean componet fit order (use 5 for large BW)
bandwidth = MSs.getBandwidth()
if bandwidth > 25e6: cc_fit_order = 5
else: cc_fit_order = 3

#################################################################
# Get online model
if sourcedb == '':
    if not os.path.exists('tgts.skydb'):
        fwhm = MSs.getListObj()[0].getFWHM(freq='min')
        phasecentre = MSs.getListObj()[0].getPhaseCentre()
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<1')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')
    sourcedb = 'tgts.skydb'

#################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS + '/' + sourcedb_basename)
    logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
    os.system('cp -r ' + sourcedb + ' ' + MS)

with w.if_todo('init_model'):
    # note: do not add MODEL_DATA or the beam is transported from DATA, while we want it without beam applied
    logger.info('Creating CORRECTED_DATA...')
    MSs.addcol('CORRECTED_DATA', 'DATA')

    logger.info('Add model to MODEL_DATA...')
    MSs.run('DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=true pre.sourcedb=$pathMS/' + sourcedb_basename,
            log='$nameMS_pre.log', commandType='DP3')
### DONE

with w.if_todo('beam'):
    logger.info('Set CORRECTED_DATA = DATA...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = DATA"', log='$nameMS_taql.log', commandType='general')

    # correct beam - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correct beam')
    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=False', log='$nameMS_beam.log', commandType='DP3')
    #MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DP3')
### DONE

with w.if_todo('calibrate'):
    # Smooth CORRECTED_DATA -> SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -c 8 -n 8 -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')
    MSs.run('BLsmooth.py -c 8 -n 8 -r -i MODEL_DATA -o MODEL_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')

    # solve amp+ph - ms:SMOOTHED_DATA
    logger.info('Solving amp+ph...')
    MSs.run('DP3 '+parset_dir+'/DP3-solG.parset msin=$pathMS sol.h5parm=$pathMS/quick.h5 \
           	sol.flagunconverged=True sol.flagdivergedonly=True sol.solint=30 sol.nchan=8', \
                log='$nameMS_solG.log', commandType='DP3')

    lib_util.run_losoto(s, 'quick', [ms+'/quick.h5' for ms in MSs.getListStr()], [parset_dir+'/losoto-plot.parset'])
### DONE

with w.if_todo('correct'):
    # correct amp+ph - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correcting amp+ph...')
    #MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA\
    #            cor.parmdb=cal-quick.h5 cor.correction=amplitude000 cor.missingantennabehavior=unit',
    #            log='$nameMS_cor.log', commandType='DP3')
    MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA\
                cor.parmdb=cal-quick.h5 cor.correction=phase000 cor.missingantennabehavior=unit',
                log='$nameMS_cor.log', commandType='DP3')

### DONE

imagename = 'img/quick'
with w.if_todo('imaging'):
    logger.info('Cleaning...')
    # make temp mask for cycle 0, in cycle 1 use the maske made from cycle 0 image
    lib_util.run_wsclean(s, 'wsclean.log', MSs.getStrWsclean(), name=imagename,
                                 size=imgsizepix, scale='60arcsec',
                                 weight='briggs -0.5', niter=100000, no_update_model_required='',
                                 parallel_gridding=2, baseline_averaging='', mgain=0.85,
                                 parallel_deconvolution=512, local_rms='', auto_threshold=4,
                                 join_channels='', fit_spectral_pol=cc_fit_order, channels_out=MSs.getChout(4.e6),
                                 deconvolution_channels=cc_fit_order)

logger.info("Done.")
