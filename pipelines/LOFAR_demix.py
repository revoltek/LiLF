#!/usr/bin/env python
# -*- coding: utf-8 -*-
# demix of a set of SBs from a given dir, output is in the local dir

import sys, os, glob
import lsmtool

###############################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-demix.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_demix','parset_dir')
data_dir = parset.get('LOFAR_demix','data_dir')
demix_skymodel = parset.get('LOFAR_demix','demix_model')  # skymodel to be converted in skydb
include_target = parset.getboolean('LOFAR_demix','include_target')

##############################################
MSs = lib_ms.AllMSs(glob.glob('mss-predemix/TC*[0-9].MS'), s)
ateams = ['VirA', 'TauA']
ateams_todemix = []
for ateam in ateams:
    sep = MSs.getListObj()[0].distBrightSource(ateam)
    if sep < 4 or sep > 15:
        logger.debug('No demix of %s (sep: %.0f deg)' % (ateam, sep))
    else:
        ateams_todemix.append(ateam)
        logger.warning('Demix of %s (sep: %.0f deg)' % (ateam, sep))

if len(ateams_todemix) < 0:
    logger.info('Notingh to demix.')
    sys.exit()

if os.path.exists('mss-predemix'):
    logger.warning('Reset mss...')
    lib_util.check_rm('mss/*MS')
else:
    logger.info('Move mss in mss-predemix...')
    os.system('mv mss mss-predemix')
    os.system('mkdir mss')

MSs = lib_ms.AllMSs(glob.glob('mss-predemix/TC*[0-9].MS'), s)

if include_target:
    if not os.path.exists('demix_combined.skymodel'):
        logger.info('Including target...')
        fwhm = MSs.getListObj()[0].getFWHM(freq='min')
        phasecentre = MSs.getListObj()[0].getPhaseCentre()
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O demix_tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
        lsm = lsmtool.load('demix_tgts.skymodel')
        lsm.remove('I<1')
        lsm.group('single', root='target')
        lsm.setColValues('LogarithmicSI', ['true']*len(lsm))
        # join with ateam skymodel
        lsm.concatenate(demix_skymodel)
        lsm.write('demix_combined.skymodel', clobber=True)
    demix_skymodel = 'demix_combined.skymodel'
    targetsource = 'target'
else:
    targetsource = ''

logger.info('Creating skydb...')
lib_util.check_rm('demix.skydb')
os.system('makesourcedb outtype="blob" format="<" in=%s out=demix.skydb' % demix_skymodel)

for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/demix.skydb')
    os.system('cp -r demix.skydb '+MS+'/demix.skydb')

logger.info('Demixing...')
MSs.run('DP3 '+parset_dir+'/DP3-demix.parset msin=$pathMS msout=mss/$nameMS.MS demixer.skymodel=$pathMS/demix.skydb \
        demixer.instrumentmodel=$pathMS/instrument_demix demixer.targetsource='+targetsource+'\
        demixer.subtractsources=\['+','.join(ateams_todemix)+'\]',
        log='$nameMS_demix.log', commandType='DP3', maxThreads=1)

logger.info("Done.")
