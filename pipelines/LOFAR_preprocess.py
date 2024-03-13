#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, re, glob, time

import lsmtool
import numpy as np
import casacore.tables as pt
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u

##########################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-preprocess')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-preprocess.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_preprocess'])))
parset_dir = parset.get('LOFAR_preprocess','parset_dir')
fix_table = parset.getboolean('LOFAR_preprocess','fix_table')
renameavg = parset.getboolean('LOFAR_preprocess','renameavg')
keep_IS = parset.getboolean('LOFAR_preprocess','keep_IS')
backup_full_res = parset.getboolean('LOFAR_preprocess','backup_full_res')
demix_sources = parset.get('LOFAR_preprocess','demix_sources') # Demix  sources in these patches (e.g. CasA or [VirA,TauA]), default: No demix
demix_skymodel = parset.get('LOFAR_preprocess','demix_skymodel') # Use non-default demix skymodel
demix_field_skymodel = parset.get('LOFAR_preprocess','demix_field_skymodel') # provide a custom target skymodel instead of online gsm model


###########################################
if os.path.exists('html.txt'):
    download_file = 'html.txt'
else:
    download_file = None # just renaming

def getName(ms):
    """
    Get new MS name based on obs name and time
    """
    # get pointing name
    with pt.table(ms+'/FIELD', readonly=True, ack=False) as t:
        code = t.getcell('CODE',0)
    if code == '':
        with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
            code = t.getcell('LOFAR_TARGET',0)[0]
    code = code.lower().replace(' ','_')
    
    # get obsid
    with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
        obsid = t.getcell('LOFAR_OBSERVATION_ID',0)

    # get freq
    with pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False) as t:
        freq = t.getcell('REF_FREQUENCY',0)

    # get time (saved in ms as MJD in seconds)
    #with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
    #    time = Time(t.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
    #    time = time.iso.replace('-','').replace(' ','').replace(':','')[0:12]

    #pattern = re.compile("^c[0-9][0-9]-.*$")
    # is survey?
    #if pattern.match(code):
    #    cycle_obs, sou = code.split('_')
    #    if not os.path.exists(cycle_obs+'/'+sou): os.makedirs(cycle_obs+'/'+sou)
    #    return cycle_obs+'/'+sou+'/'+sou+'_t'+time+'_SB'+str(lib_util.lofar_nu2num(freq/1.e6))+'.MS'
    #else:
    
    if not os.path.exists('mss/id'+obsid+'_-_'+code): os.makedirs('mss/id'+obsid+'_-_'+code)
    return 'mss/id'+obsid+'_-_'+code+'/'+code+'_SB%03i.MS' % lib_util.lofar_nu2num(freq/1.e6)

########################################
if not download_file is None:
    with w.if_todo('download'):
       with open(download_file,'r') as df:
            logger.info('Downloading...')
            downloaded = glob.glob('*MS')
            # add renamed files
            if os.path.exists('renamed.txt'):
                with open('renamed.txt','r') as flog:
                    downloaded += [line.rstrip('\n') for line in flog]

            for i, line in enumerate(df):
                ms = re.findall(r'L[0-9]*.*_SB[0-9]*_uv', line)[0]
                if ms+'.MS' in downloaded or ms+'.dppp.MS' in downloaded: continue
                if ms+'.MS' in glob.glob('*MS') or ms+'.dppp.MS' in glob.glob('*MS'): continue
                s.add('wget -nv "'+line[:-1]+'" -O - | tar -x', log=ms+'_download.log', commandType='general')
            #    print 'wget -nv "'+line[:-1]+'" -O - | tar -x'
                logger.debug('Queue download of: '+line[:-1])
            s.run(check=True, maxThreads=4)

MSs = lib_ms.AllMSs(glob.glob('*MS'), s, check_flags=False)
if len(MSs.getListStr()) == 0:
    logger.info('Done.')
    sys.exit(0)

######################################
with pt.table(MSs.getListStr()[0]+'/OBSERVATION', readonly=True, ack=False) as obs:
    t = Time(obs.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
    time = np.int(t.iso.replace('-','')[0:8])

if fix_table:
    with w.if_todo('fix_table'):
        #logger.info('Fix MS table...')
        #MSs.run('fixMS_TabRef.py $pathMS', log='$nameMS_fixms.log', commandType='python')
        # only ms created in range (2/2013->2/2014)
        if time > 20130200 and time < 20140300:
            logger.info('Fix beam table...')
            MSs.run('/home/fdg/scripts/fixinfo/fixbeaminfo $pathMS', log='$nameMS_fixbeam.log', commandType='python')

# Rescale visibilities by 1e3 if before 2014-03-19 (old correlator), and by 1e-2 otherwise
with w.if_todo('rescale_flux'):
    logger.info('Rescaling flux...')
    if time < 20140319:
        rescale_factor = 1e6
    else:
        rescale_factor = 1e-4

    for MS in MSs.getListStr():
        with pt.table(MS+'/HISTORY', readonly=False, ack=False) as hist:
            if "Flux rescaled" not in hist.getcol('MESSAGE'):
                s.add('taql "update %s set DATA = %f*DATA" && taql "insert into %s/HISTORY (TIME,MESSAGE) values (mjd(), \'Flux rescaled\')"' % (MS,rescale_factor,MS), \
                        log='taql.log', commandType='general')
    s.run(check=True)

######################################
# Avg to 4 chan and 2 sec
# Remove internationals
if renameavg:
    with w.if_todo('renameavg'):
        logger.info('Renaming/averaging...')
        with open('renamed.txt','a') as flog:
            MSs = lib_ms.AllMSs([MS for MS in glob.glob('*MS') if not os.path.exists(getName(MS))], s, check_flags=False)
            minfreq = np.min(MSs.getFreqs())
            fwhm = MSs.getListObj()[0].getFWHM(freq='min') # for radius of model
            logger.info('Min freq: %.2f MHz' % (minfreq/1e6))
            for MS in MSs.getListObj():

                if np.all(MS.getFreqs() > 168.3e6):
                    logger.warning(f'Skipping HBA above 168 MHz: deleting {MS.pathMS}')
                    lib_util.check_rm(MS.pathMS)
                    continue

                # get avg time/freq values
                nchan = MS.getNchan()
                timeint = MS.getTimeInt()

                if nchan == 1:
                    avg_factor_f = 1
                elif nchan % 2 == 0 and MSs.isHBA: # case HBA
                    avg_factor_f = int(nchan / 4)  # to 2 ch/SB
                elif nchan % 8 == 0 and minfreq < 40e6:
                    avg_factor_f = int(nchan / 8)  # to 8 ch/SB
                elif nchan % 8 == 0 and 'SPARSE' in MS.getAntennaSet():
                    avg_factor_f = int(nchan / 8)  # to 8 ch/SB
                elif nchan % 4 == 0:
                    avg_factor_f = int(nchan / 4)  # to 4 ch/SB
                elif nchan % 5 == 0:
                    avg_factor_f = int(nchan / 5)  # to 5 ch/SB
                else:
                    logger.error('Channels should be a multiple of 4 or 5.')
                    sys.exit(1)

                if keep_IS:
                     avg_factor_f = int(nchan / 16) if MSs.isHBA else int(nchan / 16) # to have the full FoV in LBA we need 32 ch/SB
                if avg_factor_f < 1: avg_factor_f = 1

                avg_factor_t = int(np.round(2/timeint)) if keep_IS else int(np.round(4/timeint)) # to 4 sec (2 for IS)
                if avg_factor_t < 1: avg_factor_t = 1

                # do not select IS baselines expect if specified
                bl_sel = '"*&"' if keep_IS else '"[CR]S*&"'

                MSout = getName(MS.pathMS)

                if avg_factor_f != 1 or avg_factor_t != 1 or demix_sources:
                    print(demix_sources)
                    # normal case: no demixing, just averaging
                    if not demix_sources:
                        logger.info('%s->%s: Average in freq (factor of %i) and time (factor of %i)...' % (MS.nameMS, MSout, avg_factor_f, avg_factor_t))

                        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin={MS.pathMS} msin.baseline={bl_sel} msout={MSout} '
                              f'msin.datacolumn=DATA avg.timestep={avg_factor_t} avg.freqstep={avg_factor_f}',
                              log=MS.nameMS+'_avg.log', commandType='DP3')
                        s.run(check=True, maxThreads=1) # limit threads to prevent I/O isssues
                    # special case: run demix and average in demixer call
                    else:
                        # check HBA - not tested
                        if MSs.isHBA:
                            logger.warning('Demixing in HBA not tested, please be careful!')
                            demix_sm_file = os.path.dirname(__file__) + '/../models/A-Team_lowres.skymodel' if not demix_skymodel else demix_skymodel
                        else:
                            demix_sm_file = os.path.dirname(__file__) + '/../models/demix_all.skymodel' if not demix_skymodel else demix_skymodel
                        # Check demix sources
                        demix_sm = lsmtool.load(demix_sm_file, beamMS=MS.pathMS)
                        demix_sm_patches = demix_sm.getPatchPositions()
                        ra, dec = MS.getPhaseCentre() # to calculate distance
                        demix_sm_sources = demix_sm_patches.keys()
                        if not isinstance(demix_sources, list):
                            demix_sources = demix_sources.replace('[','').replace(']','').split(',')
                        # check distances and that all sources in config are actually in skymodel.
                        for demix_source in demix_sources:
                            if demix_source not in demix_sm_sources:
                                logger.error(f'demix_source={demix_source} not in demix skymodel sources {demix_sm_sources}!')
                                sys.exit()
                            else:
                                coord_demix = SkyCoord(ra=demix_sm_patches[demix_source][0], dec=demix_sm_patches[demix_source][1])
                                sep = coord_demix.separation(SkyCoord(ra * u.deg, dec * u.deg)).deg
                                logger.info(f'Demix source: {demix_source} separation: {sep:.2f}deg')

                        # check target field skymodel
                        # if not provided, use GSM
                        if demix_field_skymodel == '':
                            logger.info('Include target from GSM...')
                            # get model the size of the image (radius=fwhm/2)
                            os.system('wget -O demix_tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (
                                    ra, dec, fwhm/2))  # ASTRON
                            lsm = lsmtool.load('demix_tgts.skymodel', beamMS=MS.pathMS)
                            lsm.remove('I<1')
                        else:
                            lsm = lsmtool.load(demix_field_skymodel, beamMS=MS.pathMS)
                        lsm.group('single', root='target')
                        lsm.setColValues('LogarithmicSI', ['true'] * len(lsm))
                        # join with ateam skymodel
                        lsm.concatenate(demix_sm)
                        # let's try to apply the beam, will probably be more accurate?
                        # TODO apply beam probably better, but creates bug?
                        lsm.write(f'{MS.pathMS}/demix_combined.skymodel', clobber=True, applyBeam=False)
                        # lsm.write(f'{MS.pathMS}/demix_combined.skymodel', clobber=True, applyBeam=True)
                        logger.info('%s->%s: Demix %s and average in freq (factor of %i) and time (factor of %i)...' % (MS.nameMS, MSout, demix_sources ,avg_factor_f, avg_factor_t))
                        demix_f_step = nchan
                        demix_t_step = int(np.round(8/timeint))

                        s.add(f'DP3 {parset_dir}/DP3-demix.parset msin={MS.pathMS} msin.baseline={bl_sel} msout={MSout} '
                              f'demix.skymodel={MS.pathMS}/demix_combined.skymodel demix.instrumentmodel={MS.pathMS}/instrument_demix.parmdb '
                              f'demix.targetsource=target demix.subtractsources=\[{",".join(demix_sources)}\] '
                              f'demix.freqstep={avg_factor_f} demix.timestep={avg_factor_t} '
                              f'demix.demixfreqstep={demix_f_step} demix.demixtimestep={demix_t_step}',
                              log=MS.nameMS+'_demix.log', commandType='DP3')
                        s.run(check=True, maxThreads=None) # We do not want to limit threads since this is CPU limited?
                    if backup_full_res:
                        logger.info('Backup full resolution data...')
                        if not os.path.exists('data-bkp'):
                            os.makedirs('data-bkp')
                        for MS in MSs.getListObj():
                            MS.move('data-bkp/' + MS.nameMS + '.MS', keepOrig=False, overwrite=False)
                    else:
                        lib_util.check_rm(MS.pathMS)
                    flog.write(MS.nameMS+'.MS\n') # after averaging to be sure no log is written if an error occurs
                else:
                    logger.info('%s->%s: Move data - no averaging...' % (MS.nameMS, MSout))
                    flog.write(MS.nameMS+'.MS\n') # before move or the filenmae is changed
                    MS.move(MSout)


logger.info("Done.")
