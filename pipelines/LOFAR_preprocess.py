#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, re, glob

import lsmtool
import numpy as np
import casacore.tables as pt
from astropy.coordinates import SkyCoord
import astropy.units as u

##########################################
from LiLF import lib_ms, lib_util, lib_log, lib_dd_parallel, lib_cat
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
keep_IS = parset.getboolean('LOFAR_preprocess','keep_IS') # default True
backup_full_res = parset.getboolean('LOFAR_preprocess','backup_full_res')
demix_sources = parset.get('LOFAR_preprocess','demix_sources') # demix the sources in these patches (e.g. CasA or [VirA,TauA]), default: No demix. Assumes intrinsic sky
demix_skymodel = parset.get('LOFAR_preprocess','demix_skymodel') # Use non-default demix skymodel
demix_field_skymodel = parset.get('LOFAR_preprocess','demix_field_skymodel') # provide a custom target skymodel instead of online gsm model - assumes intrinsic sky.
run_aoflagger = parset.getboolean('LOFAR_preprocess','run_aoflagger') # run aoflagger on individual subbands - do this only in rare cases where it was not done by the observatory!
tar = parset.getboolean('LOFAR_preprocess','tar') # tar the output ms
data_dir = parset.get('LOFAR_preprocess','data_dir') # directory where the data is stored

###########################################
if os.path.exists('html.txt'):
    download_file = 'html.txt'
else:
    download_file = None # just renaming

def getName(ms):
    """
    Get new MS name based on target name and obs id
    """
    # get pointing name
    with pt.table(ms+'/FIELD', readonly=True, ack=False) as t:
        code = t.getcell('CODE',0)
    if code == '':
        with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
            code = t.getcell('LOFAR_TARGET',0)[0]
    code = code.lower().replace(' ','_')
    
    # remove unecessary info for survey pointings
    if len(re.findall(r'[P|p]\d{3}\+\d{2}',code)) != 0:
        code =  re.findall(r'[P|p]\d{3}\+\d{2}',code)[0].upper()
    
    # get obsid
    with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
        obsid = t.getcell('LOFAR_OBSERVATION_ID',0)

    # get freq
    with pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False) as t:
        freq = t.getcell('REF_FREQUENCY',0)
    
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
                s.add('wget -nv --no-check-certificate "'+line[:-1]+'" -O - | tar -x', log=ms+'_download.log', commandType='general')
            #    print 'wget -nv "'+line[:-1]+'" -O - | tar -x'
                logger.debug('Queue download of: '+line[:-1])
            s.run(check=True, maxProcs=4)

MSs = lib_ms.AllMSs(glob.glob(data_dir + '*MS'), s, check_flags=False)
if len(MSs.getListStr()) == 0:
    logger.info('Done.')
    sys.exit(0)

if keep_IS and not MSs.hasIS:
    logger.debug('Keeping IS requested but IS were not recorded - switching to averaging settings for Dutch only.')
    keep_IS = False

######################################
if len(MSs.getListObj()) > 10000:
    logger.warning('Many MSs detected, using only the first to determine the observing time (for rescaling/fixtables).')
    t = MSs.getListObj()[0].get_time()
    times = [int(t.iso.replace('-','')[0:8])] * len(MSs.getListObj())
else:
    times = []
    for MS in MSs.getListObj():
        t = MS.get_time()
        times.append(int(t.iso.replace('-','')[0:8]))
    
if run_aoflagger:
    with w.if_todo('flag'):
        # Flag in an identical way to the observatory flagging
        logger.info('Flagging...')
        MSs.run('DP3 ' + parset_dir + '/DP3-flag.parset msin=$pathMS aoflagger.strategy=' + parset_dir + '/LBAdefaultwideband.lua',
            log='$nameMS_flag.log', commandType='DP3') # there might be a better way of parallelizing as this might be I/O or memory limited

if fix_table:
    with w.if_todo('fix_table'):
        #logger.info('Fix MS table...')
        #MSs.run('fixMS_TabRef.py $pathMS', log='$nameMS_fixms.log', commandType='python')
        # only ms created in range (2/2013->2/2014)
        for i, MS in enumerate(MSs.getListStr()):
            if times[i] > 20130200 and times[i] < 20140300:
                logger.info('Fix beam table...')
                s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+MS, log='fixbeam.log', commandType='python')
        s.run(check=True)

# Rescale visibilities by 1e3 if before 2014-03-19 (old correlator), and by 1e-2 otherwise
with w.if_todo('rescale_flux'):
    logger.info('Rescaling flux...')

    for i, MS in enumerate(MSs.getListStr()):
        if times[i] < 20140319:
            rescale_factor = 1e6
        else:
            rescale_factor = 1e-4

        with pt.table(MS+'/HISTORY', readonly=False, ack=False) as hist:
            if "Flux rescaled" not in hist.getcol('MESSAGE'):
                s.add('taql "update %s set DATA = %f*DATA" && taql "insert into %s/HISTORY (TIME,MESSAGE) values (mjd(), \'Flux rescaled\')"' % (MS,rescale_factor,MS), \
                        log='taql.log', commandType='general')
    s.run(check=True)

######################################
# Averaging/demixing/removing IS (if requested)
if renameavg:
    with w.if_todo('renameavg'):
        logger.info('Renaming/averaging...')
        with open('renamed.txt','a') as flog:
            MSs = lib_ms.AllMSs([MS for MS in glob.glob(data_dir + '*MS') if not os.path.exists(getName(MS))], s, check_flags=False)
            minfreq = np.min(MSs.getFreqs())
            logger.info('Min freq: %.2f MHz' % (minfreq/1e6))
            for MS in MSs.getListObj():

                antennaset = MS.getAntennaSet()

                if 'HBA' in antennaset:
                    logger.warning(f'Skipping HBA: deleting {MS.pathMS}')
                    lib_util.check_rm(MS.pathMS)
                    flog.write(MS.nameMS+'.MS\n') # after averaging to be sure no log is written if an error occurs
                    continue

                # get avg time/freq values
                nchan = MS.getNchan()
                timeint = MS.getTimeInt()
                # TODO change these lines to use MS.getAvgFactors() after running survey
                if nchan == 1:
                    avg_factor_f = 1
                #elif nchan % 2 == 0 and MSs.isHBA: # case HBA
                #    avg_factor_f = int(nchan / 4)  # to 2 ch/SB
                elif nchan % 8 == 0 and minfreq < 40e6:
                    avg_factor_f = int(nchan / 8)  # to 8 ch/SB
                elif nchan % 8 == 0 and 'SPARSE' in antennaset:
                    avg_factor_f = int(nchan / 8)  # to 8 ch/SB
                elif nchan % 4 == 0:
                    avg_factor_f = int(nchan / 4)  # to 4 ch/SB
                elif nchan % 5 == 0:
                    avg_factor_f = int(nchan / 5)  # to 5 ch/SB
                else:
                    logger.error('Channels should be a multiple of 4 or 5.')
                    sys.exit(1)

                if keep_IS:
                    # TODO HE - I think we want at least 2s 32ch/SB for IS here (at least for the target field)!
                     avg_factor_f = int(nchan / 16) # to have the full FoV in LBA we need 16 ch/SB
                if avg_factor_f < 1: avg_factor_f = 1

                avg_factor_t = int(np.round(2/timeint)) if keep_IS else int(np.round(4/timeint)) # to 4 sec (2 for IS)
                if avg_factor_t < 1: avg_factor_t = 1

                # do not select IS baselines expect if specified
                bl_sel = '"*&"' if keep_IS else '"[CR]S*&"'

                MSout = getName(MS.pathMS)

                if avg_factor_f != 1 or avg_factor_t != 1 or demix_sources:
                    # normal case: no demixing, just averaging
                    if not demix_sources:
                        logger.info('%s->%s: Average in freq (factor of %i) and time (factor of %i)...' % (MS.nameMS, MSout, avg_factor_f, avg_factor_t))
                        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin={MS.pathMS} msin.baseline={bl_sel} msout={MSout} '
                              f'msin.datacolumn=DATA avg.timestep={avg_factor_t} avg.freqstep={avg_factor_f}',
                              log=MS.nameMS+'_avg.log', commandType='DP3')
                        s.run(check=True, maxProcs=1) # limit threads to prevent I/O isssues
                    # special case: run demix and average in demixer call
                    else:
                        # check HBA - not tested
                        if MSs.isHBA:
                            logger.warning('Demixing in HBA not tested, please be careful!')
                            demix_sm_file = os.path.dirname(__file__) + '/../models/A-Team_lowres.skymodel' if not demix_skymodel else demix_skymodel
                        else:
                            demix_sm_file = os.path.dirname(__file__) + '/../models/demix_all.skymodel' if not demix_skymodel else demix_skymodel
                        demix_sm = lsmtool.load(demix_sm_file, beamMS=MS.pathMS) # load demix sm
                        ra, dec = MS.getPhaseCentre()  # to calculate distance
                        # check target field skymodel
                        # if not provided, use LOTSS DR3. If this is not available, use GSM

                        if demix_field_skymodel:
                            phasecentre = MS.getPhaseCentre()
                            fwhm = MS.getFWHM(freq='min')  # for radius of model
                            if demix_field_skymodel.upper() not in ['GSM','LOTSS','TGSS','VLSSR','NVSS','WENSS','LOTSS-DR3']:
                                sm = lsmtool.load(demix_field_skymodel, beamMS=MS.pathMS)
                            else:
                                if demix_field_skymodel.upper() == 'LOTSS-DR3':
                                    if lib_dd_parallel.check_lotss_coverage(phasecentre, fwhm/2):
                                        logger.info('Target fully in LoTSS-DR3 - start from LoTSS.')
                                        sm = lib_cat.get_LOTSS_DR3_cone_as_skymodel(phasecentre, fwhm / 2,
                                                                                'demixfield_lotss.skymodel', MS.pathMS)
                                    else:
                                        logger.info('Target not fully in LoTSS-DR3 - start from GSM.')
                                        demix_field_skymodel = 'GSM'
                                if demix_field_skymodel.upper() in ['GSM','LOTSS','TGSS','VLSSR','NVSS','WENSS']:
                                    logger.info(f'Include target from {demix_field_skymodel}...')
                                    # get model the size of the image (radius=fwhm/2)
                                    sm = lsmtool.load(demix_field_skymodel, VOPosition=phasecentre, VORadius=fwhm/2, beamMS=MS.pathMS)
                                    sm.remove('I<1')
                                    if demix_field_skymodel.upper() == 'LOTSS':
                                        sm.setColValues('I', sm.getColValues('I')/1000) # convert mJy to Jy TODO fix in LSMtool
                                        sm.setColValues('SpectralIndex', [[-0.7]]*len(sm.getColValues('I'))) # add standard spidx
                            sm.group('single', root='target')
                            sm.setColValues('LogarithmicSI', ['true'] * len(sm))
                            # apply beam to the target-field
                            sm.setColValues('I', sm.getColValues('I', applyBeam=True))
                            # join with ateam skymodel
                            sm.concatenate(demix_sm)
                            # DO NOT USE APPARENT MODEL FOR DEMIX SOURCES! THIS GIVES WORSE RESULTS IN SOME CASES!
                            sm.write(f'{MS.pathMS}/demix_combined.skymodel', clobber=True, applyBeam=False)
                            # lsm.write(f'{MS.pathMS}/demix_combined_apparent.skymodel', clobber=True, applyBeam=True)
                        else: # ignore target
                            logger.info('Ignoring target...')
                            demix_sm.write(f'{MS.pathMS}/demix_combined.skymodel', clobber=True, applyBeam=False)
                            # demix_sm.write(f'{MS.pathMS}/demix_combined_apparent.skymodel', clobber=True, applyBeam=True)
                        # lsm_app = lsmtool.load(f'{MS.pathMS}/demix_combined_apparent.skymodel') # to not calulate beam twice
                        lsm = lsmtool.load(f'{MS.pathMS}/demix_combined.skymodel') # to not calulate beam twice

                        # Get debug info about demix skymodel
                        logger.info('%s->%s: Demix and average in freq (factor of %i) and time (factor of %i)...' % (MS.nameMS, MSout, avg_factor_f, avg_factor_t))
                        demix_sm_patches = demix_sm.getPatchPositions()
                        demix_sm_sources = demix_sm_patches.keys()
                        if not isinstance(demix_sources, list):
                            demix_sources = demix_sources.replace('[','').replace(']','').split(',')
                        # check distances and that all sources in config are actually in skymodel.
                        this_ms_demix_sources = demix_sources.copy()
                        for demix_source in demix_sources:
                            if demix_source not in demix_sm_sources:
                                logger.error(f'demix_source={demix_source} not in demix skymodel sources {demix_sm_sources}!')
                                sys.exit()
                            else:
                                coord_demix = SkyCoord(ra=demix_sm_patches[demix_source][0], dec=demix_sm_patches[demix_source][1])
                                sep = coord_demix.separation(SkyCoord(ra * u.deg, dec * u.deg)).deg
                                logger.info(f'Demix source {demix_source}: sep={sep:.2f}deg')
                                # DEBUG: commented since this is kinda slow and just for information
                                # app_flux = demix_sm.getColValues('I', aggregate='sum', applyBeam=True)[demix_sm.getPatchNames() == demix_source].item() # convert to float
                                # logger.info(f'Demix source {demix_source}: sep={sep:.2f}deg, app. flux={app_flux:.2f}Jy')
                                # flux = demix_sm.getColValues('I', aggregate='sum')[demix_sm.getPatchNames() == demix_source].item() # convert to float
                                # logger.info(f'Demix source {demix_source}: sep={sep:.2f}deg, flux={flux:.2f}Jy')

                        # RUN demixing
                        demix_f_step = nchan
                        demix_t_step = int(np.round(8/timeint))

                        cmd = f"DP3 {parset_dir}/DP3-demix.parset msin={MS.pathMS} msin.baseline={bl_sel} msout={MSout} \
                              demix.skymodel={MS.pathMS}/demix_combined.skymodel demix.instrumentmodel={MS.pathMS}/instrument_demix.parmdb \
                              demix.subtractsources=[{','.join(this_ms_demix_sources)}] \
                              demix.freqstep={avg_factor_f} demix.timestep={avg_factor_t} \
                              demix.demixfreqstep={demix_f_step} demix.demixtimestep={demix_t_step}"
                        if demix_field_skymodel:
                            cmd += ' demix.targetsource=target demix.ignoretarget=False '
                        else:
                            logger.warning('You did not provide a target source. Using ignoretarget=False to deproject (not tested)')
                            cmd += 'demix.ignoretarget=False demix.targetsource="" '
                        s.add(cmd, log=MS.nameMS+'_demix.log', commandType='DP3')
                        s.run(check=True)

                    if backup_full_res:
                        logger.info('Backup full resolution data...')
                        if not os.path.exists('data-bkp'):
                            os.makedirs('data-bkp')
                        MS.move('data-bkp/' + MS.nameMS + '.MS', keepOrig=False, overwrite=False)
                    else:
                        lib_util.check_rm(MS.pathMS)
                    flog.write(MS.nameMS+'.MS\n') # after averaging to be sure no log is written if an error occurs
                else:
                    logger.info('%s->%s: Move data - no averaging...' % (MS.nameMS, MSout))
                    flog.write(MS.nameMS+'.MS\n') # before move or the filenmae is changed
                    MS.move(MSout)

                if tar:
                    logger.info(f"Tar {MSout}...")
                    s.add(f'tar cf {MSout}.tar --directory={os.path.dirname(MSout)} {os.path.basename(MSout)}', log='tar.log', commandType='general')
                    s.run(check=True)
                    lib_util.check_rm(MSout)

w.alldone()
