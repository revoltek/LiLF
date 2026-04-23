#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, re, glob

import lsmtool
import numpy as np
import casacore.tables as pt
from astropy.coordinates import SkyCoord
import astropy.units as u

##########################################
from LiLF import lib_ms, lib_util, lib_log, lib_dd_parallel, lib_cat, lib_scheduler, lib_walker


def run(step):

    log_dir = lib_log.Logger(f'pipeline-{step.kind}-{step.name}').log_dir
    s = lib_scheduler.Scheduler(dry_run=False, log_dir=log_dir)
    w = lib_walker.Walker(f'pipeline-{step.kind}-{step.name}.walker')
    logger = lib_log.logger

    input_dir = step['input']
    output_dir = step['output']
    parset_dir = step['parset_dir']
    fix_table = step['fix_table']
    keep_IS = step['keep_IS']
    backup_full_res = step['backup_full_res']
    demix_sources = step['demix_sources']
    demix_skymodel = step['demix_skymodel']
    demix_field_skymodel = step['demix_field_skymodel']
    run_aoflagger = step['run_aoflagger']
    tar = step['tar']
    validate = step['validate']


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
    
        # remove unnecessary info for survey pointings
        m = re.search(r'[Pp](\d{3}\+\d{2})', code)
        if m:
            code = 'P' + m.group(1)
    
        # get obsid
        #with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
        #    obsid = t.getcell('LOFAR_OBSERVATION_ID',0)

        # get freq
        with pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False) as t:
            freq = t.getcell('REF_FREQUENCY',0)
    
        return f'{output_dir}/{code}_SB{lib_util.lofar_nu2num(freq / 1e6):03d}.MS'

    MSs = lib_ms.AllMSs(glob.glob(os.path.join(input_dir, '*MS')), s, check_flags=False, check_consistency=validate)
    if len(MSs.getListStr()) == 0:
        raise RuntimeError(f"No MSs found in input_dir '{input_dir}'.")

    if keep_IS and not MSs.hasIS:
        logger.warning('Keeping IS requested but IS were not recorded - switching to averaging settings for Dutch only.')
        keep_IS = False

  ######################################
    time = int(MSs.getListObj()[0].get_time().iso.replace('-', '')[:8])
    antennaset = MSs.getListObj()[0].getAntennaSet()
    minfreq = np.min(MSs.getFreqs())
    nchan = MSs.getListObj()[0].getNchan()
    timeint = MSs.getListObj()[0].getTimeInt()
    logger.info(f'Starting from: nchan={nchan}, timeint={timeint:.2f}s (antennaset={antennaset})')

    if run_aoflagger:
        with w.if_todo('flag'):
            # Flag in an identical way to the observatory flagging
            logger.info('Flagging...')
            MSs.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS aoflagger.strategy={parset_dir}/LBAdefaultwideband.lua',
                log='$nameMS_flag.log', commandType='DP3')

    if fix_table:
        with w.if_todo('fix_table'):
            logger.info('Fix MS table...')
            # only ms created in range (2/2013->2/2014)
            if 20130200 < time < 20140300:
                logger.info('Fix beam table...')
                MSs.run('fixbeaminfo ' + MS, log='fixbeam.log', commandType='python')

    # Rescale visibilities by 1e3 if before 2014-03-19 (old correlator), and by 1e-2 otherwise
    with w.if_todo('rescale_flux'):
        logger.info('Rescaling flux...')
        if time < 20140319:
            rescale_factor = 1e6
        else:
            rescale_factor = 1e-4
        
        for i, MS in enumerate(MSs.getListStr()):
            with pt.table(f'{MS}/HISTORY', readonly=False, ack=False) as hist:
                if 'Flux rescaled' not in hist.getcol('MESSAGE'):
                    s.add(f'taql "update {MS} set DATA = {rescale_factor}*DATA" && '
                          f'taql "insert into {MS}/HISTORY (TIME,MESSAGE) values (mjd(), \'Flux rescaled\')"',
                          log='taql.log', commandType='general')
        s.run(check=True)

    ######################################
    # Averaging/demixing/removing IS (if requested)
    with w.if_todo('renameavg'):
        logger.info('Renaming/averaging...')
        # Build a getName() cache to avoid opening casacore tables twice per MS
        # (once in the filter and once in the loop body).
        _name_cache = {ms: getName(ms) for ms in glob.glob(os.path.join(input_dir, '*MS'))}
        with open('renamed.txt','a') as flog:
            MSs = lib_ms.AllMSs([ms for ms, dst in _name_cache.items() if not os.path.exists(dst)], s, check_flags=False)
            for MS in MSs.getListObj():
                if 'HBA' in antennaset:
                    logger.warning(f'Skipping HBA: deleting {MS.pathMS}')
                    lib_util.check_rm(MS.pathMS)
                    flog.write(MS.nameMS+'.MS\n') # after averaging to be sure no log is written if an error occurs
                    continue

                # get avg time/freq values
                logger.info(f'Processing {MS.nameMS}...')
                avg_factor_t, avg_factor_f = MS.getAvgFactors(keep_IS, minfreq)

                # do not select IS baselines expect if specified
                bl_sel = '"*&"' if keep_IS else '"[CR]S*&"'

                MSout = _name_cache[MS.pathMS]

                if avg_factor_f != 1 or avg_factor_t != 1 or demix_sources:
                    # normal case: no demixing, just averaging
                    if not demix_sources:
                        logger.info(f'{MS.nameMS}->{MSout}: Average in freq (x{avg_factor_f}) and time (x{avg_factor_t})...')
                        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin={MS.pathMS} msin.baseline={bl_sel} msout={MSout} '
                                f'msin.datacolumn=DATA avg.timestep={avg_factor_t} avg.freqstep={avg_factor_f}',
                                log=MS.nameMS+'_avg.log', commandType='DP3')
                        s.run(check=True, max_proc=1) # limit threads to prevent I/O isssues
                    else:
                        # special case: run demix and average in demixer call
                        _models_dir = os.path.join(os.path.dirname(__file__), '..', 'models')
                        demix_sm_file = demix_skymodel or os.path.join(_models_dir, 'demix_all.skymodel')
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
                        else: # ignore target
                            logger.info('Ignoring target...')
                            demix_sm.write(f'{MS.pathMS}/demix_combined.skymodel', clobber=True, applyBeam=False)

                        # Get debug info about demix skymodel
                        logger.info(f'{MS.nameMS}->{MSout}: Demix and average in freq (x{avg_factor_f}) and time (x{avg_factor_t})...')
                        demix_sm_patches = demix_sm.getPatchPositions()
                        demix_sm_sources = demix_sm_patches.keys()
                        # Normalise demix_sources to a list once, before the loop.
                        if not isinstance(demix_sources, list):
                            demix_sources = [s.strip() for s in
                                                demix_sources.strip('[]').split(',')]
                        # check distances and that all sources in config are actually in skymodel.
                        this_ms_demix_sources = demix_sources.copy()
                        for demix_source in demix_sources:
                            if demix_source not in demix_sm_sources:
                                raise RuntimeError(
                                    f'demix_source={demix_source!r} not found in demix '
                                    f'skymodel sources: {list(demix_sm_sources)}')
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
                        os.makedirs(f'{output_dir}/data-bkp', exist_ok=True)
                        MS.move(f'{output_dir}/data-bkp/{MS.nameMS}.MS', keepOrig=False, overwrite=False)
                    else:
                        lib_util.check_rm(MS.pathMS)
                    flog.write(MS.nameMS+'.MS\n') # after averaging to be sure no log is written if an error occurs
                else:
                    logger.info(f'{MS.nameMS}.MS->{MSout}: Move data - no averaging...')
                    flog.write(MS.nameMS+'.MS\n')  # before move or the filename is changed
                    MS.move(MSout)

                if tar:
                    logger.info(f"Tar {MSout}...")
                    s.add(f'tar cf {MSout}.tar --directory={os.path.dirname(MSout)} {os.path.basename(MSout)}', log='tar.log', commandType='general')
                    s.run(check=True)
                    lib_util.check_rm(MSout)

    w.alldone()
