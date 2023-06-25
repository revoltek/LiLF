#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, re, glob, time
import numpy as np
import casacore.tables as pt
from astropy.time import Time

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

                MSout = getName(MS.pathMS)
                if avg_factor_f != 1 or avg_factor_t != 1:
                    logger.info('%s->%s: Average in freq (factor of %i) and time (factor of %i)...' % (MS.nameMS, MSout, avg_factor_f, avg_factor_t))
                    if keep_IS:
                        s.add('DP3 '+parset_dir+'/DP3-avg.parset msin='+MS.pathMS+' msout='+MSout+' msin.datacolumn=DATA \
                            avg.timestep='+str(avg_factor_t)+' avg.freqstep='+str(avg_factor_f), \
                            log=MS.nameMS+'_avg.log', commandType='DP3')
                    else: # remove IS
                        s.add('DP3 '+parset_dir+'/DP3-avg.parset msin='+MS.pathMS+' msout='+MSout+' msin.datacolumn=DATA \
                            msin.baseline="[CR]S*&" \
                            avg.timestep='+str(avg_factor_t)+' avg.freqstep='+str(avg_factor_f), \
                            log=MS.nameMS+'_avg.log', commandType='DP3')
                    s.run(check=True, maxThreads=1) # limit threads to prevent I/O isssues
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
