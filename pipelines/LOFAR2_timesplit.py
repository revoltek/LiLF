#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Data preparation for selfcal, apply cal solutions
# and split SB in time and concatenate in frequency.

import sys, os, glob, re
import numpy as np
from astropy.time import Time
import casacore.tables as pt


########################################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-timesplit.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-timesplit.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR2_timesplit','parset_dir')
data_dir = parset.get('LOFAR2_timesplit','data_dir')
cal_dir = parset.get('LOFAR2_timesplit','cal_dir')
ngroups = parset.getint('LOFAR2_timesplit','ngroups')
initc = parset.getint('LOFAR2_timesplit','initc') # initial tc num (useful for multiple observation of same target)
bl2flag = parset.get('flag','stations')

#################################################
# Clean
if w.todo('clean'):
    logger.info('Cleaning...')
    lib_util.check_rm('mss*')

    w.done('clean')
### DONE

if w.todo('copy'):
    MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

    logger.info('Copy data...')
    for MS in MSs.getListObj():
        if min(MS.getFreqs()) > 30.e6:
            # overwrite=True to prevent updating the weights twice
            MS.move(MS.nameMS+'.MS', keepOrig=True, overwrite=True)

    w.done('copy')
### DONE

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

##################################################
# Find solutions to apply
h5_pa = cal_dir+'/cal-pa.h5'
h5_amp = cal_dir+'/cal-amp.h5'

assert os.path.exists(h5_pa)
assert os.path.exists(h5_amp)

####################################################
# Correct fist for BP(diag)+TEC+Clock and then for beam
if w.todo('apply'):
    
    # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign, amp, beam corrected)
    logger.info('Apply solutions (pa,amp,beam)...')
    MSs.run('DPPP '+parset_dir+'/DPPP-apply.parset msin=$pathMS \
            cor.pa.parmdb='+h5_pa+' cor.amp.parmdb='+h5_amp, log='$nameMS_apply.log', commandType='DPPP')
    

    w.done('apply')
### DONE

###################################################################################################
# Create groups
# TODO: the creation of groups should always be:
# - 1 group with the most sensitive region
# - 1 group below that region
# - 1 group above that region
# to be understood if calibration/imaging gain anything from using only the central group or all of them
groupnames = []
logger.info('Concatenating in frequency...')
timechunks = sorted(set([re.findall(r'_t[0-9]+', ms)[0][2:] for ms in MSs.getListStr() ]))

for timechunk in timechunks:
    for i, msg in enumerate(np.array_split(sorted(glob.glob('*_t'+timechunk+'_*.MS')), ngroups)):
        logger.info(i, msg)
        if ngroups == 1:
            groupname = 'mss_t%s' % timechunk
        else:
            groupname = 'mss_t%s-%02i' % (timechunk, i)
        groupnames.append(groupname)

        # skip if already done
        if not os.path.exists(groupname):
            os.makedirs(groupname)
    
            # add missing SB with a fake name not to leave frequency holes
            num_init = int(re.findall(r'\d+', msg[0])[-1])
            num_fin = int(re.findall(r'\d+', msg[-1])[-1])
            ms_name_init = msg[0]
            msg = []
            for j in range(num_init, num_fin+1):
                msg.append(ms_name_init.replace('SB%03i' % num_init, 'SB%03i' % j))
    
            # prepare concatenated time chunks (TC) - SB.MS:CORRECTED_DATA -> group#.MS:DATA (cal corr data, beam corrected, circular)
            s.add('DPPP '+parset_dir+'/DPPP-concat.parset msin="['+','.join(msg)+']"  msout='+groupname+'/'+groupname+'.MS', \
                        log=groupname+'_DPPP_concat.log', commandType='DPPP')
    s.run(check=True)

MSs = lib_ms.AllMSs( glob.glob('mss_t*/*MS'), s )

#############################################################
# Flagging on concatenated dataset - also flag low-elevation
if w.todo('flag'):

    logger.info('Flagging...')
    # MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS ant.baseline=\"' + bl2flag + '\"', \
    #                 log='$nameMS_DPPP_flag.log', commandType='DPPP')
    
    logger.info('Remove bad timestamps...')
    MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
    
    logger.info('Plot weights...')
    MSs.run('reweight.py $pathMS -v -p -a CS001LBA', log='$nameMS_weights.log', commandType='python')
    os.system('mkdir plots-weights; mv *png plots-weights')

    w.done('flag')
### DONE

#sys.exit() # for DDFacet

#####################################
# Create time-chunks
if w.todo('timesplit'):

    logger.info('Splitting in time...')
    tc = initc
    for groupname in groupnames:
        ms = groupname+'/'+groupname+'.MS'
        if not os.path.exists(ms): continue
        t = pt.table(ms, ack=False)
        starttime = t[0]['TIME']
        endtime   = t[t.nrows()-1]['TIME']
        hours = (endtime-starttime)/3600.
        logger.debug(ms+' has length of '+str(hours)+' h.')
    
        for timerange in np.array_split(sorted(set(t.getcol('TIME'))), round(hours)):
            logger.info('%02i - Splitting timerange %f %f' % (tc, timerange[0], timerange[-1]))
            t1 = t.query('TIME >= ' + str(timerange[0]) + ' && TIME <= ' + str(timerange[-1]), sortlist='TIME,ANTENNA1,ANTENNA2')
            splitms = groupname+'/TC%02i.MS' % tc
            lib_util.check_rm(splitms)
            t1.copy(splitms, True)
            t1.close()
            tc += 1
        t.close()
    
        lib_util.check_rm(ms) # remove not-timesplitted file

    w.done('timesplit')
### DONE

############################################
# put everything together
if w.todo('concat'):

    if ngroups == 1:
        lib_util.check_rm('mss')
        os.makedirs('mss')
        os.system('mv mss_t*/*MS mss')
        lib_util.check_rm('mss_t*')
    else:
        for group in range(ngroups):
            groupname = 'mss-%02i' % group
            lib_util.check_rm(groupname)
            os.makedirs(groupname)
            os.system('mv mss_t*-%02i/*MS %s' % (group, groupname))
        lib_util.check_rm('mss_t*')
    
    logger.info('Cleaning up...')
    os.system('rm -r *MS')

    w.done('concat')
### DONE

logger.info("Done.")
