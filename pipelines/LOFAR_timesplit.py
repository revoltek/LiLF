#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Data preparation for selfcal, apply cal solutions
# and split SB in time and concatenate in frequency.

import sys, os, glob, re
import numpy as np
from astropy.time import Time
import casacore.tables as pt

parset_dir = '/home/fdg/scripts/LiLF/parsets/LOFAR_timesplit'
initc = 0 # initial tc num (useful for multiple observation of same target)

# temporary
if 'LBAsurvey' in os.getcwd():
    ngroups = 1 # number of groups (totalSB/SBperFREQgroup)
    datadir = '../../download/%s/%s' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])
    soldir = 'dsk:/disks/paradata/fdg/LBAsurvey/cal_'+os.getcwd().split('/')[-2]
elif 'bootes2' in os.getcwd():
    ngroups = 1
    datadir = '../tgts-bkp/' 
    soldir = '../cals/'
else:
    ngroups = 2
    datadir = '../tgts-bkp/' 
    soldir = '../cals/'

assert os.path.isdir(soldir)

########################################################
from LiLF import lib_ms, lib_util, lib_log
lib_log.set_logger('pipeline-timesplit.logger')
logger = lib_log.logger
lib_util.check_rm('logs')
s = lib_util.Scheduler(dry = False)

#################################################
# Clear
logger.info('Cleaning...')
lib_util.check_rm('mss*')
MSs = lib_ms.AllMSs( glob.glob(datadir+'/*MS'), s )

##############################################
# Avg to 4 chan and 4 sec
# Remove internationals
# TODO: move to download pipeline
#nchan = find_nchan(mss[0])
#timeint = find_timeint(mss[0])
#if nchan % 4 != 0 and nchan != 1:
#    logger.error('Channels should be a multiple of 4.')
#    sys.exit(1)
#
#avg_factor_f = nchan / 4 # to 4 ch/SB
#if avg_factor_f < 1: avg_factor_f = 1
#avg_factor_t = int(np.round(4/timeint)) # to 4 sec
#if avg_factor_t < 1: avg_factor_t = 1
#
#if avg_factor_f != 1 or avg_factor_t != 1:
#    logger.info('Average in freq (factor of %i) and time (factor of %i)...' % (avg_factor_f, avg_factor_t))
#    for ms in mss:
#        msout = ms.replace('.MS','-avg.MS').split('/')[-1]
#        if os.path.exists(msout): lib_util.check_rm(ms)
#        s.add('DPPP '+parset_dir+'/DPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.timestep='+str(avg_factor_t)+' avg.freqstep='+str(avg_factor_f), \
#                log=msout+'_avg.log', commandType='DPPP')
#    s.run(check=True)
#    nchan = nchan / avg_factor_f
#    timeint = timeint * avg_factor_t
#else:
logger.info('Copy data...')
for MS in MSs.getListObj():
    MS.move(MS.nameMS+'.MS', keepOrig=True)

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

####################################################
# Correct fist for BP(diag)+TEC+Clock and then for beam
# Copy instrument tables
# TODO: this can now be made AFTER grouping
logger.info('Copy solutions...')
if not os.path.exists('cal-pa.h5'):
    os.system('scp -q '+soldir+'/cal-pa.h5 .')
if not os.path.exists('cal-amp.h5'):
    os.system('scp -q '+soldir+'/cal-amp.h5 .')
if not os.path.exists('cal-iono.h5'):
    os.system('scp -q '+soldir+'/cal-iono.h5 .')

# Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
logger.info('Apply solutions...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.steps=[pa] \
        cor.pa.parmdb=cal-pa.h5 cor.pa.correction=polalign', log='$nameMS_cor1.log', commandType='DPPP')

# Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
logger.info('Beam correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DPPP')

# Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
logger.info('Apply solutions...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
        cor.amp.parmdb=cal-amp.h5 cor.amp.correction=amplitudeSmooth000 cor.amp.updateweights=True\
        cor.ph.parmdb=cal-iono.h5 cor.ph.correction=clock000', log='$nameMS_cor2.log', commandType='DPPP') # TODO: clock?
#MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
#        cor.amp.parmdb=cal-amp.h5 cor.amp.correction=amplitudeSmooth000 cor.amp.updateweights=True\
#        cor.ph.parmdb=cal-iono.h5 cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log', commandType='DPPP')

# Re-set weights of flagged data to 0 - this is necessary if we want to use dysco
# due to DPPP leaving flagged weight at super-high values compared to unflagged ones
#logger.info('Set weight of flagged data to 0...')
#MSs.run('flag_weight_to_zero.py $pathMS', log='$nameMS_resetweight.log', commandType='python')

###################################################################################################
# Create groups
groupnames = []
logger.info('Concatenating in frequency...')
timechunks = set([re.findall(r'_t\d+', ms)[0][2:] for ms in MSs.getListStr() ])
for timechunk in timechunks:
    for i, msg in enumerate(np.array_split(sorted(glob.glob('*_t'+timechunk+'_*MS')), ngroups)):
        if ngroups == 1:
            groupname = 'mss_t%s' % timechunk
        else:
            groupname = 'mss_t%s-%02i' % (timechunk, i)
        groupnames.append(groupname)
        lib_util.check_rm(groupname)
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

# Flagging on concatenated dataset - also flag low-elevation
logger.info('Flagging...')
MSs = lib_ms.AllMSs( glob.glob('mss_t*/*MS'), s )
MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS', \
                log='$nameMS_DPPP_flag.log', commandType='DPPP')

# Create time-chunks
logger.info('Splitting in time...')
for groupname in groupnames:
    ms = groupname+'/'+groupname+'.MS'
    if not os.path.exists(ms): continue
    t = pt.table(ms, ack=False)
    starttime = t[0]['TIME']
    endtime   = t[t.nrows()-1]['TIME']
    hours = (endtime-starttime)/3600.
    logger.debug(ms+' has length of '+str(hours)+' h.')

    tc = initc
    for timerange in np.array_split(t.getcol('TIME'), round(hours)):
        logger.info('%02i - Splitting timerange %f %f' % (tc, timerange[0], timerange[-1]))
        t1 = t.query('TIME >= ' + str(timerange[0]) + ' && TIME <= ' + str(timerange[-1]), sortlist='TIME,ANTENNA1,ANTENNA2')
        splitms = groupname+'/TC%02i.MS' % tc
        lib_util.check_rm(splitms)
        t1.copy(splitms, True)
        t1.close()
        tc += 1
    t.close()

    lib_util.check_rm(ms) # remove not-timesplitted file

# put everything together
if ngroups == 1:
    lib_util.check_rm('mss')
    os.makedirs('mss')
    os.system('mv mss_t*/*MS mss')
    lib_util.check_rm('mss_t*')
else:
    for group in xrange(ngroups):
        groupname = 'mss-%02i' % group
        lib_util.check_rm(groupname)
        os.makedirs(groupname)
        os.system('mv mss_t*-%02i/*MS %s' % (group, groupname))
    lib_util.check_rm('mss_t*')

logger.info('Cleaning up...')
os.system('rm -r *MS')

logger.info("Done.")
