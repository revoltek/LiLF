#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Data preparation for selfcal, apply cal solutions
# and split SB in time and concatenate in frequency.

import sys, os, glob, re
import numpy as np
from astropy.time import Time
import casacore.tables as pt


if 'LBAsurvey' in os.getcwd():
    datadir = '/home/baq1889/lofar1/LBAsurvey/%s/%s' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])
    cal_dir = 'portal_lei:/disks/paradata/fdg/LBAsurvey/cal_'+os.getcwd().split('/')[-2]

########################################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-timesplit.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('timesplit','parset_dir')
data_dir = parset.get('timesplit','data_dir')
cal_dir = parset.get('timesplit','cal_dir')
ngroups = parset.getint('timesplit','ngroups')
initc = parset.getint('timesplit','initc') # initial tc num (useful for multiple observation of same target)

assert os.path.isdir(cal_dir)

#################################################
# Clean
logger.info('Cleaning...')
lib_util.check_rm('mss*')
MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

logger.info('Copy data...')
for MS in MSs.getListObj():
    if min(MS.getFreqs()) > 30.e6:
        MS.move(MS.nameMS+'.MS', keepOrig=True)

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

####################################################
# Correct fist for BP(diag)+TEC+Clock and then for beam
# Copy instrument tables
logger.info('Copy solutions...')
if not os.path.exists('cal-pa.h5'):
    os.system('scp -q '+cal_dir+'/cal-pa.h5 .')
if not os.path.exists('cal-amp.h5'):
    os.system('scp -q '+cal_dir+'/cal-amp.h5 .')
if not os.path.exists('cal-iono.h5'):
    os.system('scp -q '+cal_dir+'/cal-iono.h5 .')

# Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
logger.info('Apply solutions...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.steps=[pa] \
        cor.pa.parmdb=cal-pa.h5 cor.pa.correction=polalign', log='$nameMS_cor1.log', commandType='DPPP')

# Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
#logger.info('Apply solutions...')
#MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
#        cor.amp.parmdb=cal-amp.h5 cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
#        cor.ph.parmdb=cal-iono.h5 cor.ph.correction=clock000', log='$nameMS_cor2.log', commandType='DPPP') # TODO: clock?
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
        cor.amp.parmdb=cal-amp.h5 cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
        cor.ph.parmdb=cal-iono.h5 cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log', commandType='DPPP')

# Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
logger.info('Beam correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DPPP')

###################################################################################################
# Create groups
# TODO: the creation of groups should always be:
# - 1 group with the most sensitive region
# - 1 group below that region
# - 1 group above that region
# to be understood if calibration/imaging gain anything from using only the central group or all of them
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

#logger.info('Remove bad timestamps...')
#MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

logger.info('Plot weights...')
MSs.run('reweight.py $pathMS -v -p -a CS001LBA', log='$nameMS_weights.log', commandType='python')
os.system('mkdir plots-weights; mv *png plots-weights')

#sys.exit() # for DDFacet

# Create time-chunks
logger.info('Splitting in time...')
for groupname in groupnames:
    tc = initc
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
