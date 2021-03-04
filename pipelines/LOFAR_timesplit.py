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
parset_dir = parset.get('LOFAR_timesplit','parset_dir')
data_dir = parset.get('LOFAR_timesplit','data_dir')
cal_dir = parset.get('LOFAR_timesplit','cal_dir')
ngroups = parset.getint('LOFAR_timesplit','ngroups')
initc = parset.getint('LOFAR_timesplit','initc') # initial tc num (useful for multiple observation of same target)
bl2flag = parset.get('flag','stations')

#################################################
# Clean
with w.if_todo('clean'):
    logger.info('Cleaning...')
    lib_util.check_rm('mss*')
### DONE

with w.if_todo('copy'):
    MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

    logger.info('Copy data...')
    for MS in MSs.getListObj():
        #if min(MS.getFreqs()) > 30.e6:
        # overwrite=True to prevent updating the weights twice
        MS.move(MS.nameMS+'.MS', keepOrig=True, overwrite=True)
### DONE

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

##################################################
# Find solutions to apply
if cal_dir == '':
    obsid = MSs.getListObj()[0].getObsID()
    # try standard location
    cal_dir = glob.glob('../id%i_-_3[c|C]196' % obsid)+glob.glob('../id%i_-_3[c|C]295' % obsid)+glob.glob('../id%i_-_3[c|C]380' % obsid)
    if len(cal_dir) > 0:
        cal_dir = cal_dir[0]
    else:
        logger.error('Cannot find solutions.')
        sys.exit()
else:
    cal_dir = '../'+cal_dir

logger.info('Calibrator directory: %s' % cal_dir)
h5_pa = cal_dir+'/cal-pa.h5'
h5_amp = cal_dir+'/cal-amp.h5'
h5_iono = cal_dir+'/cal-iono.h5'
if not os.path.exists(h5_pa) or not os.path.exists(h5_amp) or not os.path.exists(h5_iono):
    logger.error("Missing solutions in %s" % cal_dir)
    sys.exit()

####################################################
# Correct fist for BP(diag)+TEC+Clock and then for beam
with w.if_todo('apply'):
    
    # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
    logger.info('Apply solutions (pa)...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS \
            cor.parmdb='+h5_pa+' cor.correction=polalign', log='$nameMS_cor1.log', commandType='DPPP')
    
    # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
    logger.info('Apply solutions (amp/ph)...')
    MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
            cor.amp.parmdb='+h5_amp+' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
            cor.ph.parmdb='+h5_iono+' cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log', commandType='DPPP')
    
    # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
    logger.info('Beam correction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DPPP')
### DONE

###################################################################################################
# Create groups
groupnames = []
logger.info('Concatenating in frequency...')
for i, msg in enumerate(np.array_split(sorted(glob.glob('*MS')), ngroups)):
   if ngroups == 1:
       groupname = 'mss'
   else:
       groupname = 'mss-%02i' % i
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

       # prepare concatenated mss - SB.MS:CORRECTED_DATA -> group#.MS:DATA (cal corr data, beam corrected)
       s.add('DPPP '+parset_dir+'/DPPP-concat.parset msin="['+','.join(msg)+']"  msout='+groupname+'/'+groupname+'.MS', \
                   log=groupname+'_DPPP_concat.log', commandType='DPPP')
       s.run(check=True)

MSs = lib_ms.AllMSs( glob.glob('mss*/*MS'), s )

#############################################################
# Flagging on concatenated dataset - also flag low-elevation
with w.if_todo('flag'):

    logger.info('Flagging...')
    MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS ant.baseline=\"' + bl2flag + '\" \
            aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.rfis',
            log='$nameMS_DPPP_flag.log', commandType='DPPP')
    
    logger.info('Remove bad timestamps...')
    MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
    
    logger.info('Plot weights...')
    MSs.run('reweight.py $pathMS -v -p -a CS001LBA', log='$nameMS_weights.log', commandType='python')
    lib_util.check_rm('plots-weights')
    os.system('mkdir plots-weights; mv *png plots-weights')
### DONE

#####################################
# Create time-chunks
with w.if_todo('timesplit'):

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
### DONE

logger.info('Cleaning up...')
os.system('rm -r *MS')

logger.info("Done.")
