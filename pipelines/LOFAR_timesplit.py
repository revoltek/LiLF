#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Data preparation for selfcal, apply cal solutions
# and split SB in time and concatenate in frequency.

import sys, os, glob, re
import numpy as np
import casacore.tables as pt

########################################################
from LiLF import lib_ms, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-timesplit')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-timesplit.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_timesplit'])))
parset_dir = parset.get('LOFAR_timesplit','parset_dir')
data_dir = parset.get('LOFAR_timesplit','data_dir')
cal_dir = parset.get('LOFAR_timesplit','cal_dir')
fillmissingedges = parset.getboolean('LOFAR_timesplit', 'fillmissingedges')
ngroups = parset.getint('LOFAR_timesplit','ngroups')
initc = parset.getint('LOFAR_timesplit','initc') # initial tc num (useful for multiple observation of same target)
apply_fr = parset.getboolean('LOFAR_timesplit','apply_fr') # also transfer the FR solutions (possibly useful if calibrator and target are close, especially for IS data.)
no_aoflagger = parset.getboolean('LOFAR_timesplit','no_aoflagger')
bl2flag = parset.get('flag','stations')
use_GNSS = parset.getboolean('LOFAR_timesplit', 'use_GNSS')
#################################################

# Clean
with w.if_todo('clean'):
    logger.info('Cleaning...')
    mss_list = glob.glob('mss*/*MS')
    if len(mss_list) > 0:
        raise ValueError(f'mss folders exist already {mss_list}! If this is the output of a previous LOFAR_timesplit.py run and you want to re-run LOFAR_timesplit.py, then delete them manually.')
### DONE

with w.if_todo('copy'):
    #if not os.path.exists(data_dir):
    #    os.system('mkdir '+data_dir)
    #    for ms in sorted(glob.glob('*.MS')):
    #        os.system(f'cp -r {ms} {data_dir}')
    for tarfile in glob.glob(data_dir + '/*tar'):
        if not os.path.exists(tarfile.replace('.tar','')):
            s.add(f'tar xf {tarfile} --one-top-level={data_dir}', log='tar.log', commandType='general')
    if len(s.action_list) > 0:
        logger.info('Untar files...')
        s.run(check=True, maxProcs=5)

    MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )

    logger.info('Copy data...')
    for MS in MSs.getListObj():
        # if min(MS.getFreqs()) > 30.e6:
        # overwrite=True to prevent updating the weights twice
        MS.move(MS.nameMS+'.MS', keepOrig=True, overwrite=True)
### DONE

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )

##################################################
# Find solutions to apply
obsid = MSs.getListObj()[0].getObsID()
if cal_dir == '':
    # try standard location
    cal_dir = glob.glob('../id%i_-_*3[c|C]196' % obsid)+glob.glob('../id%i_-_*3[c|C]295' % obsid)+glob.glob('../id%i_-_*3[c|C]380' % obsid)
    if len(cal_dir) > 0:
        cal_dir = cal_dir[0]
    else:
        logger.error('Cannot find solutions.')
        sys.exit()
else:
    if not cal_dir[0] == '/': # if not abolute path
        cal_dir = '../'+cal_dir
    # cal_dir can either be a path to the directory containing multiple calibrator observation or one specifies an exact directory.
    subdirs = glob.glob(f'{cal_dir}/id{obsid}_-_*3[c|C]196')+glob.glob(f'{cal_dir}/id{obsid}_-_*3[c|C]295')+glob.glob(f'{cal_dir}/id{obsid}_-_*3[c|C]380')
    if len(subdirs) > 0:
        logger.warning('Multiple cal dirs found (using the first):', subdirs)
        cal_dir = subdirs[0]

logger.info('Calibrator directory: %s' % cal_dir)
h5_pa = cal_dir+'/cal-pa.h5'
h5_bp = cal_dir+'/cal-bp.h5'
h5_iono = cal_dir+'/cal-iono.h5'
h5_iono_cs = cal_dir+'/cal-iono-cs.h5'
if not os.path.exists(h5_pa) or not os.path.exists(h5_bp) or not os.path.exists(h5_iono) or not os.path.exists(h5_iono_cs):
    logger.error("Missing solutions in %s" % cal_dir)
    sys.exit()

####################################################
# Correct fist for PA+beam+BP(diag)+TEC+Clock and then for beam
with w.if_todo('apply'):

    # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
    logger.info('Apply solutions (pa)...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
            cor.parmdb={h5_pa} cor.correction=polalign', log='$nameMS_corPA.log', commandType='DP3')
    
    # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
    logger.info('Beam correction...')
    MSs.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_corBEAM.log', commandType='DP3')
    if use_GNSS:
        # Correct gps-tec concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('TEC correction (GPS)...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb={cal_dir}/cal-gps-tec.h5 \
                      cor.correction=tec000', log='$nameMS_cor-gps-tec.log', commandType="DP3")
        # Correct TEC concat_all:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('dTEC correction (fitted)...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb={cal_dir}/cal-dtec.h5 \
                      cor.correction=tec000', log='$nameMS_cor-dtec.log', commandType="DP3")
        # Correct FR concat_all.MS:CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Faraday rotation pre-correction (GPS)...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb={cal_dir}/cal-gps-rm.h5 \
                        cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
    # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, beam corrected+reweight, calibrator corrected+reweight)
    logger.info('Iono correction...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={h5_iono_cs} msin.datacolumn=CORRECTED_DATA \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={h5_iono} msin.datacolumn=CORRECTED_DATA \
                cor.correction=phase000', log='$nameMS_corIONO.log', commandType="DP3")

    logger.info('BP correction...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={h5_bp} msin.datacolumn=CORRECTED_DATA \
                    cor.correction=amplitudeSmooth cor.updateweights=True',
                    log='$nameMS_corBP.log', commandType="DP3")
    if apply_fr:
        logger.info('FR correction...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={cal_dir}/cal-fr.h5 msin.datacolumn=CORRECTED_DATA \
                        cor.correction=rotationmeasure000 ', log='$nameMS_corBP.log', commandType="DP3")

### DONE

###################################################################################################
# Create groups
groupnames = []
for i, msg in enumerate(np.array_split(sorted(glob.glob('*MS')), ngroups)):
    if ngroups == 1:
        groupname = 'mss'
    else:
        groupname = 'mss-%02i' % i
    if MSs.hasIS:
        groupname += '-IS'
    groupnames.append(groupname)

    # skip if already done
    if not os.path.exists(groupname):
        logger.info('Concatenating in frequency...')
        os.makedirs(groupname)

        if fillmissingedges:
            # add missing SB with a fake name not to leave frequency holes
            min_nu = pt.table(MSs.getListStr()[0], ack=False).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MIN']
            max_nu = pt.table(MSs.getListStr()[0], ack=False).OBSERVATION[0]['LOFAR_OBSERVATION_FREQUENCY_MAX']
        else:
            min_nu = min(MSs.getFreqs())/1e6
            max_nu = max(MSs.getFreqs())/1e6
        num_init = lib_util.lofar_nu2num(min_nu)+1  # +1 because FREQ_MIN/MAX somewhat have the lowest edge of the SB freq
        num_fin = lib_util.lofar_nu2num(max_nu)+1
        prefix = re.sub('SB[0-9]*.MS','',msg[0])
        msg = []
        for j in range(num_init, num_fin+1):
            msg.append(prefix+'SB%03i.MS' % j)

        # prepare concatenated mss - SB.MS:CORRECTED_DATA -> group#.MS:DATA (cal corr data, beam corrected)
        s.add('DP3 '+parset_dir+'/DP3-concat.parset msin="['+','.join(msg)+']" msin.missingdata=True msin.orderms=False \
               msout='+groupname+'/'+groupname+'-temp.MS', log=groupname+'_DP3_concat.log', commandType='DP3')
        s.run(check=True)

        # We need a number of channels that is - after averaging to the final dutch wide-field resolution - divisable by 48.check that nchan is divisible by 48 - necessary in dd pipeline; discard high freq unused channels
        nchan_init = MSs.getListObj()[0].getNchan()*len(msg)
        final_freqres_dutch = 0.048828e6 if 'OUTER' in MSs.getListObj()[0].getAntennaSet() else 0.024414e6
        freqres = MSs.getListObj()[0].getChanband()
        averaging_factor = int(round(final_freqres_dutch / freqres))
        nchan = nchan_init - nchan_init % 48*averaging_factor
        logger.info('Reducing total channels: %ich -> %ich)' % (nchan_init, nchan))
        s.add(f'DP3 {parset_dir}/DP3-concat.parset msin={groupname}/{groupname}-temp.MS msin.datacolumn=DATA msin.nchan={nchan} msout={groupname}/{groupname}.MS',
              log=groupname+'_DP3_concat.log', commandType='DP3')
        s.run(check=True)
        # delete temporary MSs:
        os.system(f'rm -r {groupname}/{groupname}-temp.MS')

MSs = lib_ms.AllMSs( glob.glob('mss*/*MS'), s )

#############################################################
# Flagging on concatenated dataset - also flag low-elevation
with w.if_todo('flag'):
    logger.info('Flagging...')
    flag_strat = '/HBAdefaultwideband.lua' if MSs.isHBA else '/LBAdefaultwideband.lua'
    MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS ant.baseline=\"' + bl2flag + '\" \
            aoflagger.strategy='+parset_dir+flag_strat,
            log='$nameMS_DP3_flag.log', commandType='DP3')

    if MSs.hasIS:
        # flagonmindata code cannot handle larger MS - anyway shouldn't be needed if we have the minvisratio in the solves?
        logger.warning('Skip flagonmindata for data with IS present...')
    else:
        logger.info('Remove bad timestamps...')
        MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    try:
        logger.info('Plot weights...')
        MSs.run('reweight.py $pathMS -v -p -a %s' % (MSs.getListObj()[0].getAntennas()[0]),
                log='$nameMS_weights.log', commandType='python')
        lib_util.check_rm('plots-weights')
        os.system('mkdir plots-weights; mv *png plots-weights')
    except RuntimeError:
        logger.warning('Plotting weights failed... continue.')
### DONE

#####################################
# Create time-chunks
with w.if_todo('timesplit'):

    logger.info('Splitting in time...')
    tc = initc
    for groupname in groupnames:
        ms = groupname + '/' + groupname + '.MS'
        if not os.path.exists(ms):
            continue

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
        lib_util.check_rm(ms)  # remove not-timesplitted file
### DONE
# If we have IS present, also split out the averaged dutch baselines for DDparallel processing
if MSs.hasIS:
    for groupname in groupnames:
        MSs = lib_ms.AllMSs( glob.glob(groupname+'/*MS'), s )
        groupname_dutch = groupname.replace("-IS","")
        with w.if_todo('avgdutch'):
            if not os.path.exists(groupname_dutch):
                os.system(f'mkdir {groupname_dutch}')
            MS = lib_ms.AllMSs([glob.glob(data_dir+'/*MS')[0]], s).getListObj()[0]
            avg_factor_t, avg_factor_f = MS.getAvgFactors(keep_IS=False)
            MSs.run(f'DP3 {parset_dir}/DP3-avgdutch.parset msin=$pathMS msout={groupname_dutch}/$nameMS.MS avg.freqstep={avg_factor_f} avg.timestep={avg_factor_t}',
                              log=MS.nameMS+'_avg.log', commandType='DP3')

logger.info('Cleaning up...')
os.system('rm -r *MS')

w.alldone()
