#!/usr/bin/env python3

import os, sys, pickle, glob, re
import numpy as np
from losoto import h5parm

# get input
caldirs = sorted(glob.glob('id*'))
badobs_picklefile = '/homes/fdg/storage/scripts/survey/LoLSS-bkp/badobsids.pickle'
re_sun_pattern = r".*Sun distance: (\d+).*"
re_jup_pattern = r".*Jupiter distance: (\d+).*"
re_elev_pattern = r".*Elev: (\d+).*"

allvalstd = []; allvaldynrange = []; alltimes = [];  allids = []; badids = []; allsundist = []; alljupdist = []; allelev = []
for caldir in caldirs:

    loggerfile = glob.glob(caldir+'/pipeline-cal*logger')[0]
    with open(loggerfile, "r") as f:
        for ln in f:
            if "Sun distance" in ln:
                match_sun = re.search(re_sun_pattern, ln)
                sun_dist = int(match_sun.group(1))  # Extract the first captured group

                match_jup = re.search(re_jup_pattern, ln)
                jup_dist = int(match_jup.group(1))  # Extract the first captured group

                match_elev = re.search(re_elev_pattern, ln)
                elev = int(match_elev.group(1))  # Extract the first captured group
 
    h5parmFile = caldir+'/cal-bp.h5'
    if not os.path.exists(h5parmFile):
        print(f'{h5parmFile} - Not completed - BAD!')
        continue
    else:
        #print(f"Working on: {h5parmFile}")
        h5 = h5parm.h5parm(h5parmFile, readonly=True)
        obsid = int(h5parmFile[2:].split('_-_')[0])
        solset = h5.getSolset('sol000')
        soltab = solset.getSoltab('amplitudeRes')
        
        vals = soltab.getValues(retAxesVals=False)
        weights = soltab.getValues(weight=True, retAxesVals=False)
        times = soltab.getAxisValues('time')
        if np.all(weights == 0):
            print(f'{h5parmFile} - All flagged - BAD!')
            badids.append(obsid)
            valstd = np.nan
            valdynrange = np.nan
        elif np.sum(weights == 0) > 2.*np.sum(weights != 0):
            print(f'{h5parmFile} - > 66% flagged - BAD!')
            badids.append(obsid)
            valstd = np.nan
            valdynrange = np.nan
        else:
            valstd = np.std(vals[weights!=0])
            valmin = np.min(vals[weights!=0])
            valmax = np.max(vals[weights!=0])
            valdynrange = valmax-valmin
    
            if (valstd > 0.15 and valdynrange > 1):
                print(f'{h5parmFile} - Std all residuals: {valstd} - Dyn range: {valdynrange} (elev: {elev}, sun dist: {sun_dist}, jup dist: {jup_dist}) - BAD!')
                badids.append(obsid)
            else:
                print(f'{h5parmFile} - Std all residuals: {valstd} - Dyn range: {valdynrange} (elev: {elev}, sun dist: {sun_dist}, jup dist: {jup_dist})')

   
    allsundist.append(sun_dist)
    alljupdist.append(jup_dist)
    allelev.append(elev)
    allvalstd.append(valstd)
    allvaldynrange.append(valdynrange)
    alltimes.append(times[0])
    allids.append(obsid)
    h5.close()

pickle.dump(badids, open(badobs_picklefile, 'wb'))
pickle.dump({'time':np.array(alltimes),'std':np.array(allvalstd),'dynrange':np.array(allvaldynrange),'id':np.array(allids),'elev':np.array(allelev),\
        'sundist':allsundist,'jupdist':alljupdist}, open('stats.pickle', 'wb'))
