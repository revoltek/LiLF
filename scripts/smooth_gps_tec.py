#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "J. Boxelaar"

import numpy as np
import losoto.h5parm as losoto # type: ignore

def open_sols(path, soltab="phase000", constrain=False, apply_flags=False):
    sols = losoto.h5parm(path, readonly = True)
    solset = sols.getSolset(solset = "sol000")
    #print(solset.getSoltabNames())
    
    soltab = solset.getSoltab(soltab = soltab)  
    phases, data = soltab.getValues()
    if constrain:
        phase_shift = phases.copy()
        phase_shift[phases > np.pi] -= 2*np.pi
        phase_shift[phases < -np.pi] += 2*np.pi
        phases = phase_shift
    
    #print(data.keys())
    #print(phases.shape)
    if apply_flags:
        flags, __ = soltab.getValues(weight=True)
        phases[flags < 0.5] = np.nan
        print("flagged data percentage:", np.sum(flags < 0.5)/len(flags))
        return phases, data, flags
    else:
        return phases, data
    
def smooth_phases_in_time(phases, data, window=5):
    from scipy.interpolate import make_smoothing_spline
    
    smooth_phases = np.zeros_like(phases)
    for i, ant in enumerate(data["ant"]):
        #dtec = (phases[data['ant']==ant,:]-phases[data['ant']=="CS002LBA",:]).squeeze()
        dtec = phases[data['ant']==ant,:].squeeze()
        spl = make_smoothing_spline(data['time'], dtec, lam=1e8)
        smooth_phases[data['ant']==ant,:] = spl(data['time'])
        
        #import matplotlib.pyplot as plt
        #plt.plot(data['time'], dtec, label="dtec")
        #plt.plot(data['time'], spl(data['time']), label="smoothed dtec")
        #plt.legend()
        #plt.savefig(f"smoothed_dtec_{ant}.png")
        #plt.close()

    return smooth_phases

def write_solutions(solspath, soltype="phase"):
    soltabname = f"{soltype}000"
    phases, data, weights = open_sols(solspath, soltab=soltabname, constrain=False, apply_flags=True)
    new_solspath = f"{solspath[:-3]}-smooth.h5"
 
    import os
    try: os.remove(new_solspath)
    except: pass
    os.system(f"cp {solspath} {new_solspath}")
    
    sols = losoto.h5parm(new_solspath, readonly = False)
    solset = sols.getSolset(solset = "sol000")
    solset.getSoltab(soltab = soltabname).delete()
    
    smooth_phases = smooth_phases_in_time(phases, data)  

    #weights = np.ones_like(phases)
    names = list(data.keys())
    vals = [data[key] for key in names]
    solset.makeSoltab(soltype=soltype, soltabName=soltabname, axesNames=names, axesVals=vals, vals=smooth_phases, weights=weights)
    sols.close()


if __name__ == "__main__":
    import sys, os
    write_solutions(sys.argv[1], soltype=sys.argv[2])
    os.system(f"mv {sys.argv[1][:-3]}-smooth.h5 {sys.argv[1]}") # overwrite original file
