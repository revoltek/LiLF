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
    return phases, data

def write_noflags_solutions(solspath, soltype="phase"):
    soltabname = f"{soltype}000"
    phases, data = open_sols(solspath, soltab=soltabname, constrain=False, apply_flags=False)
    new_solspath = f"{solspath[:-3]}-noflag.h5"

    import os
    try: os.remove(new_solspath)
    except: pass
    os.system(f"cp {solspath} {new_solspath}")
    
    sols = losoto.h5parm(new_solspath, readonly = False)
    solset = sols.getSolset(solset = "sol000")
    soltab = solset.getSoltab(soltab = soltabname)
    __, data = soltab.getValues()
    
    soltab.delete()

    weights = np.ones_like(phases)
    names = list(data.keys())
    vals = [data[key] for key in names]
    solset.makeSoltab(soltype=soltype, soltabName=soltabname, axesNames=names, axesVals=vals, vals=phases, weights=weights)
    sols.close()


if __name__ == "__main__":
    import sys
    write_noflags_solutions(sys.argv[1], soltype=sys.argv[2])
    #solspath = "/local/work/j.boxelaar/data/deepfields/cals/LOFAR_CAL-bkp/cal-fr.h5"
