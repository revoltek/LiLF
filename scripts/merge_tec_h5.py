#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "J. Boxelaar"

import numpy as np
import losoto.h5parm as losoto # type: ignore
import argparse

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('path_in', help='(str) h5 parm paths', type=str, nargs='+')
    parser.add_argument('--soltype', help='(str) phase or amplitude', default='tec', type=str)
    parser.add_argument('--path_out', help='out path', default=None)
    parser.add_argument('--mode', help='(str) add or subtract', default='add', type=str)
    return parser.parse_args()

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
    
def merge_h5_parms(solspaths: list[str], soltabname:str='tec000'):
    #
    # IMPORTANT: Assumes all h5 parms are made for the same MS.  (equal freq and time range and sizes)
    #
    base_tec, base_data = open_sols(solspaths[0], soltab=soltabname)
    base_tec = np.zeros_like(base_tec)
    base_weights = np.ones_like(base_tec)
    base_antaxis = list(base_data.keys()).index("ant")
    
    for path in solspaths:
        tec, data, weights = open_sols(path, soltab=soltabname, apply_flags=True)
        #assert base_data['time'][0] == data['time'][0], "Time ranges do not match"
        #assert base_data['freq'][0] == data['freq'][0], "Frequency ranges do not match"
        
        if args.mode == 'subtract':
            sign = -1
        elif args.mode == 'add':
            sign = 1
        else:
            raise ValueError("Mode must be either 'add' or 'subtract'")
        
        antaxis = list(data.keys()).index("ant")
        for i, ant in enumerate(base_data['ant']):
            if ant not in data['ant']:
                print(f"Antenna {ant} not found in {path}")
                continue
            
            if antaxis == 0:
                tec_vals = tec[data['ant'] == ant, :].squeeze()
                weights_vals = weights[data['ant'] == ant, :].squeeze()
            else:
                tec_vals = tec[:, data['ant'] == ant].squeeze()
                weights_vals = weights[:, data['ant'] == ant].squeeze()
            
            if base_antaxis == 0:
                base_tec[i, :] += sign * tec_vals
                base_weights[i,:] = np.logical_and(base_weights[i,:], weights_vals).astype(int)
            else:
                base_tec[:, i] += sign * tec_vals
                base_weights[:, i] = np.logical_and(base_weights[:,i], weights_vals).astype(int)
                   
    return base_tec, base_data, base_weights
        

def write_solutions(solspaths, soltype="phase", path_out=None):
    soltabname = f"{soltype}000"
    #phases, data, weights = open_sols(solspath, soltab=soltabname, constrain=False, apply_flags=True)
    if path_out is None:
        new_solspath = f"{solspaths[0][:-3]}-m.h5"
    else:
        new_solspath = path_out
 
    import os
    try: os.remove(new_solspath)
    except: pass
    os.system(f"cp {solspaths[0]} {new_solspath}")
    
    sols = losoto.h5parm(new_solspath, readonly = False)
    solset = sols.getSolset(solset = "sol000")
    solset.getSoltab(soltab = soltabname).delete()
    
    tec, data, weights = merge_h5_parms(solspaths, soltabname=soltabname)  

    #weights = np.ones_like(phases)
    names = list(data.keys())
    vals = [data[key] for key in names]
    solset.makeSoltab(soltype=soltype, soltabName=soltabname, axesNames=names, axesVals=vals, vals=tec, weights=weights)
    sols.close()


if __name__ == "__main__":
    args = arguments()
    write_solutions(args.path_in, soltype=args.soltype, path_out=args.path_out)
