import sys, os
import losoto.h5parm as losoto
import numpy as np

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
        #print("flagged data percentage:", np.sum(flags < 0.5)/len(flags))
    return phases, data

def add_dir_to_sols(path, new_path):
    # obscure way to avoid overwriting the original file if the new path is the same as the old path, since h5parm doesn't support in-place editing
    if new_path == path:
        final_path = path
        new_path = path + '_temp.h5'
    else:
        final_path = None
    
    os.system(f'cp {path} {new_path}')
    
    refsols = losoto.h5parm(path, readonly = True)
    refsolset = refsols.getSolset(solset = "sol000")
    refsoltab_names = refsolset.getSoltabNames()
    refsoltab = refsolset.getSoltab(soltab = refsoltab_names[0])
    __, refdata = refsoltab.getValues()
    refsols.close()
    
    sols = losoto.h5parm(new_path, readonly = False)
    solset = sols.getSolset(solset = "sol000")
    if not 'dir' in refdata:
        for soltab_name in refsoltab_names:
            soltab = solset.getSoltab(soltab = soltab_name)
            phases, data = soltab.getValues()
            weights, __ = soltab.getValues(weight=True)
            data.update({'dir': np.array(['direction'], dtype='U128')})
            soltab.delete()
             
            phases = phases.reshape(phases.shape + (1,))
            weights = weights.reshape(weights.shape + (1,))

            solset.makeSoltab(
                soltype=soltab_name[:-3], 
                soltabName=soltab_name, 
                axesNames=data.keys(), 
                axesVals=[data[key] for key in data.keys()], 
                vals=phases, 
                weights=weights
            )
    sols.close()
    
    if final_path is not None:
        os.system(f'cp {new_path} {final_path}')
        os.system(f'rm {new_path}')

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        new_path = sys.argv[1]
    else:
        new_path = sys.argv[2]
        
    add_dir_to_sols(sys.argv[1], new_path)

