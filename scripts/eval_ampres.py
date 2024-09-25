#!/usr/bin/env python3

import os, sys
import numpy as np
from losoto import h5parm

# get input
h5parmFile = 'cal-bp.h5'
soltabName = 'amplitudeSmooth'

h5 = h5parm.h5parm(h5parmFile, readonly=True)
solset = h5.getSolset('sol000')
soltab = solset.getSoltab(soltabName)

vals = soltab.getValues(retAxesVals=False)
weights = soltab.getValues(weight=True, retAxesVals=False)

sumall = np.std(vals[weights!=0])
print(f'Std all residuals: {sumall}')
h5.close()