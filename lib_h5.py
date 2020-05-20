import os, sys
import numpy as np
from losoto.h5parm import h5parm, Soltab
from LiLF.lib_log import logger

def repoint(h5parmFile, dirname):
    """
    rename the pointing direction of an h5parm from 'pointing' to 'dirname'
    """
    
    logger.info('%s: update dir name "pointing" -> "%s".' % (h5parmFile, dirname))
    dirname = '%s' % dirname

    # open h5parm
    h5 = h5parm(h5parmFile, readonly=False)
    ss = h5.getSolset('sol000')

    # rename each axes (must re-create the array as otherwise it truncates the dir name to the previous max string length)
    for tab in ss.getSoltabs():
        if 'dir' in tab.getAxesNames():
            tab.obj._v_file.remove_node('/'+tab.getAddress(), 'dir')
            tab.obj._v_file.create_array('/'+tab.getAddress(), 'dir', obj=[dirname.encode()])

    # rename directions table
    sourceTable = ss.obj.source
    direction = sourceTable[0][1]
    sourceTable[0] = (dirname,direction)

    # write h5parm
    sourceTable.close()
    h5.close()

