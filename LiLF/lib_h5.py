import os, sys
import numpy as np
from losoto.h5parm import h5parm
from LiLF.lib_log import logger

def repoint(h5parmFile, dirname, solsetname='sol000'):
    """
    rename the pointing direction of an h5parm from 'pointing' to 'dirname'
    """
    
    dirname = '%s' % dirname

    # open h5parm
    h5 = h5parm(h5parmFile, readonly=False)
    ss = h5.getSolset(solsetname)

    # rename each axes (must re-create the array as otherwise it truncates the dir name to the previous max string length)
    for tab in ss.getSoltabs():
        if 'dir' in tab.getAxesNames():
            tab.obj._v_file.remove_node('/'+tab.getAddress(), 'dir')
            tab.obj._v_file.create_array('/'+tab.getAddress(), 'dir', obj=[dirname.encode()])

    # rename directions table
    sourceTable = ss.obj.source
    direction = sourceTable[0][1]
    logger.info('%s: update dir name "%s" -> "%s".' % (h5parmFile,sourceTable[0][0], dirname))
    sourceTable[0] = (dirname,direction)

    # write h5parm
    sourceTable.close()
    h5.close()

def addpol(h5parmFile, soltabname, solsetname='sol000'):
    """
    add pol axes on a soltab
    """
    # open h5parm
    logger.info('%s: add pol axis to %s.' % (h5parmFile, soltabname))
    h5 = h5parm(h5parmFile, readonly=False)
    ss = h5.getSolset(solsetname)
    st = ss.getSoltab(soltabname)

    if 'pol' in st.getAxesNames():
        h5.close()
        logger.warning('%s: polarisation axis already present in %s.' % (h5parmFile, soltabname))
        return
    
    # create values for new soltab
    typ = st.getType()
    axesNames = st.getAxesNames()+['pol']
    axesVals = [st.getAxisValues(axisName) for axisName in st.getAxesNames()]+[np.array(['XX','YY'])]
    vals = st.getValues(retAxesVals = False)
    vals = np.array([vals,vals])
    vals = np.moveaxis(vals, 0, -1)
    weights = st.getValues(weight = True, retAxesVals = False)
    weights = np.array([weights,weights])
    weights = np.moveaxis(weights, 0, -1)

    # remove old soltab
    st.delete()

    # make new soltab
    soltabout = ss.makeSoltab(soltype = typ, soltabName = soltabname, axesNames=axesNames, \
                axesVals=axesVals, vals=vals, weights=weights)

    # write h5parm
    h5.close()


def adddir(h5parmFile, soltabname, solsetname='sol000', dirname='[pointing]'):
    """
    add dir axis on a soltab
    """
    # open h5parm
    logger.info('%s: add dir axis to %s.' % (h5parmFile, soltabname))
    h5 = h5parm(h5parmFile, readonly=False)
    ss = h5.getSolset(solsetname)
    st = ss.getSoltab(soltabname)

    if 'dir' in st.getAxesNames():
        h5.close()
        logger.warning('%s: direction axis already present in %s.' % (h5parmFile, soltabname))
        return

    # create values for new soltab
    typ = st.getType()
    axesNames = ['dir'] + st.getAxesNames()
    if dirname not in ss.getSou().keys():
        logger.error(f'Direction {dirname} is not in Solset!')
    axesVals = [np.array([dirname])] + [st.getAxisValues(axisName) for axisName in st.getAxesNames()]
    vals = st.getValues(retAxesVals=False)
    vals = np.array([vals])
    # vals = np.moveaxis(vals, 0, -1)
    weights = st.getValues(weight=True, retAxesVals=False)
    weights = np.array([weights])
    # weights = np.moveaxis(weights, 0, -1)
    # remove old soltab
    st.delete()

    # make new soltab
    soltabout = ss.makeSoltab(soltype=typ, soltabName=soltabname, axesNames=axesNames, \
                              axesVals=axesVals, vals=vals, weights=weights)

    # write h5parm
    h5.close()
