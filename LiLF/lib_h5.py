import os, sys
import numpy as np
import lsmtool
from astropy.coordinates import SkyCoord
from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes
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
    stokesI = False

    if 'pol' in st.getAxesNames():
        pol = st.getAxisValues('pol')
        if len(pol) == 1 and pol[0] == 'I':
            logger.warning('%s: Stokes I polarisation present in %s. Setting I -> XX,YY' % (h5parmFile, soltabname))
            stokesI = True
        else:
            h5.close()
            logger.warning('%s: polarisation axis already present in %s.' % (h5parmFile, soltabname))
            return
    
    # create values for new soltab
    typ = st.getType()

    vals = st.getValues(retAxesVals=False)
    weights = st.getValues(weight=True, retAxesVals=False)
    if stokesI: # extend I to XX,YY pol
        axesNames = st.getAxesNames()
        axesVals = []
        for axisName in st.getAxesNames():
            if axisName == 'pol':
                axesVals.append(['XX','YY'])
            else:
                axesVals.append(st.getAxisValues(axisName))
        vals = np.repeat(vals,2,-1)
        weights = np.repeat(weights,2,-1)
    else:
        axesNames = st.getAxesNames()+['pol']
        axesVals = [st.getAxisValues(axisName) for axisName in st.getAxesNames()]+[np.array(['XX','YY'])]
        vals = np.array([vals,vals])
        vals = np.moveaxis(vals, 0, -1)
        weights = np.array([weights,weights])
        weights = np.moveaxis(weights, 0, -1)

    # remove old soltab
    st.delete()

    # make new soltab
    soltabout = ss.makeSoltab(soltype = typ, soltabName = soltabname, axesNames=axesNames, \
                axesVals=axesVals, vals=vals, weights=weights)

    # write h5parm
    h5.close()


def reorder_axes(h5parmFile, order, soltabname, solsetname='sol000'):
    """
    reorder axis of a soltab
    """
    # open h5parm
    logger.info('%s: reorder axes: %s.' % (h5parmFile, soltabname))
    h5 = h5parm(h5parmFile, readonly=False)
    ss = h5.getSolset(solsetname)
    st = ss.getSoltab(soltabname)

    if st.getAxesNames() == order:
        logger.info('Axes are already in requested order.')
        h5.close()
        return None
    elif not set(order) == set(st.getAxesNames()):
        logger.error(f'Soltab axes {st.getAxesNames()} cannot be mapped to requested order {order}')
        h5.close()
        sys.exit()

    vals = st.getValues(retAxesVals=False)
    weights = st.getValues(weight=True, retAxesVals=False)
    vals = reorderAxes(vals, st.getAxesNames(), order)
    weights = reorderAxes(weights, st.getAxesNames(), order)

    # create values for new soltab
    typ = st.getType()
    axesVals = [st.getAxisValues(axisName) for axisName in order]

    # remove old soltab
    st.delete()

    # make new soltab
    ss.makeSoltab(soltype=typ, soltabName=soltabname, axesNames=order, \
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

def point_h5dirs_to_skymodel(h5, skymodel):
    """
    Change the [ra,dec] of h5parm directions
    Parameters
    ----------
    h5: input file
    skymodel: skymodel from which we take patch names and dirs
    """
    h5 = h5parm(h5, readonly=False)
    ss = h5.getSolset('sol000')
    sourceTable = ss.obj.source
    skymodel = lsmtool.load(skymodel)
    patches = skymodel.getPatchPositions()
    dirnames = patches.keys()

    # check which of the directions to rename is in solset
    for dirname in dirnames:
        direction = patches[dirname]
        # iterate sourcetable
        repointed = False
        for i, (dirname_ss, direction_ss) in enumerate(sourceTable[:]):
            if str(dirname_ss, 'utf-8') == f'[{dirname}]': # convert bytestring
                ra, dec = direction[0].rad, direction[1].rad
                if ra > np.pi:
                    ra -= 2*np.pi  # map to -pi, pi
                sourceTable[i] = (sourceTable[i][0], np.array([direction[0].rad, direction[1].rad]))
                logger.info(f'Repointing dir {dirname} {direction_ss} -> {[ra, dec]}')
                repointed = True

        if not repointed:
            raise ValueError(f'{dirname} not found in solset! Available directions: {",".join(ss.getSou().keys())}.')
    h5.close()

def get_colsest_dir(h5, dir):
    """

    Parameters
    ----------
    h5
    dir

    Returns
    -------
    d
    """
    h5 = h5parm(h5)
    directions = h5.getSolset('sol000').getSou()
    dir = SkyCoord(dir[0], dir[1], unit='deg')
    dir_seps = dict()

    for name, (ra, dec) in directions.items():
        dir_seps[name] = dir.separation(SkyCoord(np.degrees(ra), np.degrees(dec), unit='deg'))

    def key_of_min(d):
        return min(d, key=d.get)

    closest = key_of_min(dir_seps)
    print(closest, dir_seps[closest], dir_seps)
    h5.close()
    return closest
