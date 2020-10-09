#! /usr/bin/env python
"""
Script to combine two h5parms
"""
import argparse
from argparse import RawTextHelpFormatter
from factor.lib import miscellaneous as misc
from losoto.h5parm import h5parm
import scipy.interpolate as si
import sys
import os


def main(h5parm1, h5parm2, outh5parm, solset1='sol000', solset2='sol000', add_values=False, add_soltab='tec000'):
    """
    Combines two h5parms

    Parameters
    ----------
    h5parm1 : str
        Filenames of h5parm 1
    h5parm2 : str
        Filenames of h5parm 2
    outh5parm : str
        Filename of the output h5parm
    solset1 : str, optional
        Name of solset for h5parm1
    solset2 : str, optional
        Name of solset for h5parm2
    add_values : bool, optional
        If True, add values of the two h5parms, resampling on time if necessary (input
        h5parms must have the same axes)
    add_soltab : str, optional
        Name of soltab values to add
    """
    add_values = misc.string2bool(add_values)
    h1 = h5parm(h5parm1)
    h2 = h5parm(h5parm2)
    ss1 = h1.getSolset(solset=solset1)
    ss2 = h2.getSolset(solset=solset2)
    if os.path.exists(outh5parm):
        os.remove(outh5parm)
    ho = h5parm(outh5parm, readonly=False)

    sso = ho.makeSolset(solsetName = 'sol000', addTables=False)

    if add_values:
        # Figure out which one has finest time grid and interpolate other onto that grid,
        # then add
        st1 = ss1.getSoltab(add_soltab)
        st2 = ss2.getSoltab(add_soltab)
        axis_names = st1.getAxesNames()  # assume both have same axes
        time_ind = axis_names.index('time')
        if st1.val.shape[time_ind] > st2.val.shape[time_ind]:
            f = si.interp1d(st2.time, st2.val, axis=time_ind, kind='nearest', fill_value='extrapolate')
            vals = f(st1.time) + st1.val
            ss1.obj._f_copy_children(sso.obj, recursive=True, overwrite=True)
            sto = sso.getSoltab(add_soltab)
            sto.setValues(vals)
        else:
            f = si.interp1d(st1.time, st1.val, axis=time_ind, kind='nearest', fill_value='extrapolate')
            vals = f(st2.time) + st2.val
            ss2.obj._f_copy_children(sso.obj, recursive=True, overwrite=True)
            sto = sso.getSoltab(add_soltab)
            sto.setValues(vals)
    else:
        # Just copy over both solsets
        ss1.obj._f_copy_children(sso.obj, recursive=True, overwrite=True)
        ss2.obj._f_copy_children(sso.obj, recursive=True, overwrite=True)

    h1.close()
    h2.close()
    ho.close()


if __name__ == '__main__':
    descriptiontext = "Combine two h5parms.\n"

    parser = argparse.ArgumentParser(description=descriptiontext, formatter_class=RawTextHelpFormatter)
    parser.add_argument('h51', help='name of input h5 1')
    parser.add_argument('h52', help='name of input h5 2')
    parser.add_argument('outh5', help='name of the output h5')
    args = parser.parse_args()

    main(args.h51, args.h52, args.outh5)
