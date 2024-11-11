#!/usr/bin/env python3

"""
h5parm merger script for merging h5parms containing phase and/or amplitude solutions for calibrating LOFAR observations.

This script is a standalone script, such that it is easy to copy-paste it to other pipelines.

For examples to run this script from the command line, please have a look at:
https://github.com/jurjen93/lofar_helpers/blob/master/h5_merger_examples.md

It is a work in progress, so for any bug reports or suggestions for improvements:
please contact jurjendejong@strw.leidenuniv.nl (or find me on LOFAR slack pages)
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from casacore import tables as ct
from collections import OrderedDict
from glob import glob
from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes
from numpy import zeros, ones, round, unique, array_equal, append, where, isfinite, complex128, expand_dims, \
    pi, array, all, exp, angle, sort, sum, finfo, take, diff, equal, take, transpose, cumsum, insert, abs, asarray, \
    newaxis, argmin
import os
import re
from scipy.interpolate import interp1d
import sys
import tables
import warnings
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord
from astropy import units as u

warnings.filterwarnings('ignore')

__all__ = ['merge_h5', 'output_check', 'move_source_in_sourcetable', 'h5_check']

diagdiag_math = """
Diagonal times Diagonal
-------------------------------------------------

 /Axx  0 \     /Bxx  0 \     /Axx*Bxx   0   \ 
|         | X |         | = |                |
 \ 0  Ayy/     \ 0  Byy/     \   0   Ayy*Byy/

-------------------------------------------------
 """

diagfull_math = """
Diagonal times Fulljones
-------------------------------------------------

 /Axx  0 \     /Bxx Bxy\     /Axx*Bxx  Axx*Bxy\ 
|         | X |         | = |                  |
 \ 0  Ayy/     \Byx Byy/     \Ayy*Byx  Ayy*Byy/

-------------------------------------------------
 """

fulldiag_math = """
 Fulljones times Diagonal
-------------------------------------------------

 /Axx Axy \     /Bxx  0 \     /Axx*Bxx  Axy*Byy\ 
|          | X |         | = |                  |
 \Ayx  Ayy/     \ 0  Byy/     \Ayx*Bxx  Ayy*Byy/

-------------------------------------------------
"""

doublefulljones_math = """
 Fulljones times Fulljones
-------------------------------------------------

 /Axx Axy \     /Bxx  Bxy \     /Axx*Bxx + Axy*Byx  Axy*Byy + Axx*Bxy\ 
|          | X |           | = |                                      |
 \Ayx  Ayy/     \ Byx  Byy/     \Ayx*Bxx + Ayy*Byx  Ayy*Byy + Ayx*Bxy/

-------------------------------------------------
"""

lin2circ_math = """
Convert linear polarization to circular polarization
-----------------------------
RR = XX - iXY + iYX + YY
RL = XX + iXY + iYX - YY
LR = XX - iXY - iYX - YY
LL = XX + iXY - iYX + YY
-----------------------------
"""

circ2lin_math = """
Convert circular polarization to linear polarization
-----------------------------
XX = RR + RL + LR + LL
XY = iRR - iRL + iLR - iLL
YX = -iRR - iRL + iLR + iLL
YY = RR - RL - LR + LL
-----------------------------
"""

debug_message = 'Please contact jurjendejong@strw.leidenuniv.nl.'


def remove_numbers(inp):
    """
    Remove numbers from string (keep only letters)

    :param inp: string input
    """

    return "".join(re.findall("[a-zA-z]+", inp))


def make_utf8(inp):
    """
    Convert input to utf8 instead of bytes

    :param inp: string input
    """

    try:
        inp = inp.decode('utf8')
        return inp
    except (UnicodeDecodeError, AttributeError):
        return inp


def overwrite_table(T, solset, table, values, title=None):
    """
    Create table for given solset, opened with the package tables.
    Best to use for antenna or source table.

    :param T: Table opeend with tables
    :param solset: solution set of the table (f.ex. sol000)
    :param table: table name (f.ex. antenna or source)
    :param values: new values
    :param title: title of new table
    """

    try:
        T.root
    except:
        sys.exit(
            'ERROR: Create table failed. Given table is not opened with the package "tables" (https://pypi.org/project/tables/).')

    if 'sol' not in solset:
        print('WARNING: Usual input have sol*** as solset name.')

    ss = T.root._f_get_child(solset)
    ss._f_get_child(table)._f_remove()
    if table == 'source':
        values = array(values, dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
        title = 'Source names and directions'
    elif table == 'antenna':
        title = 'Antenna names and positions'
        values = array(values, dtype=[('name', 'S16'), ('position', '<f4', (3,))])
    else:
        try:  # check if numpy structure
            values.shape
        except:
            values = array(values)
    T.create_table(ss, table, values, title=title)

    return


def copy_antennas_from_MS_to_h5(MS, h5, solset):
    """
    Copy the antennas from an MS to an h5 file

    :param MS: measurement set
    :param h5: h5 file
    :param solset: solution set name
    """

    # Read MS antennas
    t = ct.table(MS + "::ANTENNA", ack=False)
    new_antlist = t.getcol('NAME')
    new_antpos = t.getcol('POSITION')
    antennas_ms = list(zip(new_antlist, new_antpos))
    t.close()

    # Open H5 antenna table
    T = tables.open_file(h5, 'r+')
    ss = T.root._f_get_child(solset)
    ants_h5 = [a.decode('utf8') if type(a) != str else a for a in
               T.root._f_get_child(solset)._f_get_child(list(ss._v_groups.keys())[0]).ant[:]]

    # Write antenna table
    if ants_h5 == new_antlist:
        overwrite_table(T, solset, 'antenna', antennas_ms, title=None)
    else:
        new_antennas = list(zip(ants_h5, [[0., 0., 0.]] * len(ants_h5)))
        for n, ant in enumerate(antennas_ms):
            ant_ms = make_utf8(ant[0])
            # add antenna from ms if not in H5
            if ant_ms in list(ants_h5):
                new_antennas[ants_h5.index(ant_ms)] = ant
        ss.antenna._f_remove()
        T.create_table(ss, 'antenna', array(new_antennas, dtype=[('name', 'S16'), ('position', '<f4', (3,))]),
                       title='Antenna names and positions')
    T.close()

    return


def find_closest_indices(arr1, arr2):
    """
    Index mapping between two arrays where each index in arr1 corresponds to the index of the closest value in arr2.

    Parameters:
        arr1 (np.array): The first array.
        arr2 (np.array): The second array, where we find the closest value to each element of arr1.

    Returns:
        np.array: An array of indices from arr2, corresponding to each element in arr1.
    """
    # Convert lists to NumPy arrays if not already
    arr1 = asarray(arr1)
    arr2 = asarray(arr2)

    # Calculate the absolute differences between each element of arr1 and all elements in arr2
    # The resulting matrix will have shape (len(arr1), len(arr2))
    diff_matrix = abs(arr1[:, newaxis] - arr2)

    # Find the index of the minimum value in arr2 for each element in arr1
    # np.argmin will return the indices of the closest values along axis 1 (across columns)
    closest_indices = argmin(diff_matrix, axis=1)

    return closest_indices


def take_numpy_axis(numpy_array, axis, idx):
    """
    Take specific axis from multidimensional array

    :param numpy_array: numpy array
    :param axis: axis
    :param idx: indices

    :return numpy array
    """

    if axis == 0:
        return numpy_array[idx, ...]
    elif axis == 1:
        return numpy_array[:, idx, ...]
    elif axis == 2:
        return numpy_array[:, :, idx, ...]
    elif axis == 3:
        return numpy_array[:, :, :, idx, ...]
    elif axis == 4:
        return numpy_array[:, :, :, :, idx, ...]
    else:
        sys.exit('ERROR: typical index number does not exceed 4, but you asked for ' + str(axis) +
                 '\nOr this script is outdated or you made a mistake.')


def has_integer(input):
    """
    Check if string has integer

    :param input: input string

    :return: return boolean (True/False) if integer in string
    """

    try:
        for s in str(input):
            if s.isdigit():
                return True
        return False
    except:  # dangerous but ok for now ;-)
        return False


def coordinate_distance(c1, c2):
    """
    Find distance between sources

    :param c1: first coordinate
    :param c2: second coordinate
    """
    if sys.version_info.major == 2:
        print("WARNING: --min_distance function does not work in Python 2")
        return 9999

    if abs(c1[0]) < 2 * pi and abs(c1[1]) < pi / 2 and abs(c2[0]) < 2 * pi and abs(c2[1]) < pi / 2:
        c1 = SkyCoord(c1[0], c1[1], unit='radian', frame='icrs')  # your coords
        c2 = SkyCoord(c2[0], c2[1], unit='radian', frame='icrs')
    else:
        c1 = SkyCoord(c1[0], c1[1], unit='degree', frame='icrs')  # your coords
        c2 = SkyCoord(c2[0], c2[1], unit='degree', frame='icrs')
    return c1.separation(c2).to(u.degree).value


class MergeH5:
    """Merge multiple h5 tables"""

    def __init__(self, h5_out, h5_tables=None, ms_files=None, h5_time_freq=None, convert_tec=True,
                 merge_all_in_one=False, solset='sol000', filtered_dir=None, no_antenna_crash=None,
                 freq_concat=None, time_concat=None):
        """
        :param h5_out: name of merged output h5 table
        :param h5_tables: h5 tables to merge, can be both list and string
        :param ms_files: read time and frequency from measurement set
        :param h5_time_freq: read time and frequency from h5
        :param convert_tec: convert TEC to phase or not
        :param merge_all_in_one: merge all in one direction
        :param solset: solset name
        :param filtered_dir: directions to filter (needs to be list with indices)
        :param no_antenna_crash: do not crash if antenna tables are not the same between h5s
        :param freq_concat: merging tables with different frequencies
        """

        # output name
        self.h5name_out = h5_out

        # for now this is standard sol000, might change in future version
        self.solset = solset

        # read in ms files
        if type(ms_files) == list:
            self.ms = ms_files
        elif type(ms_files) == str:
            self.ms = glob(ms_files)
        else:
            self.ms = []

        # read in h5 solution files
        if type(h5_tables) == list:
            self.h5_tables = h5_tables
        elif type(h5_tables) == str:
            self.h5_tables = glob(h5_tables)
        else:
            print('No h5 table given. Use all h5 tables in current folder.')
            self.h5_tables = tuple(glob('*.h5'))

        # get time and freq axis
        if type(h5_time_freq) == bool and h5_time_freq == True:
            if len(self.ms) > 0:
                print('Ignore MS for time and freq axis, as --h5_time_freq is given.')
            self.ax_time = array([])
            self.ax_freq = array([])
            for h5_name in self.h5_tables:
                h5 = tables.open_file(h5_name)
                for solset in h5.root._v_groups.keys():
                    ss = h5.root._f_get_child(solset)
                    for soltab in ss._v_groups.keys():
                        st = ss._f_get_child(soltab)
                        axes = make_utf8(st.val.attrs['AXES']).split(',')
                        if 'time' in axes:
                            time = st._f_get_child('time')[:]
                            self.ax_time = sort(unique(append(self.ax_time, time)))
                        else:
                            print('No time axes in ' + h5_name + '/' + solset + '/' + soltab)
                        if 'freq' in axes:
                            freq = st._f_get_child('freq')[:]
                            self.ax_freq = sort(unique(append(self.ax_freq, freq)))
                        else:
                            print('No freq axes in ' + h5_name + '/' + solset + '/' + soltab)
                h5.close()
        elif type(h5_time_freq) == str:
            if len(self.ms) > 0:
                print('Ignore MS for time and freq axis, as --h5_time_freq is given.')
            print('Take the time and freq from the following h5 solution file:\n' + h5_time_freq)
            T = tables.open_file(h5_time_freq)
            self.ax_time = T.root.sol000.phase000.time[:]
            self.ax_freq = T.root.sol000.phase000.freq[:]
            T.close()

        # use ms files for available information
        elif len(self.ms) > 0:
            print('Take the time and freq from the following measurement sets:\n' + '\n'.join(self.ms))
            self.ax_time = array([])
            self.ax_freq = array([])
            for m in self.ms:
                t = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + m + '::SPECTRAL_WINDOW')
                self.ax_freq = append(self.ax_freq, t.getcol('CHAN_FREQ')[0])
                t.close()

                t = ct.table(m)
                self.ax_time = append(self.ax_time, t.getcol('TIME'))
                t.close()
            self.ax_time = array(sorted(unique(self.ax_time)))
            self.ax_freq = array(sorted(unique(self.ax_freq)))

        # if no ms files, use the longest time and frequency resolution from h5 tables
        else:
            print(
                'No MS or h5 file given for time/freq axis.\nUse frequency and time axis by combining all input h5 tables.')
            self.ax_time = array([])
            self.ax_freq = array([])
            for h5_name in self.h5_tables:
                h5 = tables.open_file(h5_name)
                for solset in h5.root._v_groups.keys():
                    ss = h5.root._f_get_child(solset)
                    for soltab in ss._v_groups.keys():
                        st = ss._f_get_child(soltab)
                        axes = make_utf8(st.val.attrs['AXES']).split(',')
                        if 'time' in axes:
                            time = st._f_get_child('time')[:]
                            self.ax_time = sort(unique(append(self.ax_time, time)))
                        else:
                            print('No time axes in ' + h5_name + '/' + solset + '/' + soltab)
                        if 'freq' in axes:
                            freq = st._f_get_child('freq')[:]
                            self.ax_freq = sort(unique(append(self.ax_freq, freq)))
                        else:
                            print('No freq axes in ' + h5_name + '/' + solset + '/' + soltab)
                h5.close()

        # get polarization output axis and check number of error and tec tables in merge list
        self.polarizations, polarizations = [], []
        self.doublefulljones = False
        self.tecnum, self.errornum = 0, 0  # to average in case of multiple tables
        for n, h5_name in enumerate(self.h5_tables):
            h5 = tables.open_file(h5_name)
            if 'phase000' in h5.root.sol000._v_children.keys() \
                    and 'pol' in make_utf8(h5.root.sol000.phase000.val.attrs["AXES"]).split(','):
                polarizations = h5.root.sol000.phase000.pol[:]

                # having two fulljones solution files to merge --> we will use a matrix multiplication for this type of merge
                if len(polarizations) == len(self.polarizations) == 4:
                    self.doublefulljones = True
                    self.fulljones_phases = OrderedDict()
                    self.fulljones_amplitudes = OrderedDict()

            # take largest polarization list/array
            if len(polarizations) > len(self.polarizations):
                if type(self.polarizations) == list:
                    self.polarizations = polarizations
                else:
                    self.polarizations = polarizations.copy()

            if 'tec000' in h5.root.sol000._v_children.keys():
                self.tecnum += 1
            if 'error000' in h5.root.sol000._v_children.keys():
                self.errornum += 1

            h5.close()

        # check if fulljones
        if len(self.polarizations) == 4:
            self.fulljones = True
        elif len(self.polarizations) > 4 or len(self.polarizations) == 3:
            sys.exit('Output solutions cannot have ' + str(len(self.polarizations)) + ' solutions')
        else:
            self.fulljones = False

        print('Output polarization:\n' + str(self.polarizations))
        self.poldim = len(self.polarizations)

        # validation checks
        if len(self.ax_freq) == 0:
            sys.exit('ERROR: Cannot read frequency axis from input MS set or input H5.')
        if len(self.ax_time) == 0:
            sys.exit('ERROR: Cannot read time axis from input MS or input H5.')
        if not self.have_same_antennas and not no_antenna_crash:
            sys.exit('ERROR: Antenna tables are not the same')

        # convert tec to phase?
        self.convert_tec = convert_tec
        self.merge_all_in_one = merge_all_in_one
        if filtered_dir:
            self.filtered_dir = filtered_dir
        else:
            self.filtered_dir = None

        # possible solution axis in order used for our merging script
        self.solaxnames = ['pol', 'dir', 'ant', 'freq', 'time']

        # directions in an ordered dictionary
        self.directions = OrderedDict()
        if len(self.directions) > 1 and self.doublefulljones:
            sys.exit(
                "ERROR: Merging not compatitable with multiple directions and double fuljones merge")  # TODO: update

        self.freq_concat = freq_concat
        self.time_concat = time_concat

    @property
    def have_same_antennas(self):
        """
        Compare antenna tables with each other.
        These should be the same.

        :return: boolean if antennas are the same (True/False).
        """

        for h5_name1 in self.h5_tables:
            H_ref = tables.open_file(h5_name1)
            for solset1 in H_ref.root._v_groups.keys():
                ss1 = H_ref.root._f_get_child(solset1)
                antennas_ref = ss1.antenna[:]
                for soltab1 in ss1._v_groups.keys():
                    if (len(antennas_ref['name']) != len(ss1._f_get_child(soltab1).ant[:])) or \
                            (not all(antennas_ref['name'] == ss1._f_get_child(soltab1).ant[:])):
                        # print(type(ss1._f_get_child(soltab1).ant[:] ))
                        message = '\n'.join(['\nMismatch in antenna tables in ' + h5_name1,
                                             'Antennas from ' + '/'.join([solset1, 'antenna']),
                                             str(antennas_ref['name']),
                                             'Antennas from ' + '/'.join([solset1, soltab1, 'ant']),
                                             ','.join(ss1._f_get_child(soltab1).ant[:].astype(str))])

                        print(message)
                        H_ref.close()
                        return False
                    for soltab2 in ss1._v_groups.keys():
                        if (len(ss1._f_get_child(soltab1).ant[:]) !=
                            len(ss1._f_get_child(soltab2).ant[:])) or \
                                (not all(ss1._f_get_child(soltab1).ant[:] ==
                                         ss1._f_get_child(soltab2).ant[:])):
                            message = '\n'.join(['\nMismatch in antenna tables in ' + h5_name1,
                                                 'Antennas from ' + '/'.join([solset1, soltab1, 'ant']),
                                                 ss1._f_get_child(soltab1).ant[:],
                                                 'Antennas from ' + '/'.join([solset1, soltab2, 'ant']),
                                                 ss1._f_get_child(soltab2).ant[:]])
                            print(message)
                            H_ref.close()
                            return False
                for h5_name2 in self.h5_tables:
                    H = tables.open_file(h5_name2)
                    for solset2 in H.root._v_groups.keys():
                        ss2 = H.root._f_get_child(solset2)
                        antennas = ss2.antenna[:]
                        if (len(antennas_ref['name']) != len(antennas['name'])) \
                                or (not all(antennas_ref['name'] == antennas['name'])):
                            message = '\n'.join(
                                ['\nMismatch between antenna tables from ' + h5_name1 + ' and ' + h5_name2,
                                 'Antennas from ' + h5_name1 + ' and ',
                                 str(antennas_ref['name']),
                                 'Antennas from ' + h5_name2 + ':',
                                 str(antennas['name'])])
                            print(message)
                            H.close()
                            H_ref.close()
                            return False
                    H.close()
            H_ref.close()

        return True

    def concat(self, values, soltab, axes, ax_name):
        """
        Concat instead of interpolation

        :param values: values (phase or amplitude)
        :param soltab: amplitude or phase soltab
        :param axes: h5 axes (to be concatted)
        :param ax_name: over time or freq

        :return: concattenated
        """

        if ax_name == 'freq':
            axes_new = self.ax_freq
        elif ax_name == 'time':
            axes_new = self.ax_time
        else:
            sys.exit("ERROR: Should not arrive here")

        shape = list(values.shape)
        shape[self.axes_current.index(ax_name)] = len(axes_new)
        values_tmp = zeros(shape)
        if 'pol' in self.axes_current and 'amplitude' in soltab:
            if self.axes_current.index('pol') == 0:
                values_tmp[-1, ...] = 1
                values_tmp[0, ...] = 1
            elif self.axes_current.index('pol') == 1:
                values_tmp[:, -1, ...] = 1
                values_tmp[:, 0, ...] = 1
            elif self.axes_current.index('pol') == 2:
                values_tmp[:, :, -1, ...] = 1
                values_tmp[:, :, 0, ...] = 1
            elif self.axes_current.index('pol') == 3:
                values_tmp[:, :, :, -1, ...] = 1
                values_tmp[:, :, :, 0, ...] = 1
            elif self.axes_current.index('pol') == 4:
                values_tmp[:, :, :, :, -1, ...] = 1
                values_tmp[:, :, :, :, 0, ...] = 1

        idx = find_closest_indices(axes, axes_new)

        if self.axes_current.index(ax_name) == 0:
            values_tmp[idx, ...] = values[:]

        elif self.axes_current.index(ax_name) == 1:
            values_tmp[:, idx, ...] = values[:]

        elif self.axes_current.index(ax_name) == 2:
            values_tmp[:, :, idx, ...] = values[:]

        elif self.axes_current.index(ax_name) == 3:
            values_tmp[:, :, :, idx, ...] = values[:]

        elif self.axes_current.index(ax_name) == 4:
            values_tmp[:, :, :, :, idx, ...] = values[:]

        return values_tmp

    def _unpack_h5(self, st, solset, soltab):
        """
        Unpack, check, and reorder the values from the h5 table to merge.

        :param st: solution table
        :param solset: solset name
        :param soltab: soltab name

        :return: values, time axis, frequency axis
        """

        # Check if there is a polarization axes
        if 'pol' in st.getAxesNames():
            print("polarization is in {solset}/{soltab}".format(solset=solset, soltab=soltab))
        else:
            print("polarization is not in {solset}/{soltab}".format(solset=solset, soltab=soltab))

        # Get time axes
        time_axes = st.getAxisValues('time')

        # Get freq axes
        if 'freq' in st.getAxesNames():
            freq_axes = st.getAxisValues('freq')
        else:
            freq_axes = self.ax_freq

        # Do checks
        if (self.ax_time[0] > time_axes[-1] or time_axes[0] > self.ax_time[-1]) and (
                len(time_axes) > 1 and len(self.ax_time) > 1):
            sys.exit("ERROR: Time axes of h5 and MS are not overlapping.\n"
                     "SUGGESTION: add --h5_time_freq=true if you want to use the input h5 files to construct the time and freq axis.")
        if (self.ax_freq[0] > freq_axes[-1] or freq_axes[0] > self.ax_freq[-1]) and (
                len(freq_axes) > 1 and len(self.ax_freq) > 1) and not self.freq_concat:
            sys.exit("ERROR: Frequency axes of h5 and MS are not overlapping.\n"
                     "SUGGESTION: add --h5_time_freq=true if you want to use the input h5 files to construct the time and freq axis.")
        if float(soltab[-3:]) > 0:
            sys.exit("ERROR: {soltab} does not end on 000".format(soltab=soltab))
        for av in self.axes_final:
            if av in st.getAxesNames() and st.getAxisLen(av) == 0:
                print("No {av} in {solset}/{soltab}".format(av=av, solset=solset, soltab=soltab))

        if len(st.getAxesNames()) != len(st.getValues()[0].shape):
            sys.exit('ERROR: Axes ({axlen}) and Value dimensions ({vallen}) are not equal'.format(
                axlen=len(st.getAxesNames()), vallen=len(st.getValues()[0].shape)))

        print('Value shape before --> {values}'.format(values=st.getValues()[0].shape))

        # Reorder and add dir
        axes_current = [an for an in self.solaxnames if an in st.getAxesNames()]
        if 'dir' in st.getAxesNames():
            values = reorderAxes(st.getValues()[0], st.getAxesNames(), axes_current)
        else:
            print('Add dir axis.')
            origin_values = st.getValues()[0]
            if 'dir' not in axes_current:
                if 'pol' not in axes_current:
                    axes_current.insert(0, 'dir')
                else:
                    axes_current.insert(1, 'dir')
            values = reorderAxes(origin_values.reshape(origin_values.shape + (1,)), st.getAxesNames() + ['dir'],
                                 axes_current)

        # current axes for reordering of axes
        self.axes_current = [an for an in self.solaxnames if an in st.getAxesNames()]

        if 'dir' not in self.axes_current:
            if 'pol' in self.axes_current:
                self.axes_current.insert(1, 'dir')
            else:
                self.axes_current.insert(0, 'dir')

        # remove invalid values
        values = self.remove_invalid_values(st.getType(), values, self.axes_current)

        if remove_numbers(st.getType()) == 'tec':
            # add frequencies
            if 'freq' not in st.getAxesNames():
                ax = self.axes_final.index('freq') - len(self.axes_final)
                values = expand_dims(values, axis=ax)
                self.axes_current.insert(-1, 'freq')
                valuestmp = values
                for _ in range(len(self.ax_freq) - 1):
                    values = append(values, valuestmp, axis=-2)

            if self.convert_tec:
                # convert tec to phase
                shape = [1 for _ in range(values.ndim)]
                shape[self.axes_current.index('freq')] = -1
                values = self.tecphase_conver(values, self.ax_freq.reshape(shape))
                soltab = 'phase'

        # expand pol dimensions
        if self.fulljones:
            if 'pol' not in self.axes_current:
                values = self._expand_poldim(values, 4, remove_numbers(soltab), False)
                self.axes_current.insert(0, 'pol')
            elif st.getAxisLen('pol') != 4:
                values = self._expand_poldim(values, 4, remove_numbers(soltab), True)
        elif 'pol' in self.axes_current and 'pol' in self.axes_final:
            if st.getAxisLen('pol') > self.phases.shape[0]:
                self.phases = self._expand_poldim(self.phases, st.getAxisLen('pol'), remove_numbers(soltab), True)
            elif len(self.polarizations) > st.getAxisLen('pol'):
                values = self._expand_poldim(values, len(self.polarizations), remove_numbers(soltab), True)
        elif 'pol' not in self.axes_current and 'pol' in self.axes_final:
            values = self._expand_poldim(values, len(self.polarizations), remove_numbers(soltab), False)
            self.axes_current.insert(0, 'pol')
        elif 'pol' in self.axes_current and 'pol' not in self.axes_final:
            self.phases = self._expand_poldim(self.phases, st.getAxisLen('pol'), remove_numbers(soltab), False)
            self.axes_final.insert(0, 'pol')
        elif 'pol' not in self.axes_current and 'pol' not in self.axes_final and len(
                self.polarizations) > 0:
            if len(self.phases.shape) != 5:
                self.phases = self._expand_poldim(self.phases, len(self.polarizations), remove_numbers(soltab), False)
                if 'pol' not in self.axes_final:
                    self.axes_final.insert(0, 'pol')
            values = self._expand_poldim(values, len(self.polarizations), remove_numbers(soltab), False)
            self.axes_final.insert(0, 'pol')

        # time interpolation
        if self.time_concat:
            values = self.concat(values, st.getType(), time_axes, 'time')
        else:
            values = self._interp_along_axis(values, time_axes, self.ax_time,
                                             self.axes_current.index('time'))

        # freq interpolation
        if self.freq_concat:
            values = self.concat(values, st.getType(), freq_axes, 'freq')
        elif remove_numbers(st.getType()) != 'tec' and remove_numbers(st.getType()) != 'error':
            values = self._interp_along_axis(values, freq_axes, self.ax_freq,
                                             self.axes_current.index('freq'))

        return values

    def _sort_soltabs(self, soltabs):
        """
        Sort solution tables.
        This is important to run the steps and add directions according to our algorithm.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dont touch this part if you dont have to. ;-)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        :param soltabs: solutions tables

        :return: sorted phase, tec, amplitude, rotation, error tables
        """

        soltabs = set(soltabs)
        if self.convert_tec:
            tp_phasetec = [li for li in soltabs if 'tec' in li or 'phase' in li]
            tp_amplitude = [li for li in soltabs if 'amplitude' in li]
            tp_rotation = [li for li in soltabs if 'rotation' in li]
            tp_error = [li for li in soltabs if 'error' in li]
            return [sorted(tp_amplitude, key=lambda x: float(x[-3:])),
                    sorted(tp_rotation, key=lambda x: float(x[-3:])),
                    sorted(sorted(tp_phasetec), key=lambda x: float(x[-3:])),
                    sorted(tp_error, key=lambda x: float(x[-3:]))]
        else:
            tp_phase = [li for li in soltabs if 'phase' in li]
            tp_tec = [li for li in soltabs if 'tec' in li]
            tp_amplitude = [li for li in soltabs if 'amplitude' in li]
            tp_rotation = [li for li in soltabs if 'rotation' in li]
            tp_error = [li for li in soltabs if 'error' in li]
            return [sorted(tp_phase, key=lambda x: float(x[-3:])),
                    sorted(tp_tec, key=lambda x: float(x[-3:])),
                    sorted(tp_amplitude, key=lambda x: float(x[-3:])),
                    sorted(tp_rotation, key=lambda x: float(x[-3:])),
                    sorted(tp_error, key=lambda x: float(x[-3:]))]

    def get_allkeys(self):
        """
        Get all solution sets, solutions tables, and ax names in lists.
        This returns an order that is optimized for this code.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dont touch this part if you dont have to. ;-)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        """

        self.all_soltabs, self.all_solsets, self.all_axes, self.ant = [], [], [], []

        for h5_name in self.h5_tables:
            h5 = h5parm(h5_name)
            for solset in h5.getSolsetNames():
                self.all_solsets += [solset]
                ss = h5.getSolset(solset)
                for n, soltab in enumerate(ss.getSoltabNames()):
                    self.all_soltabs += [soltab]
                    st = ss.getSoltab(soltab)
                    self.all_axes += ['/'.join([solset, soltab, an]) for an in st.getAxesNames()]
                    if n == 0:
                        self.ant = st.getAxisValues('ant')  # check if same for all h5
                    elif list(self.ant) != list(st.getAxisValues('ant')):
                        sys.exit('ERROR: antennas not the same.')
            h5.close()
        self.all_soltabs = self._sort_soltabs(self.all_soltabs)
        self.all_solsets = set(self.all_solsets)
        self.all_axes = set(self.all_axes)
        return self

    def _make_template_values(self, soltab):
        """
        Make template numpy array to fill up during merging.

        :param soltab: solution table name
        :param st: solution table itself
        """

        num_dir = max(len(self.directions), 1)

        if 'amplitude' in soltab and len(self.polarizations) > 0:
            self.amplitudes = ones(
                (max(len(self.polarizations), 1), num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))
            if len(self.polarizations) == 4:
                self.amplitudes[1, ...] = 0.
                self.amplitudes[2, ...] = 0.
        else:
            self.amplitudes = ones((num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))

        if 'phase' in soltab and len(self.polarizations) > 0:
            self.phases = zeros(
                (max(len(self.polarizations), 1), num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))
        elif 'rotation' in soltab:
            self.phases = zeros((num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))
        else:
            self.phases = zeros((num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))

        if 'tec' in soltab and not self.convert_tec and len(self.polarizations) > 0:
            self.tec = zeros(
                (max(len(self.polarizations), 1), num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))
        elif 'tec' in soltab and not self.convert_tec:
            self.tec = zeros((num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))

        if 'error' in soltab and len(self.polarizations) > 0:
            self.error = zeros(
                (max(len(self.polarizations), 1), num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))
        elif 'error' in soltab:
            self.error = zeros((num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))

        if self.doublefulljones:
            self.fullgains = ones((4, num_dir, len(self.ant), len(self.ax_freq), len(self.ax_time)))

        self.n = len(self.directions)  # direction number reset

        return self

    @staticmethod
    def tecphase_conver(tec, freqs):
        """
        convert tec to phase --> See Equation 1 from Sweijen et al. 2022

        :param tec: TEC
        :param freqs: frequencies

        :return: tec phase converted values
        """

        return -8.44797245e9 * tec / freqs

    @staticmethod
    def _interp_along_axis(x, interp_from, interp_to, axis, fill_value='extrapolate'):
        """
        Interpolate along axis.

        :param x: frequency or time axis. Must be equal to the length of interp_from.
        :param interp_from: interpolate from this axis.
        :param interp_to: interpolate to this axis
        :param axis: interpolation axis

        :return: the interpolated result
        """

        # axis with length 1
        if len(interp_from) == 1:
            new_vals = x
            for _ in range(len(interp_to) - 1):
                new_vals = append(new_vals, x, axis=axis)
        else:
            interp_vals = interp1d(interp_from, x, axis=axis, kind='nearest', fill_value=fill_value, bounds_error=False)
            new_vals = interp_vals(interp_to)
        return new_vals

    def get_model_h5(self, solset, soltab):
        """
        Get model (clean) h5 table.

        :param solset: solution set name (sol000)
        :param soltab: solution table name
        """

        if '000' in soltab:

            # make template
            for h5_name in self.h5_tables:

                h5_to_merge = h5parm(h5_name)

                if solset not in h5_to_merge.getSolsetNames():
                    h5_to_merge.close()
                    sys.exit('ERROR ' + solset + ' does not exist in ' + h5_name)
                else:
                    ss = h5_to_merge.getSolset(solset)
                if soltab not in ss.getSoltabNames():
                    h5_to_merge.close()
                    continue  # use other h5 table with soltab
                else:
                    st = ss.getSoltab(soltab)

                if not self.convert_tec or (self.convert_tec and 'tec' not in soltab):
                    self._make_template_values(soltab)
                    self.axes_final = [an for an in self.solaxnames if an in st.getAxesNames()]
                    if len(self.polarizations) > 1 and 'pol' not in st.getAxesNames():
                        self.axes_final.insert(0, 'pol')
                elif 'tec' in soltab and self.convert_tec:
                    for st_group in self.all_soltabs:
                        if soltab in st_group and ('phase000' not in st_group and 'phase{n}'.format(n=soltab[-3:])):
                            self._make_template_values(soltab)
                            self.axes_final = [an for an in self.solaxnames if an in st.getAxesNames()]
                            if len(self.polarizations) > 1 and 'pol' not in st.getAxesNames():
                                self.axes_final = ['pol'] + self.axes_final

                # add dir if missing
                if 'dir' not in self.axes_final:
                    if 'pol' in self.axes_final:
                        self.axes_final.insert(1, 'dir')
                    else:
                        self.axes_final.insert(0, 'dir')

                h5_to_merge.close()
                break

        return self

    @staticmethod
    def get_number_of_directions(st):
        """
        Get number of directions in solution table

        :param st: solution table

        :return: number of directions
        """

        if 'dir' in st.getAxesNames():
            dir_index = st.getAxesNames().index('dir')
            return st.getValues()[0].shape[dir_index]
        else:
            return 1

    def add_direction(self, source):
        """
        Add direction to dictionary

        :param source: source direction
        """

        self.directions.update(source)
        self.directions = OrderedDict(sorted(self.directions.items()))

        return self

    @staticmethod
    def remove_invalid_values(soltab, values, axlist):
        """
        Correct invalid values in tables (infinite, nan).

        Only finite values allowed.

        :param soltab: solution table name
        :param values: numpy array values
        :param axlist: list with axes

        :return: new values
        """

        if 'phase' in soltab or 'tec' in soltab:
            values[~isfinite(values)] = 0.
        elif 'amplitude' in soltab:
            if 'pol' in axlist:
                if axlist.index('pol') == 0:
                    values[(~isfinite(values))] = 0.
                    values[0, ...][(~isfinite(values[0, ...])) | (values[0, ...] == 0.)] = 1.
                    values[-1, ...][(~isfinite(values[-1, ...])) | (values[-1, ...] == 0.)] = 1.
                elif axlist.index('pol') - len(axlist) == -1:
                    values[(~isfinite(values))] = 0.
                    values[..., 0][(~isfinite(values[..., 0])) | (values[..., 0] == 0.)] = 1.
                    values[..., -1][(~isfinite(values[..., -1])) | (values[..., -1] == 0.)] = 1.
                else:
                    sys.exit('ERROR: polarization found at unexpected axis ' + str(axlist))
            else:
                values[(~isfinite(values)) | (values == 0.)] = 1.

        return values

    @staticmethod
    def _expand_poldim(values, dim_pol, type, haspol):
        """
        Add extra polarization dimensions

        :param values: values which need to get a polarization
        :param dim_pol: number of dimensions
        :param type: phase or amplitude
        :param haspol: has polarization in input values

        :return: input values with extra polarization axis
        """

        if dim_pol == 3 or dim_pol > 4:
            sys.exit('ERROR: Invalid number of polarizations.'
                     '\nOnly 1, 2, or 4 are allowed. '
                     'This corresponds to XX, (XY, YX,) YY or RR, (RL, LR,) LL.')

        if dim_pol == 4:
            print("Convert matrix to fulljones for merge")

        if values.ndim < 5 and not haspol:
            if type == 'amplitude':
                values_new = ones((dim_pol,) + values.shape)
                if dim_pol == 4:
                    values_new[1, ...] = 0
                    values_new[2, ...] = 0
            elif type in ['phase', 'error', 'tec']:
                values_new = zeros((dim_pol,) + values.shape)
            else:
                sys.exit('ERROR: Only type in [amplitude, phase] allowed [%s requested].' % str(type))
            for i in range(dim_pol):
                if dim_pol == 4:
                    if i in [0, 3]:
                        values_new[i, ...] = values
                else:
                    values_new[i, ...] = values

        elif values.shape[0] in [1, 2] and dim_pol in [2, 4] and haspol:
            if type == 'amplitude':
                values_new = ones((dim_pol,) + values.shape[1:])
                if dim_pol == 4:
                    values_new[1, ...] = 0
                    values_new[2, ...] = 0
            elif type in ['phase', 'error', 'tec']:
                values_new = zeros((dim_pol,) + values.shape[1:])
            else:
                sys.exit('ERROR: Only type in [amplitude, phase] allowed [%s requested].' % str(type))
            values_new[0, ...] = values[0, ...]
            if values.shape[0] == 2:
                values_new[-1, ...] = values[1, ...]
            elif values.shape[0] == 1:
                values_new[-1, ...] = values[0, ...]
        else:
            print('WARNING: No pol ax dimension changed.')
            return values

        return values_new

    def merge_tables(self, solset=None, soltab=None, min_distance=0.):
        """
        Merge solution files

        :param solset: solution set name
        :param soltab: solution table name
        """

        # flag necessary to track merge order (important because of matrix commutations)
        fulljones_done = False

        # loop over all tables
        for h5num, h5_name in enumerate(self.h5_tables):

            h5 = h5parm(h5_name)
            if solset not in h5.getSolsetNames():
                h5.close()
                continue

            ss = h5.getSolset(solset)

            if soltab not in ss.getSoltabNames():
                h5.close()
                continue

            st = ss.getSoltab(soltab)

            print('Solution table from {table}'.format(table=h5_name.split('/')[-1]))
            num_dirs = self.get_number_of_directions(st)  # number of directions
            print('This table has {dircount} direction(s)'.format(dircount=num_dirs))

            # get values from h5 file
            table_values = self._unpack_h5(st, solset, soltab)

            for dir_idx in range(num_dirs):  # loop over all directions in h5

                if self.filtered_dir != None:
                    if len(self.filtered_dir) > 0 and dir_idx not in self.filtered_dir:
                        continue

                # different axes if pol in output
                if 'pol' in self.axes_final:
                    # print(dir_idx, ss.getSou())
                    values = table_values[:, dir_idx, ...]
                    # if dir_idx == 0:
                    #     values *= 0.0
                else:
                    values = table_values[dir_idx, ...]

                # get source coordinates
                dirs = ss.getSou()
                if 'dir' in list(dirs.keys())[0].lower() and list(dirs.keys())[0][-1].isnumeric():
                    dirs = OrderedDict(sorted(dirs.items()))
                elif len(dirs) > 1 and (
                        sys.version_info.major == 2 or (sys.version_info.major == 3 and sys.version_info.minor < 6)):
                    print(
                        'WARNING: Order of source directions from h5 table might not be ordered. This is an old Python issue.'
                        '\nSuggest to switch to Python 3.6 or higher')
                # print('diridx',dir_idx,list(dirs.keys())[dir_idx])

                # coordinate list
                source_coords = dirs[list(dirs.keys())[dir_idx]]
                # print(dirs)
                # print(self.directions)
                print(22, source_coords, [list(sv) for sv in self.directions.values()])

                if self.merge_all_in_one and self.n == 1:
                    idx = 0
                    print('Merging direction {:f},{:f} with previous direction'.format(*source_coords))
                    if abs(self.directions['Dir00'][0]) > 0 and abs(self.directions['Dir00'][1]) > 0:
                        self.add_direction({'Dir00': source_coords})  # 0.0 coordinate bug
                        print('Adding new direction {:f},{:f}'.format(*source_coords))
                elif any([array_equal(source_coords, list(sv)) for sv in self.directions.values()]):
                    # Direction already exists, add to the existing solutions.
                    print('Direction {:f},{:f} already exists. Adding to this direction.'.format(*source_coords))
                    # Matching on 5 decimals rounding
                    idx = list([[round(l[0], 5), round(l[1], 5)] for l in self.directions.values()]). \
                        index([round(source_coords[0], 5), round(source_coords[1], 5)])
                elif any([coordinate_distance(source_coords, list(sv)) <= min_distance for sv in
                          self.directions.values()]):
                    # Direction closer than minimal allowed distance.
                    md = 999
                    for i, d in enumerate(self.directions.values()):
                        if coordinate_distance(source_coords, list(d)) < md:
                            md = coordinate_distance(source_coords, list(d))
                            idx = i
                    print('Direction {:f},{:f}'.format(*source_coords) + ' is closer than minimal distance of '
                          + str(min_distance) + ' to ' + str(list(self.directions.values())[idx])
                          .replace('[', '').replace(']', '').replace('  ', ','))

                else:  # new direction
                    if abs(source_coords[0]) > 0 and abs(source_coords[1]) > 0:
                        # print(f'Adding new direction Dir{self.n:02d} {source_coords[0]},{source_coords[1]}')
                        print(f'Create from existing direction {list(dirs.keys())[dir_idx]} {source_coords[0]},{source_coords[1]}')
                    idx = self.n
                    self.add_direction({list(dirs.keys())[dir_idx]: source_coords})
                    if not self.merge_all_in_one:
                        self.n += 1
                    if self.n > 1:  # for self.n==1 --> dont have to do anything
                        if (st.getType() in ['tec', 'phase', 'rotation'] and self.convert_tec) \
                                or (st.getType() in ['phase', 'rotation'] and not self.convert_tec):
                            shape = list(self.phases.shape)
                            dir_index = self.phases.ndim - 4
                            # print(dir_index, shape)
                            if dir_index < 0:
                                sys.exit('ERROR: Missing dir axes.')
                            if self.n > shape[dir_index]:
                                shape[dir_index] = 1
                                self.phases = append(self.phases, zeros(shape),
                                                     axis=dir_index)  # add clean phase to merge with
                        elif st.getType() == 'amplitude':
                            shape = list(self.amplitudes.shape)
                            dir_index = self.amplitudes.ndim - 4
                            if dir_index < 0:
                                sys.exit('ERROR: Missing dir axes.')
                            if self.n > shape[dir_index]:
                                shape[dir_index] = 1
                                self.amplitudes = append(self.amplitudes, ones(shape),
                                                         axis=dir_index)  # add clean gain to merge with
                        elif st.getType() == 'tec' and not self.convert_tec:
                            shape = list(self.tec.shape)
                            dir_index = self.tec.ndim - 4
                            if dir_index < 0:
                                sys.exit('ERROR: Missing dir axes.')
                            if self.n > shape[dir_index]:
                                shape[dir_index] = 1
                                self.tec = append(self.tec, zeros(shape),
                                                  axis=dir_index)  # add clean phase to merge with
                        elif st.getType() == 'error':
                            dir_index = self.error.ndim - 4
                            if dir_index < 0:
                                sys.exit('ERROR: Missing dir axes.')
                            if self.n > shape[dir_index]:
                                shape[dir_index] = 1
                                self.error = append(self.error, zeros(shape),
                                                    axis=dir_index)  # add clean phase to merge with

                # Add the solution table axis --> tec, phase, rotation, amplitude, or error
                if (
                        st.getType() == 'tec' and self.convert_tec) or st.getType() == 'phase' or st.getType() == 'rotation':

                    # add values
                    if self.fulljones and not self.doublefulljones:

                        if 'pol' in st.getAxesNames() and not fulljones_done:
                            if st.getAxisLen('pol') == 4:
                                self.phases[1, idx, ...] = self.phases[0, idx, ...] + values[
                                    1, ...]  # Axx * Bxy
                                self.phases[2, idx, ...] = self.phases[-1, idx, ...] + values[
                                    2, ...]  # Ayy * Byx
                                fulljones_done = True
                                # print(diagfull_math)

                        elif fulljones_done:
                            self.phases[1, idx, ...] += values[-1, ...]  # Axy * Byy
                            self.phases[2, idx, ...] += values[0, ...]  # Ayx * Bxx
                            # print(fulldiag_math)

                        self.phases[0, idx, ...] += values[0, ...]  # Axx * Bxx
                        self.phases[-1, idx, ...] += values[-1, ...]  # Ayy * Byy

                    elif not self.doublefulljones:
                        if 'pol' in self.axes_current:
                            import numpy as np
                            # print('putting at', idx, np.std(values[0]), values.shape)
                            self.phases[:, idx, ...] += values[:, ...]
                        else:
                            self.phases[idx, ...] += values[...]
                        # print(diagdiag_math)

                    elif self.doublefulljones:  # save table for in later double fulljones merge
                        self.fulljones_phases.update({h5_name: values})

                    else:
                        sys.exit("ERROR: Broken code. Should not be here.")

                elif st.getType() == 'tec' and not self.convert_tec:

                    # add values
                    if len(self.polarizations) > 0:
                        # average tec
                        self.tec[:, idx, ...] += values[:, ...] / self.tecnum
                    else:
                        # average tec
                        self.tec[idx, ...] += values[...] / self.tecnum

                elif st.getType() == 'error':

                    # add values
                    if len(self.polarizations) > 0:
                        # average error
                        self.error[:, idx, ...] += values[:, ...] / self.errornum
                    else:
                        # average error
                        self.error[idx, ...] += values[...] / self.errornum


                elif st.getType() == 'amplitude':

                    # add values
                    if self.fulljones and not self.doublefulljones:

                        if 'pol' in st.getAxesNames() and not fulljones_done:
                            if st.getAxisLen('pol') == 4:
                                self.amplitudes[1, idx, ...] = self.amplitudes[0, idx, ...] * values[
                                    1, ...]  # Axx * Bxy
                                self.amplitudes[2, idx, ...] = self.amplitudes[-1, idx, ...] * values[
                                    2, ...]  # Ayy * Byx
                                fulljones_done = True
                                # print(diagfull_math)

                        elif fulljones_done:
                            self.amplitudes[1, idx, ...] *= values[-1, ...]  # Axy * Byy
                            self.amplitudes[2, idx, ...] *= values[0, ...]  # Ayx * Bxx
                            # print(fulldiag_math)

                        self.amplitudes[0, idx, ...] *= values[0, ...]  # Axx * Bxx
                        self.amplitudes[-1, idx, ...] *= values[-1, ...]  # Ayy * Byy

                    elif not self.doublefulljones:
                        if 'pol' in self.axes_current:
                            self.amplitudes[:, idx, ...] *= values[:, ...]
                        else:
                            self.amplitudes[idx, ...] *= values[...]
                        # print(diagdiag_math)

                    elif self.doublefulljones:  # save table for doublefulljones merge
                        self.fulljones_amplitudes.update({h5_name: values})

                    else:
                        sys.exit("ERROR: Broken code. Should not be here.")
            # print(np.std(self.phases, axis=(0,2,3,4)))
            h5.close()

        return self

    def matrix_multiplication(self):
        """
        Matrix multiplication for double fulljones matrices
        """

        def reshape_dir(arr):
            """
            Reshape with dir in it
            """
            shape = list(arr.shape)
            shape.insert(1, 1)
            return arr.reshape(shape)

        # print(doublefulljones_math)

        tables = []
        for h5, t in self.fulljones_phases.items():
            if h5 in self.fulljones_amplitudes.keys():
                tables.append(self.fulljones_amplitudes[h5].astype(complex128) * exp(t * 1j))
            else:
                out = exp(t * 1j)
                out[1] = 0.
                out[2] = 0.
                tables.append(out)

        pol_ax = self.axes_final.index('pol')
        output_table = self.amplitudes.astype(complex128) * exp(self.phases * 1j)
        print("Calculating...")

        # Set template matrix
        Axx = take(output_table, indices=[0], axis=pol_ax)
        Axy = zeros(Axx.shape).astype(complex128)
        Ayx = zeros(Axx.shape).astype(complex128)
        Ayy = take(output_table, indices=[3], axis=pol_ax)

        # if sys.version_info.major > 2:
        #     print('0%', end='...')
        # else:
        #     print('Start merging ...')

        # Looping over tables to construct output
        for n, input_table in enumerate(tables):
            # if sys.version_info.major > 2:
            #     print(str(int((n+1)*100/len(tables)))+'%', end='...')

            Bxx = reshape_dir(take(input_table, indices=[0], axis=pol_ax))
            Bxy = reshape_dir(take(input_table, indices=[1], axis=pol_ax))
            Byx = reshape_dir(take(input_table, indices=[2], axis=pol_ax))
            Byy = reshape_dir(take(input_table, indices=[3], axis=pol_ax))

            # Necessary to make sure Axx, .. Ayy are not overwritten during calculations
            nAxx = Axx * Bxx + Axy * Byx
            nAxy = Axy * Byy + Axx * Bxy
            nAyx = Ayx * Bxx + Ayy * Byx
            nAyy = Ayy * Byy + Ayx * Bxy

            # Set final new values
            Axx = nAxx.copy()
            Axy = nAxy.copy()
            Ayx = nAyx.copy()
            Ayy = nAyy.copy()

        print('Done\n')

        # Construct final output matrix
        if pol_ax == 0:
            output_table[0, ...] = Axx
            output_table[1, ...] = Axy
            output_table[2, ...] = Ayx
            output_table[3, ...] = Ayy
        else:
            sys.exit('ERROR: double fulljones merge error --> polarization axis is not 0')

        # Convert output matrix to phases and amplitudes
        self.phases = angle(output_table)
        self.amplitudes = abs(output_table)

        return self

    def _DP3_order(self, soltab):
        """
        Reorder the axes in the same way as DP3 output because that is needed for other LOFAR software

        :param soltab: Solution table

        :return: New DP3 axes
        """

        if 'pol' in self.axes_final and len(self.axes_final) == 5:
            DP3_axes = ['time', 'freq', 'ant', 'dir', 'pol']
        elif 'pol' not in self.axes_final and len(self.axes_final) == 4:
            DP3_axes = ['time', 'ant', 'dir', 'freq']
            if self.phases.ndim == 5:
                self.phases = self.phases[0]
        else:
            DP3_axes = []

        if 'phase' in soltab or ('tec' in soltab and self.convert_tec) or 'rotation' in soltab:
            self.phases = reorderAxes(self.phases, self.axes_final, DP3_axes)
        elif 'amplitude' in soltab:
            self.amplitudes = reorderAxes(self.amplitudes, self.axes_final, DP3_axes)
        elif 'tec' in soltab and not self.convert_tec:
            self.tec = reorderAxes(self.tec, self.axes_final, DP3_axes)
        elif 'error' in soltab:
            self.error = reorderAxes(self.error, self.axes_final, DP3_axes)

        return DP3_axes

    def reorder_directions(self):
        """
        This method will be called when the user is using Python 2, as there was a bug in the direction that can be resolved
        with this extra step.
        """

        H = tables.open_file(self.h5name_out, 'r+')
        for solset in H.root._v_groups.keys():
            ss = H.root._f_get_child(solset)
            # Not needed if only 1 source
            if len(ss.source[:]) == 1:
                H.close()
                return self
            # Problem when source table is empty
            elif len(ss.source[:]) == 0:
                H.close()
                sys.exit('ERROR: No sources in output table ' + '/'.join([solset, 'source']))
            # No reordering needed
            elif all(ss.source[:]['name'] == ss._f_get_child(list(ss._v_groups.keys())[0]).dir[:]):
                H.close()
                return self
            # Reordering needed
            else:
                sources = array(sort(ss.source[:]), dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
                overwrite_table(H, solset, 'source', sources)
                for soltab in ss._v_groups.keys():
                    st = ss._f_get_child(soltab)
                    st.dir._f_remove()
                    H.create_array(st, 'dir', array(sources['name'], dtype='|S5'))
        H.close()
        return self

    def reduce_memory_source(self):
        """
        Store the data in 136 bytes per directions.
        Python 3 saves it automatically in more than that number (for yet some unknown reason).
        """

        T = tables.open_file(self.h5name_out, 'r+')
        for solset in T.root._v_groups.keys():
            ss = T.root._f_get_child(solset)
            if ss.source[:][0].nbytes > 140:
                overwrite_table(T, solset, 'source', ss.source[:])
        T.close()
        return self

    @staticmethod
    def keep_new_sources(current_sources, new_sources):
        """
        Remove sources from new_sources that are already in current_sources

        :param current_sources: current sources that need to be compared with new_sources
        :param new_sources: new sources to be add

        :return: New unique sources
        """

        current_sources_dir = [source[0].decode('UTF-8') for source in current_sources]
        current_sources_coor = [source[1] for source in current_sources]
        new_sources = [source for source in new_sources if source[0] not in current_sources_dir]
        del_index = []
        for coor in current_sources_coor:
            for n, source in enumerate(new_sources):
                if round(coor[0], 5) == round(source[1][0], 5) and round(coor[1], 5) == round(source[1][1], 5):
                    del_index.append(n)

        return [source for i, source in enumerate(new_sources) if i not in del_index]

    def create_new_dataset(self, solset, soltab):
        """
        Create a new dataset in the h5 table

        :param solset: solution set name
        :param soltab: solution table name
        """

        if len(self.directions.keys()) == 0:  # return if no directions
            return self

        self.h5_out = h5parm(self.h5name_out, readonly=False)
        if solset in self.h5_out.getSolsetNames():
            solsetout = self.h5_out.getSolset(solset)
        else:
            solsetout = self.h5_out.makeSolset(solset)

        new_sources = self.keep_new_sources(solsetout.obj.source[:], list(self.directions.items()))

        if len(new_sources) > 0:
            solsetout.obj.source.append(new_sources)

        axes_vals = {'dir': list(self.directions.keys()),
                     'ant': self.ant,
                     'freq': self.ax_freq,
                     'time': self.ax_time}

        DP3_axes = self._DP3_order(soltab)  # reorder the axis to DP3 style

        if 'pol' in self.axes_final:
            axes_vals.update({'pol': self.polarizations})
            self.axes_final = DP3_axes
        elif len(self.axes_final) == 4 and len(DP3_axes) > 0:
            self.axes_final = DP3_axes

        # right order vals
        axes_vals = [v[1] for v in sorted(axes_vals.items(), key=lambda pair: self.axes_final.index(pair[0]))]
        shape_axes_vals = tuple(len(s) for s in axes_vals)

        # make new solution table
        if 'phase' in soltab:
            self.phases = self.remove_invalid_values('phase', self.phases, self.axes_final)
            if shape_axes_vals != self.phases.shape:
                self.phases = reorderAxes(self.phases, self.axes_current, self.axes_final)
            weights = ones(self.phases.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('phase', axesNames=self.axes_final, axesVals=axes_vals, vals=self.phases,
                                 weights=weights)
        elif 'amplitude' in soltab:
            self.amplitudes = self.remove_invalid_values('amplitude', self.amplitudes, self.axes_final)
            if shape_axes_vals != self.amplitudes.shape:
                self.amplitudes = reorderAxes(self.amplitudes, self.axes_current, self.axes_final)
            weights = ones(self.amplitudes.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('amplitude', axesNames=self.axes_final, axesVals=axes_vals, vals=self.amplitudes,
                                 weights=weights)
        elif 'tec' in soltab and not self.convert_tec:
            self.tec = self.remove_invalid_values('tec', self.tec, self.axes_final)
            if shape_axes_vals != self.tec.shape:
                self.tec = reorderAxes(self.tec, self.axes_current, self.axes_final)
            weights = ones(self.tec.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('tec', axesNames=self.axes_final,
                                 axesVals=axes_vals,
                                 vals=self.tec, weights=weights)
        elif 'error' in soltab:
            self.error = self.remove_invalid_values('error', self.error, self.axes_final)
            if shape_axes_vals != self.error.shape:
                self.error = reorderAxes(self.error, self.axes_current, self.axes_final)
            weights = ones(self.error.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('error', axesNames=self.axes_final,
                                 axesVals=axes_vals,
                                 vals=self.error, weights=weights)

        print('DONE: {solset}/{soltab}'.format(solset=solset, soltab=soltab))

        self.h5_out.close()

        return self

    def add_empty_directions(self, add_directions=None):
        """
        Add default directions (phase all zeros, amplitude all ones)

        :param add_directions: list with directions
        """

        if not add_directions:
            return self

        h5 = h5parm(self.h5name_out, readonly=True)
        filetemp = self.h5name_out.replace('.h5', '_temph5merger.h5')
        h5_temp = h5parm(filetemp, readonly=False)
        solset = h5.getSolset(self.solset)
        solsettemp = h5_temp.makeSolset(self.solset)
        if type(add_directions[0]) == list:
            sources = list([source[1] for source in solset.obj.source[:]]) + add_directions
        else:
            sources = list([source[1] for source in solset.obj.source[:]]) + [add_directions]
        if sys.version_info.major > 2:
            sources = [(bytes('Dir' + str(n).zfill(2), 'utf-8'), list(ns)) for n, ns in enumerate(sources)]
        else:
            sources = [(bytes('Dir' + str(n).zfill(2)), list(ns)) for n, ns in enumerate(sources)]
        if len(sources) > 0:
            solsettemp.obj.source.append(sources)

        for st in h5.getSolset(self.solset).getSoltabNames():
            solutiontable = h5.getSolset(self.solset).getSoltab(st)
            axes = solutiontable.getValues()[1]
            values = solutiontable.getValues()[0]
            axes['dir'] = [ns[0] for ns in sources]
            dir_index = solutiontable.getAxesNames().index('dir')
            new_shape = list(values.shape)
            last_idx = new_shape[dir_index]
            new_idx = last_idx + len(add_directions) - 1
            new_shape[dir_index] = new_idx

            if 'phase' in st:
                values_new = zeros(tuple(new_shape))
            elif 'amplitude' in st:
                values_new = ones(tuple(new_shape))
            else:
                values_new = zeros(tuple(new_shape))

            if dir_index == 0:
                values_new[0:last_idx, ...] = values
            elif dir_index == 1:
                values_new[:, 0:last_idx, ...] = values
            elif dir_index == 2:
                values_new[:, :, 0:last_idx, ...] = values
            elif dir_index == 3:
                values_new[:, :, :, 0:last_idx, ...] = values
            elif dir_index == 4:
                values_new[:, :, :, :, 0:last_idx, ...] = values

            weights = ones(values_new.shape)
            solsettemp.makeSoltab(remove_numbers(st), axesNames=list(axes.keys()), axesVals=list(axes.values()),
                                  vals=values_new,
                                  weights=weights)

            print('Default directions added for ' + self.solset + '/' + st)
            print('Shape change: ' + str(values.shape) + ' ---> ' + str(values_new.shape))

        h5.close()
        h5_temp.close()

        os.system('rm ' + self.h5name_out + ' && mv ' + filetemp + ' ' + self.h5name_out)

        return self

    def change_pol(self, single=False, nopol=False):
        """
        Change polarization dimension
        (Standard output is poldim=1)

        :param single: if True --> leave a single pole such that values have shape=(..., 1), if False --> remove pol-axis entirely
        :param nopol: if True --> no polarization in output
        """

        T = tables.open_file(self.h5name_out, 'r+')

        # DP3 axes order
        if self.poldim == 5 or single:
            output_axes = ['time', 'freq', 'ant', 'dir', 'pol']
        else:
            output_axes = ['time', 'ant', 'dir', 'freq']

        for solset in T.root._v_groups.keys():
            ss = T.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)

                if self.poldim == 0:
                    print('No polarization in ' + solset + '/' + soltab)
                    if nopol:
                        continue
                    elif single:
                        T.create_array(st, 'pol', array([b'I'], dtype='|S2'))

                elif self.poldim == 1:
                    st.pol._f_remove()
                    T.create_array(st, 'pol', array([b'I'], dtype='|S2'))
                elif self.poldim > 1 and (single or nopol):
                    for axes in ['val', 'weight']:
                        if not all(st._f_get_child(axes)[..., 0] == \
                                   st._f_get_child(axes)[..., -1]):
                            sys.exit('WARNING: ' + '/'.join([soltab, axes]) +
                                     ' does not have the same values for XX and YY polarization.'
                                     '\nERROR: No polarization reduction will be done.'
                                     '\nERROR: Do not use --no_pol or --single_pol.')

                for axes in ['val', 'weight']:
                    if self.poldim > 0 and single:
                        newval = st._f_get_child(axes)[:, :, :, :, 0:1]
                        st.pol._f_remove()
                        T.create_array(st, 'pol', array([b'I'], dtype='|S2'))
                    elif nopol and self.poldim > 0:
                        newval = st._f_get_child(axes)[:, :, :, :, 0]
                        try:
                            st.pol._f_remove()
                        except:
                            pass
                    elif self.poldim == 0 and single:
                        newval = expand_dims(st._f_get_child(axes)[:], len(st._f_get_child(axes).shape))
                        try:
                            st.pol._f_remove()
                        except:
                            pass
                        T.create_array(st, 'pol', array([b'I'], dtype='|S2'))
                    else:
                        continue

                    valtype = str(st._f_get_child(axes).dtype)
                    if '16' in valtype:
                        atomtype = tables.Float16Atom()
                    elif '32' in valtype:
                        atomtype = tables.Float32Atom()
                    elif '64' in valtype:
                        atomtype = tables.Float64Atom()
                    else:
                        atomtype = tables.Float64Atom()

                    old_axes = [ax for ax in st._f_get_child(axes).attrs['AXES'].decode('utf8').split(',') + ['pol']
                                if ax in output_axes]
                    idx = [old_axes.index(ax) for ax in output_axes]
                    newval = transpose(newval, idx)
                    st._f_get_child(axes)._f_remove()
                    T.create_array(st, axes, newval.astype(valtype), atom=atomtype)
                    if sys.version_info.major > 2:
                        st._f_get_child(axes).attrs['AXES'] = bytes(','.join(output_axes), 'utf-8')
                    else:
                        st._f_get_child(axes).attrs['AXES'] = bytes(','.join(output_axes))

                print('Value shape after changing poldim --> ' + str(st._f_get_child('val')[:].shape))
                print('Polarization --> ' + str(st.pol[:]))

        T.close()

        return self

    def add_h5_antennas(self):
        """
        Add antennas to output table from H5 list.
        """

        print('Add antenna table from ' + self.h5_tables[0])
        T = tables.open_file(self.h5_tables[0])
        antennas = T.root.sol000.antenna[:]
        T.close()
        H = tables.open_file(self.h5name_out, 'r+')
        for solset in H.root._v_groups.keys():
            overwrite_table(H, solset, 'antenna', antennas)
        H.close()

        return self

    def add_ms_antennas(self, keep_h5_interstations=None):
        """
        Add antennas from MS.

        :param keep_h5_interstations: keep h5 international stations
        """

        print('Add antenna table from ' + self.ms[0])
        if len(self.ms) == 0:
            sys.exit("ERROR: Measurement set needed to add antennas. Use --ms.")

        t = ct.table(self.ms[0] + "::ANTENNA", ack=False)
        try:
            ms_antlist = [n.decode('utf8') for n in t.getcol('NAME')]
        except AttributeError:
            ms_antlist = t.getcol('NAME')
        ms_antpos = t.getcol('POSITION')
        ms_antennas = array([list(zip(*(ms_antlist, ms_antpos)))], dtype=[('name', 'S16'), ('position', '<f4', (3,))])
        t.close()

        H = tables.open_file(self.h5name_out, 'r+')

        for solset in H.root._v_groups.keys():
            ss = H.root._f_get_child(solset)

            F = tables.open_file(self.h5_tables[0])
            h5_antennas = F.root._f_get_child(solset).antenna[:]
            F.close()

            for soltab in ss._v_groups.keys():
                print(soltab)
                st = ss._f_get_child(soltab)
                attrsaxes = st.val.attrs['AXES']
                antenna_index = attrsaxes.decode('utf8').split(',').index('ant')
                h5_antlist = [v.decode('utf8') for v in list(st.ant[:])]

                if keep_h5_interstations:  # keep international stations if these are not in MS
                    new_antlist = [station for station in ms_antlist if 'CS' in station] + \
                                  [station for station in h5_antlist if 'ST' not in station]
                    all_antennas = [a for a in unique(append(ms_antennas, h5_antennas), axis=0) if a[0] != 'ST001']
                    antennas_new = [all_antennas[[a[0].decode('utf8') for a in all_antennas].index(na)] for na in
                                    new_antlist]  # sorting
                else:
                    new_antlist = ms_antlist
                    antennas_new = ms_antennas

                st.ant._f_remove()
                overwrite_table(H, solset, 'antenna', antennas_new)
                H.create_array(st, 'ant', array(list(new_antlist), dtype='|S16'))

                try:
                    superstation_index = h5_antlist.index('ST001')
                    superstation = True
                except ValueError:
                    superstation = False
                    print('No super station (ST001) in antennas.')

                for axes in ['val', 'weight']:
                    assert axes in list(
                        st._v_children.keys()), axes + ' not in .root.' + solset + '.' + soltab + ' (not in axes)'
                    h5_values = st._f_get_child(axes)[:]
                    shape = list(h5_values.shape)
                    shape[antenna_index] = len(new_antlist)
                    ms_values = zeros(shape)

                    for idx, antenna in enumerate(new_antlist):
                        if antenna in h5_antlist:
                            idx_h5 = h5_antlist.index(antenna)
                            if antenna_index == 0:
                                ms_values[idx, ...] += h5_values[idx_h5, ...]
                            elif antenna_index == 1:
                                ms_values[:, idx, ...] += h5_values[:, idx_h5, ...]
                            elif antenna_index == 2:
                                ms_values[:, :, idx, ...] += h5_values[:, :, idx_h5, ...]
                            elif antenna_index == 3:
                                ms_values[:, :, :, idx, ...] += h5_values[:, :, :, idx_h5, ...]
                            elif antenna_index == 4:
                                ms_values[:, :, :, :, idx, ...] += h5_values[:, :, :, :, idx_h5, ...]
                        elif 'CS' in antenna and superstation:  # core stations
                            if antenna_index == 0:
                                ms_values[idx, ...] += h5_values[superstation_index, ...]
                            elif antenna_index == 1:
                                ms_values[:, idx, ...] += h5_values[:, superstation_index, ...]
                            elif antenna_index == 2:
                                ms_values[:, :, idx, ...] += h5_values[:, :, superstation_index, ...]
                            elif antenna_index == 3:
                                ms_values[:, :, :, idx, ...] += h5_values[:, :, :, superstation_index, ...]
                            elif antenna_index == 4:
                                ms_values[:, :, :, :, idx, ...] += h5_values[:, :, :, :, superstation_index, ...]
                        elif antenna not in h5_antlist and ('amplitude' in soltab or axes == 'weight') \
                                and 'RS' not in antenna and 'CS' not in antenna:
                            if axes == 'val':
                                print('Add ' + antenna + ' to output H5 from MS')
                            if antenna_index == 0:
                                ms_values[idx, ...] = 1
                            elif antenna_index == 1:
                                ms_values[:, idx, ...] = 1
                            elif antenna_index == 2:
                                ms_values[:, :, idx, ...] = 1
                            elif antenna_index == 3:
                                ms_values[:, :, :, idx, ...] = 1
                            elif antenna_index == 4:
                                ms_values[:, :, :, :, idx, ...] = 1

                    valtype = str(st._f_get_child(axes).dtype)
                    if '16' in valtype:
                        atomtype = tables.Float16Atom()
                    elif '32' in valtype:
                        atomtype = tables.Float32Atom()
                    elif '64' in valtype:
                        atomtype = tables.Float64Atom()
                    else:
                        atomtype = tables.Float64Atom()

                    st._f_get_child(axes)._f_remove()
                    H.create_array(st, axes, ms_values.astype(valtype), atom=atomtype)
                    st._f_get_child(axes).attrs['AXES'] = attrsaxes
                print('Value shape after --> ' + str(st.val.shape))

        H.close()

        return self

    def add_template(self):
        """
        Make template for phase000 and/or amplitude000 if missing.
        """

        H = tables.open_file(self.h5name_out, 'r+')
        soltabs = list(H.root.sol000._v_groups.keys())
        if 'amplitude000' not in soltabs:
            if 'phase000' in soltabs:
                self.amplitudes = ones(H.root.sol000.phase000.val.shape)
                if len(self.polarizations) == 4:
                    self.amplitudes[..., 1] = 0.
                    self.amplitudes[..., 2] = 0.
                # print("HE: Skip amplitude000 template!")
                self.create_new_dataset('sol000', 'amplitude')
        if 'phase000' not in soltabs:
            if 'amplitude000' in soltabs:
                self.phases = zeros(H.root.sol000.amplitude000.val.shape)
                self.create_new_dataset('sol000', 'phase')
        H.close()

        return self

    def flag_stations(self):
        """
        Propagate flagged station input to output.
        """

        H = tables.open_file(self.h5name_out, 'r+')
        for input_h5 in self.h5_tables:
            T = tables.open_file(input_h5)
            for solset in T.root._v_groups.keys():
                ss = T.root._f_get_child(solset)
                for soltab in ss._v_groups.keys():
                    st = ss._f_get_child(soltab)
                    weight = st.weight
                    antennas_in = list(st.ant[:])
                    axes = make_utf8(weight.attrs['AXES']).split(',')
                    if 'tec' in soltab:
                        soltab = 'phase' + soltab[-3:]
                    axes_output = make_utf8(
                        H.root._f_get_child(solset)._f_get_child(soltab).weight.attrs['AXES']).split(',')
                    ant_index = axes.index('ant')
                    weight = reorderAxes(weight, axes, [i for i in axes_output if i in axes])
                    for a in range(weight.shape[ant_index]):
                        antenna = antennas_in[a]
                        if sum(take(weight, [a], ant_index)) == 0.:
                            if soltab in list(H.root._f_get_child(solset)._v_groups.keys()):
                                st_out = H.root._f_get_child(solset)._f_get_child(soltab)
                                if antenna in list(st_out.ant[:]):
                                    if ant_index == 0:
                                        st_out.weight[a, ...] = 0.
                                    elif ant_index == 1:
                                        st_out.weight[:, a, ...] = 0.
                                    elif ant_index == 2:
                                        st_out.weight[:, :, a, ...] = 0.
                                    elif ant_index == 3:
                                        st_out.weight[:, :, :, a, ...] = 0.
                                    elif ant_index == 4:
                                        st_out.weight[:, :, :, :, a, ...] = 0.

            T.close()
        H.close()
        return self

    def add_weights(self):
        """
        Upsample weights (propagate flags to weights)

        """

        print("\nPropagating weights in:")

        H = tables.open_file(self.h5name_out, 'r+')
        for solset in H.root._v_groups.keys():
            ss = H.root._f_get_child(solset)
            for n, soltab in enumerate(ss._v_groups.keys()):
                print(soltab + ', from:')
                st = ss._f_get_child(soltab)
                shape = st.val.shape
                weight_out = ones(shape)
                axes_new = make_utf8(st.val.attrs["AXES"]).split(',')
                for m, input_h5 in enumerate(self.h5_tables):

                    print(input_h5)
                    T = tables.open_file(input_h5)
                    if soltab not in list(T.root._f_get_child(solset)._v_groups.keys()):
                        T.close()
                        continue
                    st2 = T.root._f_get_child(solset)._f_get_child(soltab)
                    axes = make_utf8(st2.val.attrs["AXES"]).split(',')
                    weight = st2.weight[:]
                    weight = reorderAxes(weight, axes, [a for a in axes_new if a in axes])

                    # important to do the following in case the input tables are not all different directions
                    m = min(weight_out.shape[axes.index('dir')] - 1, m)
                    if self.merge_all_in_one:
                        m = 0

                    newvals = self._interp_along_axis(weight, st2.time[:], st.time[:], axes_new.index('time'),
                                                      fill_value=1.).astype(int)
                    newvals = self._interp_along_axis(newvals, st2.freq[:], st.freq[:], axes_new.index('freq'),
                                                      fill_value=1.).astype(int)

                    if weight.shape[-2] != 1 and len(weight.shape) == 5:
                        print("Merge multi-dir weights")
                        if weight.shape[-2] != weight_out.shape[-2]:
                            print(weight.shape, weight_out.shape)
                            sys.exit("ERROR: multi-dirs do not have equal shape.")
                        if 'pol' not in axes and 'pol' in axes_new:
                            newvals = expand_dims(newvals, axis=axes_new.index('pol'))
                        for n in range(weight_out.shape[-1]):
                            print(weight_out.shape, newvals.shape)
                            weight_out[..., n] *= newvals[..., -1]

                    elif weight.ndim != weight_out.ndim:  # not the same value shape
                        if 'pol' not in axes and 'pol' in axes_new:
                            newvals = expand_dims(newvals, axis=axes_new.index('pol'))
                            if newvals.shape[-1] != weight_out.shape[-1]:
                                temp = ones(shape)
                                for n in range(temp.shape[-1]):
                                    temp[..., m, n] *= newvals[..., 0, 0]
                                newvals = temp
                            weight_out *= newvals
                        else:
                            sys.exit('ERROR: Upsampling of weights bug due to unexpected missing axes.\n axes from '
                                     + input_h5 + ': ' + str(axes) + '\n axes from '
                                     + self.h5name_out + ': ' + str(axes_new) + '.\n'
                                     + debug_message)
                    elif set(axes) == set(axes_new) and 'pol' in axes:  # same axes
                        pol_index = axes_new.index('pol')
                        if weight_out.shape[pol_index] == newvals.shape[pol_index]:  # same pol numbers
                            weight_out[:, :, :, m, ...] *= newvals[:, :, :, 0, ...]
                        else:  # not the same polarization axis
                            if newvals.shape[pol_index] != newvals.shape[-1]:
                                sys.exit('ERROR: Upsampling of weights bug due to polarization axis mismatch.\n'
                                         + debug_message)
                            if newvals.shape[pol_index] == 1:  # new values have only 1 pol axis
                                for i in range(weight_out.shape[pol_index]):
                                    weight_out[:, :, :, m, i] *= newvals[:, :, :, 0, 0]
                            elif newvals.shape[pol_index] == 2 and weight_out.shape[pol_index] == 1:
                                for i in range(newvals.shape[pol_index]):
                                    weight_out[:, :, :, m, 0] *= newvals[:, :, :, 0, i]
                            elif newvals.shape[pol_index] == 2 and weight_out.shape[pol_index] == 4:
                                weight_out[:, :, :, m, 0] *= newvals[:, :, :, 0, 0]
                                weight_out[:, :, :, m, 1] *= newvals[:, :, :, 0, 0] * newvals[:, :, :, 0, -1]
                                weight_out[:, :, :, m, 2] *= newvals[:, :, :, 0, 0] * newvals[:, :, :, 0, -1]
                                weight_out[:, :, :, m, -1] *= newvals[:, :, :, 0, -1]
                            else:
                                sys.exit('ERROR: Upsampling of weights bug due to unexpected polarization mismatch.\n'
                                         + debug_message)
                    elif set(axes) == set(axes_new) and 'pol' not in axes:  # same axes but no pol
                        dirind = axes_new.index('dir')
                        if weight_out.shape[dirind] != newvals.shape[dirind] and newvals.shape[dirind] == 1:
                            if len(weight_out.shape) == 4:
                                weight_out[:, :, m, :] *= newvals[:, :, 0, :]
                            elif len(weight_out.shape) == 5:
                                weight_out[:, :, :, m, :] *= newvals[:, :, :, 0, :]
                        elif weight_out.shape[dirind] != newvals.shape[dirind]:
                            sys.exit(
                                'ERROR: Upsampling of weights because same direction exists multiple times in input h5 (verify and/or remove --propagate_flags)')
                        else:
                            weight_out *= newvals

                    else:
                        sys.exit('ERROR: Upsampling of weights bug due to unexpected missing axes.\n axes from '
                                 + input_h5 + ': ' + str(axes) + '\n axes from '
                                 + self.h5name_out + ': ' + str(axes_new) + '.\n'
                                 + debug_message)

                    T.close()
                st.weight[:] = weight_out

        H.close()

        print('\n')
        return self

    def format_tables(self):
        """
        Format direction tables (making sure dir and sources are the same).
        """

        H = tables.open_file(self.h5name_out, 'r+')
        for solset in H.root._v_groups.keys():
            ss = H.root._f_get_child(solset)
            sources = ss.source[:]['name']
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)
                dirs = st._f_get_child('dir')[:]
                if len(sources[:]) > len(dirs):
                    difference = list(set(sources) - set(dirs))

                    newdir = list(st._f_get_child('dir')[:]) + difference
                    st._f_get_child('dir')._f_remove()
                    H.create_array(st, 'dir', array(newdir).astype('|S5'))

                    dir_ind = st.val.attrs['AXES'].decode('utf8').split(',').index('dir')

                    for axes in ['val', 'weight']:
                        axs = st._f_get_child(axes).attrs['AXES']
                        newval = st._f_get_child(axes)[:]

                        shape = list(newval.shape)
                        for _ in difference:
                            shape[dir_ind] = 1
                            if 'amplitude' in soltab or axes == 'weight':
                                newval = append(newval, ones(shape),
                                                axis=dir_ind)
                            else:
                                newval = append(newval, zeros(shape),
                                                axis=dir_ind)

                        valtype = str(st._f_get_child(axes).dtype)
                        if '16' in valtype:
                            atomtype = tables.Float16Atom()
                        elif '32' in valtype:
                            atomtype = tables.Float32Atom()
                        elif '64' in valtype:
                            atomtype = tables.Float64Atom()
                        else:
                            atomtype = tables.Float64Atom()

                        st._f_get_child(axes)._f_remove()
                        H.create_array(st, axes, newval.astype(valtype), atom=atomtype)
                        st._f_get_child(axes).attrs['AXES'] = axs
                elif len(sources[:]) < len(dirs):
                    output_check(self.h5name_out)
        H.close()

        return self


def _create_h5_name(h5_name):
    """
    Correct such that output has always .h5 as extension

    :param h5_name: h5 input file name

    :return: h5 corrected name with .h5 extension
    """

    if '.h5' != h5_name[-3:]:
        h5_name += '.h5'
    return h5_name


def _change_solset(h5, solset_in, solset_out, delete=True, overwrite=True):
    """
    This function is to change the solset numbers.

    :param h5: h5 input file
    :param solset_in: solution set name input
    :param solset_out: solution set name output
    :param delete: delete solution set input
    :param overwrite: overwrite if solution set out already exists

    1) Copy solset_in to solset_out (overwriting if overwrite==True)
    2) Delete solset_in if delete==True
    """

    H = tables.open_file(h5, 'r+')
    H.root._f_get_child(solset_in)._f_copy(H.root, newname=solset_out, overwrite=overwrite, recursive=True)
    print('Succesfully copied ' + solset_in + ' to ' + solset_out)
    if delete:
        H.root._f_get_child(solset_in)._f_remove(recursive=True)
        print('Removed ' + solset_in + ' in output')
    H.close()

    return


def output_check(h5):
    """
    Validate output by checking and comparing tables.

    :param h5: h5 file

    :return: True if output passed
    """

    print('\nChecking output...')

    H = tables.open_file(h5)

    # check number of solset
    assert len(list(H.root._v_groups.keys())) == 1, \
        'More than 1 solset in ' + str(list(H.root._v_groups.keys())) + '. Only 1 is allowed for h5_merger.py.'

    for solset in H.root._v_groups.keys():

        # check sol00.. name
        assert 'sol' in solset, solset + ' is a wrong solset name, should be sol***'
        ss = H.root._f_get_child(solset)

        # check antennas
        antennas = ss.antenna
        assert antennas.attrs.FIELD_0_NAME == 'name', 'No name in ' + '/'.join([solset, 'antenna'])
        assert antennas.attrs.FIELD_1_NAME == 'position', 'No coordinate in ' + '/'.join([solset, 'antenna'])

        # check sources
        sources = ss.source
        assert sources.attrs.FIELD_0_NAME == 'name', 'No name in ' + '/'.join([solset, 'source'])
        assert sources.attrs.FIELD_1_NAME == 'dir', 'No coordinate in ' + '/'.join([solset, 'source'])

        for soltab in ss._v_groups.keys():
            st = ss._f_get_child(soltab)
            assert st.val.shape == st.weight.shape, \
                'weight ' + str(st.weight.shape) + ' and values ' + str(st.val.shape) + ' do not have same shape'

            # check if pol and/or dir are missing
            for pd in ['pol', 'dir']:
                assert not (st.val.ndim == 5 and pd not in list(st._v_children.keys())), \
                    '/'.join([solset, soltab, pd]) + ' is missing'

            # check if freq, time, and ant arrays are missing
            for fta in ['freq', 'time', 'ant']:
                assert fta in list(st._v_children.keys()), \
                    '/'.join([solset, soltab, fta]) + ' is missing'

            # check if val and weight have AXES
            for vw in ['val', 'weight']:
                assert 'AXES' in st._f_get_child(vw).attrs._f_list("user"), \
                    'AXES missing in ' + '/'.join([solset, soltab, vw])

            # check if dimensions of values match with length of arrays
            for ax_index, ax in enumerate(st.val.attrs['AXES'].decode('utf8').split(',')):
                assert st.val.shape[ax_index] == len(st._f_get_child(ax)[:]), \
                    ax + ' length is not matching with dimension from val in ' + '/'.join([solset, soltab, ax])

                # check if ant and antennas have equal sizes
                if ax == 'ant':
                    assert len(antennas[:]) == len(st._f_get_child(ax)[:]), \
                        '/'.join([solset, 'antenna']) + ' and ' + '/'.join(
                            [solset, soltab, ax]) + ' do not have same length'

                # check if dir and sources have equal sizes
                if ax == 'dir':
                    assert len(sources[:]) == len(st._f_get_child(ax)[:]), \
                        '/'.join([solset, 'source']) + ' and ' + '/'.join(
                            [solset, soltab, ax]) + ' do not have same length'

            # check if phase and amplitude have same shapes
            for soltab1 in ss._v_groups.keys():
                if ('phase' in soltab or 'amplitude' in soltab) and ('phase' in soltab1 or 'amplitude' in soltab1):
                    st1 = ss._f_get_child(soltab1)
                    assert st.val.shape == st1.val.shape, \
                        '/'.join([solset, soltab, 'val']) + ' shape: ' + str(st.weight.shape) + \
                        '/'.join([solset, soltab1, 'val']) + ' shape: ' + str(st1.weight.shape)

    H.close()

    print('Awesome! Output has all necessary information and correct dimensions.')

    return True


class PolChange:
    """
    This Python class helps to convert polarization from linear to circular or vice versa.
    """

    def __init__(self, h5_in, h5_out):
        """
        :param h5_in: h5 input name
        :param h5_out: h5 output name
        """

        self.h5in_name = h5_in
        self.h5out_name = h5_out
        self.h5_in = h5parm(h5_in, readonly=True)
        self.h5_out = h5parm(h5_out, readonly=False)
        self.axes_names = ['time', 'freq', 'ant', 'dir', 'pol']

    @staticmethod
    def lin2circ(G):
        """
        Convert linear polarization to circular polarization

        RR = XX - iXY + iYX + YY
        RL = XX + iXY + iYX - YY
        LR = XX - iXY - iYX - YY
        LL = XX + iXY - iYX + YY

        :param G: Linear polarized Gain

        :return: Circular polarized Gain
        """

        RR = (G[..., 0] + G[..., -1])
        RL = (G[..., 0] - G[..., -1])
        LR = (G[..., 0] - G[..., -1])
        LL = (G[..., 0] + G[..., -1])

        if G.shape[-1] == 4:
            RR += 1j * (G[..., 2] - G[..., 1])
            RL += 1j * (G[..., 2] + G[..., 1])
            LR -= 1j * (G[..., 2] + G[..., 1])
            LL += 1j * (G[..., 1] - G[..., 2])

        RR /= 2
        RL /= 2
        LR /= 2
        LL /= 2

        G_new = zeros(G.shape[0:-1] + (4,)).astype(complex128)

        G_new[..., 0] += RR
        G_new[..., 1] += RL
        G_new[..., 2] += LR
        G_new[..., 3] += LL

        G_new = where(abs(G_new) < 10 * finfo(float).eps, 0, G_new)

        return G_new

    @staticmethod
    def circ2lin(G):
        """
        Convert circular polarization to linear polarization

        XX = RR + RL + LR + LL
        XY = iRR - iRL + iLR - iLL
        YX = -iRR - iRL + iLR + iLL
        YY = RR - RL - LR + LL

        :param G: Circular polarized Gain

        :return: linear polarized Gain
        """

        XX = (G[..., 0] + G[..., -1])
        XY = 1j * (G[..., 0] - G[..., -1])
        YX = 1j * (G[..., -1] - G[..., 0])
        YY = (G[..., 0] + G[..., -1])

        if G.shape[-1] == 4:
            XX += (G[..., 2] + G[..., 1])
            XY += 1j * (G[..., 2] - G[..., 1])
            YX += 1j * (G[..., 2] - G[..., 1])
            YY -= (G[..., 1] + G[..., 2])

        XX /= 2
        XY /= 2
        YX /= 2
        YY /= 2

        G_new = zeros(G.shape[0:-1] + (4,)).astype(complex128)

        G_new[..., 0] += XX
        G_new[..., 1] += XY
        G_new[..., 2] += YX
        G_new[..., 3] += YY

        G_new = where(abs(G_new) < 10 * finfo(float).eps, 0, G_new)

        return G_new

    @staticmethod
    def add_polarization(values, dim_pol):
        """
        Add extra polarization if there is no polarization

        :param values: values which need to get a polarization
        :param dim_pol: number of dimensions

        :return: input values with extra polarization axis
        """

        values_new = ones(values.shape + (dim_pol,))
        for i in range(dim_pol):
            values_new[..., i] = values

        return values_new

    def create_template(self, soltab):
        """
        Make template of the gains with only ones

        :param soltab: solution table (phase, amplitude)
        """

        self.G, self.axes_vals = array([]), OrderedDict()
        for ss in self.h5_in.getSolsetNames():
            for st in self.h5_in.getSolset(ss).getSoltabNames():
                solutiontable = self.h5_in.getSolset(ss).getSoltab(st)
                if soltab in st:
                    try:
                        if 'pol' in solutiontable.getAxesNames():
                            values = reorderAxes(solutiontable.getValues()[0],
                                                 solutiontable.getAxesNames(),
                                                 self.axes_names)
                            self.G = ones(values.shape).astype(complex128)
                        else:
                            values = reorderAxes(solutiontable.getValues()[0],
                                                 solutiontable.getAxesNames(),
                                                 self.axes_names[0:-1])
                            self.G = ones(values.shape + (2,)).astype(complex128)
                    except:
                        sys.exit('ERROR: Received ' + str(solutiontable.getAxesNames()) +
                                 ', but expect at least [time, freq, ant, dir] or [time, freq, ant, dir, pol]')

                    self.axes_vals = {'time': solutiontable.getAxisValues('time'),
                                      'freq': solutiontable.getAxisValues('freq'),
                                      'ant': solutiontable.getAxisValues('ant'),
                                      'dir': solutiontable.getAxisValues('dir'),
                                      'pol': ['XX', 'XY', 'YX', 'YY']}
                    break

        print('Value shape {soltab} before --> {shape}'.format(soltab=soltab, shape=self.G.shape))

        return self

    def add_tec(self, solutiontable):
        """
        Add TEC

        :param solutiontable: the solution table for the TEC
        """

        tec_axes_names = [ax for ax in self.axes_names if solutiontable.getAxesNames()]
        tec = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), tec_axes_names)
        if 'freq' in solutiontable.getAxesNames():
            axes_vals_tec = {'time': solutiontable.getAxisValues('time'),
                             'freq': solutiontable.getAxisValues('freq'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'dir': solutiontable.getAxisValues('dir')}
        else:
            axes_vals_tec = {'dir': solutiontable.getAxisValues('dir'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'time': solutiontable.getAxisValues('time')}
        if 'pol' in solutiontable.getAxesNames():
            if tec.shape[-1] == 2:
                axes_vals_tec.update({'pol': ['XX', 'YY']})
            elif tec.shape[-1] == 4:
                axes_vals_tec.update({'pol': ['XX', 'XY', 'YX', 'YY']})
        axes_vals_tec = [v[1] for v in
                         sorted(axes_vals_tec.items(), key=lambda pair: self.axes_names.index(pair[0]))]
        self.solsetout.makeSoltab('tec', axesNames=tec_axes_names, axesVals=axes_vals_tec, vals=tec,
                                  weights=ones(tec.shape))

        return self

    def create_new_gain_table(self, lin2circ, circ2lin):
        """
        Create new gain tables with polarization conversion.

        :param lin2circ: boolean for linear to circular conversion
        :param circ2lin: boolean for circular to linear conversion
        """

        for ss in self.h5_in.getSolsetNames():

            self.solsetout = self.h5_out.makeSolset(ss)

            for n, st in enumerate(self.h5_in.getSolset(ss).getSoltabNames()):
                solutiontable = self.h5_in.getSolset(ss).getSoltab(st)

                print('{ss}/{st} from {h5}'.format(ss=ss, st=st, h5=self.h5in_name))
                if 'phase' in st:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(),
                                             self.axes_names)
                        self.G *= exp(values * 1j)
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(),
                                             self.axes_names[0:-1])
                        self.G *= exp(self.add_polarization(values, 2) * 1j)

                elif 'amplitude' in st:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(),
                                             self.axes_names)
                        self.G *= values
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(),
                                             self.axes_names[0:-1])
                        self.G *= self.add_polarization(values, 2)

                elif 'tec' in st:
                    self.add_tec(solutiontable)
                else:
                    print("WARNING: didn't include {st} in this h5_merger.py version yet".format(st=st) +
                          "\nPlease create a ticket on github if this needs to be changed")

                if n == 0:
                    weights = ones(values.shape)
                weight = solutiontable.getValues(weight=True, retAxesVals=False)
                if 'pol' in solutiontable.getAxesNames():
                    weight = reorderAxes(weight, solutiontable.getAxesNames(), self.axes_names)
                else:
                    weight = reorderAxes(weight, solutiontable.getAxesNames(), self.axes_names[0:-1])
                weights *= weight

            if lin2circ:
                print(lin2circ_math)
                G_new = self.lin2circ(self.G)
            elif circ2lin:
                print(circ2lin_math)
                G_new = self.circ2lin(self.G)
            else:
                sys.exit('ERROR: No conversion given.')
            print('Value shape after --> {shape}'.format(shape=G_new.shape))

            phase = angle(G_new)
            amplitude = abs(G_new)

            # upsample weights
            if phase.shape != weights.shape:
                new_weights = ones(phase.shape)
                for s in range(new_weights.shape[-1]):
                    if len(weights.shape) == 4:
                        new_weights[..., s] *= weights[:]
                    elif len(weights.shape) == 5:
                        new_weights[..., s] *= weights[..., 0]
                weights = new_weights

            self.axes_vals = [v[1] for v in
                              sorted(self.axes_vals.items(), key=lambda pair: self.axes_names.index(pair[0]))]

            self.solsetout.makeSoltab('phase', axesNames=self.axes_names, axesVals=self.axes_vals, vals=phase,
                                      weights=weights)
            print('Created new phase solutions')

            self.solsetout.makeSoltab('amplitude', axesNames=self.axes_names, axesVals=self.axes_vals, vals=amplitude,
                                      weights=weights)
            print('Created new amplitude solutions')

            # convert the polarization names, such that it is clear if the h5 is in circular or linear polarization
            for solset in self.h5_out.getSolsetNames():
                ss = self.h5_out.getSolset(solset)
                for soltab in ss.getSoltabNames():
                    st = ss.getSoltab(soltab)
                    if 'pol' in st.getAxesNames():
                        pols = st.getAxisValues('pol')
                        if len(pols) == 2:
                            if lin2circ:
                                st.setAxisValues('pol', ['RR', 'LL'])
                            elif circ2lin:
                                st.setAxisValues('pol', ['XX', 'YY'])
                        if len(pols) == 4:
                            if lin2circ:
                                st.setAxisValues('pol', ['RR', 'RL', 'LR', 'LL'])
                            elif circ2lin:
                                st.setAxisValues('pol', ['XX', 'XY', 'YX', 'YY'])

        self.h5_in.close()
        self.h5_out.close()

        return self

    def add_antenna_source_tables(self):
        """
        Add antenna and source table to output file
        """

        T = tables.open_file(self.h5in_name)
        H = tables.open_file(self.h5out_name, 'r+')

        for solset in T.root._v_groups.keys():
            ss = T.root._f_get_child(solset)
            overwrite_table(H, solset, 'antenna', ss.antenna[:])
            overwrite_table(H, solset, 'source', ss.source[:])

        T.close()
        H.close()

        return


def h5_check(h5):
    """
    With this function you can print a summary from the h5 solution file

    :param h5: h5 solution file
    """

    print(
        '%%%%%%%%%%%%%%%%%%%%%%%%%\nSOLUTION FILE CHECK START\n%%%%%%%%%%%%%%%%%%%%%%%%%\n\nSolution file name:\n' + h5)
    H = tables.open_file(h5)
    solsets = list(H.root._v_groups.keys())
    print('\nFollowing solution sets in ' + h5 + ':\n' + '\n'.join(solsets))
    for solset in solsets:
        ss = H.root._f_get_child(solset)
        soltabs = list(ss._v_groups.keys())
        print('\nFollowing solution tables in ' + solset + ':\n' + '\n'.join(soltabs))
        print('\nFollowing stations in ' + solset + ':\n' + ', '.join(
            [make_utf8(a) for a in list(ss.antenna[:]['name'])]))
        print('\nFollowing sources in ' + solset + ':\n' + '\n'.join(
            [make_utf8(a['name']) + '-->' + str([a['dir'][0], a['dir'][1]]) for a in list(ss.source[:])]))
        for soltab in soltabs:
            st = ss._f_get_child(soltab)
            for table in ['val', 'weight']:
                axes = make_utf8(st._f_get_child(table).attrs["AXES"])
                print('\n' + '/'.join([solset, soltab, table]) + ' axes:\n' +
                      axes)
                print('/'.join([solset, soltab, table]) + ' shape:\n' +
                      str(st._f_get_child(table).shape))
                if table == 'weight':
                    weights = st._f_get_child(table)[:]
                    element_sum = 1
                    for w in weights.shape: element_sum *= w
                    flagged = round(sum(weights == 0.) / element_sum * 100, 2)
                    print('/'.join([solset, soltab, table]) + ' flagged:\n' + str(flagged) + '%')
                if 'pol' in axes:
                    print('/'.join([solset, soltab, 'pol']) + ':\n' + ','.join(
                        [make_utf8(p) for p in list(st._f_get_child('pol')[:])]))
                if 'time' in axes:
                    time = st._f_get_child('time')[:]
                    # print('/'.join([solset, soltab, 'time']) + ' start:\n' + str(time[0]))
                    # print('/'.join([solset, soltab, 'time']) + ' end:\n' + str(time[-1]))
                    if len(st._f_get_child('time')[:]) > 1:
                        print('/'.join([solset, soltab, 'time']) + ' time resolution:\n' + str(
                            diff(st._f_get_child('time')[:])[0]))
                if 'freq' in axes:
                    freq = st._f_get_child('freq')[:]
                    # print('/'.join([solset, soltab, 'freq']) + ' start:\n' + str(freq[0]))
                    # print('/'.join([solset, soltab, 'freq']) + ' end:\n' + str(freq[-1]))
                    if len(st._f_get_child('freq')[:]) > 1:
                        print('/'.join([solset, soltab, 'freq']) + ' resolution:\n' + str(
                            diff(st._f_get_child('freq')[:])[0]))

    H.close()
    print('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nSOLUTION FILE CHECK FINISHED\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

    return


def _checknan_input(h5):
    """
    Check h5 on nan or 0 values

    :param h5: h5parm file
    """

    H = tables.open_file(h5)

    try:
        amp = H.root.sol000.amplitude000.val[:]
        print('Amplitude:')
        print(amp[(~isfinite(amp)) | (amp == 0.)])
        del amp
    except:
        print('H.root.amplitude000 does not exist')

    try:
        phase = H.root.sol000.phase000.val[:]
        print('Phase:')
        print(phase[(~isfinite(phase)) | (phase == 0.)])
        del phase
    except:
        print('H.root.sol000.phase000 does not exist')

    H.close()

    return


def _degree_to_radian(d):
    """
    Convert degree to radio

    :param d: value in degrees

    :return: return value in radian
    """

    return pi * d / 180


def move_source_in_sourcetable(h5, overwrite=False, dir_idx=None, dra_degrees=0, ddec_degrees=0):
    """
    Change source table for specific direction

    :param overwrite: overwrite input file. If False -> add replace .h5 with _update.h5
    :param dir_idx: directions index
    :param ra_degrees: change right ascension degrees
    :param dec_degrees: change declination degrees
    """

    if not overwrite:
        os.system('cp ' + h5 + ' ' + h5.replace('.h5', '_upd.h5'))
        h5 = h5.replace('.h5', '_upd.h5')
    H = tables.open_file(h5, 'r+')
    sources = H.root.sol000.source[:]
    sources[dir_idx][1][0] += _degree_to_radian(dra_degrees)
    sources[dir_idx][1][1] += _degree_to_radian(ddec_degrees)
    overwrite_table(H, 'sol000', 'source', sources)
    H.close()

    return


def check_freq_overlap(h5_tables):
    """
    Verify if the frequency bands between the h5 tables overlap
    :param h5_tables: h5parm input tables
    """
    for h51 in h5_tables:
        H = tables.open_file(h51)
        for h52 in h5_tables:
            F = tables.open_file(h52)
            try:
                if h51 != h52 and F.root.sol000.phase000.freq[:].max() < H.root.sol000.phase000.freq[:].min():
                    print(
                        "WARNING: frequency bands between " + h51 + " and " + h52 + " do not overlap, you might want to use "
                                                                                    "--freq_concat to merge over different frequency bands")
            except:
                try:
                    if h51 != h52 and F.root.sol000.amplitude000.freq[:].max() < H.root.sol000.amplitude000.freq[
                                                                                 :].min():
                        print(
                            "WARNING: frequency bands between " + h51 + " and " + h52 + " do not overlap, you might want to use "
                                                                                        "--freq_concat to merge over different frequency bands")
                except:
                    pass
            F.close()
        H.close()
    return


def check_time_overlap(h5_tables):
    """
    Verify if the time slots between the h5 tables overlap

    :param h5_tables: h5parm input tables
    """
    for h51 in h5_tables:
        H = tables.open_file(h51)
        for h52 in h5_tables:
            F = tables.open_file(h52)
            try:
                if h51 != h52 and F.root.sol000.phase000.time[:].max() < H.root.sol000.phase000.time[:].min():
                    print(
                        "WARNING: time slots between " + h51 + " and " + h52 + " do not overlap, might result in interpolation issues")
            except:
                try:
                    if h51 != h52 and F.root.sol000.amplitude000.time[:].max() < H.root.sol000.amplitude000.time[
                                                                                 :].min():
                        print(
                            "WARNING: time slots between " + h51 + " and " + h52 + " do not overlap, might result in interpolation issues")
                except:
                    pass
            F.close()
        H.close()
    return


def running_mean(nparray, avgfactor):
    """
    Running mean over numpy axis

    :param nparray: numpy array
    :param avgfactor: averaging factor

    :return: running mean
    """
    cs = cumsum(insert(nparray, 0, 0))

    return (cs[avgfactor:] - cs[:-avgfactor]) / float(avgfactor)


def merge_h5(h5_out=None, h5_tables=None, ms_files=None, h5_time_freq=None, convert_tec=True, merge_all_in_one=False,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None, use_solset='sol000',
             filtered_dir=None, add_cs=None, add_ms_stations=None, check_output=None, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, freq_concat=None, time_concat=None,
             no_antenna_crash=None,
             output_summary=None, min_distance=0.):
    """
    Main function that uses the class MergeH5 to merge h5 tables.

    :param h5_out (string): h5 table name out
    :param h5_tables (string or list): h5 tables to merge
    :param ms_files (string or list): ms files
    :param h5_time_freq (str or boolean): h5 file to take freq and time axis from
    :param freq_av (int): averaging factor frequency axis
    :param time_av (int): averaging factor time axis
    :param convert_tec (boolean): convert TEC to phase or not
    :param merge_all_in_one: merge all in one direction
    :param lin2circ: boolean for linear to circular conversion
    :param circ2lin: boolean for circular to linear conversion
    :param add_directions: add default directions by giving a list of directions (coordinates)
    :param single_pol: only one polarization
    :param no_pol: no polarization output
    :param use_solset: use specific solset number
    :param filtered_dir: filter a specific list of directions from h5 file. Only lists allowed.
    :param add_cs: use MS to replace super station with core station
    :param add_ms_stations: return only stations from Measurement set
    :param check_output: check if output has all correct output information
    :param check_flagged_station: check if complete input stations are flagged, if so flag same stations in output
    :param propagate_flags: interpolate weights and return in output file
    :param freq_concat: concat freq blocks
    :param time_concat: concat time blocks
    :param no_antenna_crash: do not crash if antenna tables are not the same between h5s
    :param output_summary: print solution file output
    :param min_distance: Minimal coordinate distance between sources for merging directions (in degrees). If smaller, directions will be merged together.
    """

    tables.file._open_files.close_all()
    print(no_antenna_crash)

    if type(h5_tables) == str:
        h5_tables = glob(h5_tables)

    if h5_out is None:
        for h5check in h5_tables:
            h5_check(h5check)
        return

    print('\n##################################\nSTART MERGE HDF5 TABLES FOR LOFAR\n##################################'
          '\n\nMerging the following tables:\n' + '\n'.join(h5_tables) + '\n')

    # Make sure that the h5 output file is not in the input list (as that will cause conflicts)
    if h5_out in h5_tables:
        sys.exit('ERROR: output h5 file cannot be in your input h5 files.\n'
                 'Change your --h5_out and --h5_tables input.')
    elif h5_out.split('/')[-1] in [f.split('/')[-1] for f in glob(h5_out)]:
        os.system('rm {}'.format(h5_out))

    # If alternative solset number is given, we will make a temp h5 file that has the alternative solset number because the code runs on sol000 (will be cleaned up in the end)
    if use_solset != 'sol000':
        for h5_ind, h5 in enumerate(h5_tables):
            temph5 = h5.replace('.h5', '_temph5merger.h5')
            print('Using different solset. Make temporary h5 file: ' + temph5)
            os.system('cp ' + h5 + ' ' + temph5)
            _change_solset(temph5, use_solset, 'sol000')
            h5_tables[h5_ind] = temph5

    #################################################
    #################### MERGING ####################
    #################################################

    # Check if frequencies from h5_tables overlap
    if not freq_concat:
        check_freq_overlap(h5_tables)
    if not time_concat:
        check_time_overlap(h5_tables)
    if freq_concat and time_concat:
        sys.exit(
            "ERROR: Cannot do both time and frequency concat (ask Jurjen for assistance to implement this feature)")

    # Merge class setup
    merge = MergeH5(h5_out=h5_out, h5_tables=h5_tables, ms_files=ms_files, convert_tec=convert_tec,
                    merge_all_in_one=merge_all_in_one, h5_time_freq=h5_time_freq, filtered_dir=filtered_dir,
                    no_antenna_crash=no_antenna_crash, freq_concat=freq_concat, time_concat=time_concat)

    # Time averaging
    if time_av:
        print("Time averaging with factor " + str(time_av))
        if len(merge.ax_time) >= time_av:
            merge.ax_time = running_mean(merge.ax_time, time_av)[::int(time_av)]
        else:
            print("WARNING: Time averaging factor larger than time axis, no averaging applied")

    # Freq averaging
    if freq_av:
        print("Frequency averaging with factor " + str(freq_av))
        if len(merge.ax_freq) >= freq_av:
            merge.ax_freq = running_mean(merge.ax_freq, freq_av)[::int(freq_av)]
        else:
            print("WARNING: Frequency averaging factor larger than frequency axis, no averaging applied")

    # Get all keynames
    merge.get_allkeys()

    # Merging
    for st_group in merge.all_soltabs:
        if len(st_group) > 0:
            for st in st_group:
                merge.get_model_h5('sol000', st)
                merge.merge_tables('sol000', st, min_distance)
            if not merge.doublefulljones:
                if merge.convert_tec and (('phase' in st_group[0]) or ('tec' in st_group[0])):
                    # make sure tec is merged in phase only (if convert_tec==True)
                    merge.create_new_dataset('sol000', 'phase')
                else:
                    merge.create_new_dataset('sol000', st)
    if merge.doublefulljones:
        merge.matrix_multiplication()
        merge.create_new_dataset('sol000', 'phase')
        merge.create_new_dataset('sol000', 'amplitude')

    # Make sure direction tables are in the same format
    merge.format_tables()

    # Propagate weight flags from input into output
    if propagate_flags:
        merge.add_weights()

    # Add antennas
    if (add_cs or add_ms_stations) and len(merge.ms) == 0:
        sys.exit('ERROR: --add_cs and --add_ms_stations need an MS, given with --ms.')
    if add_cs:
        merge.add_ms_antennas(keep_h5_interstations=True)
    elif add_ms_stations:
        merge.add_ms_antennas(keep_h5_interstations=False)
    else:
        merge.add_h5_antennas()

    # Close all remaining open h5 files (should be unnecessary)
    tables.file._open_files.close_all()

    # If amplitude000 or phase000 are missing --> add a template for these
    #merge.add_template()

    # Add mock direction
    if add_directions:
        merge.add_empty_directions(add_directions)

    # Reorder directions for Python 2
    if sys.version_info.major == 2:
        merge.reorder_directions()

    # Check if station weights are fully flagged in input and flag in output as well
    if check_flagged_station and not propagate_flags:
        merge.flag_stations()

    # Check table source size
    merge.reduce_memory_source()

    # Change polarization axis if necessary
    if no_pol:
        print('Remove polarization')
        merge.change_pol(nopol=True)
    elif single_pol or merge.poldim == 0:
        print('Make a single polarization')
        merge.change_pol(single=True)

    #################################################
    ############ POLARIZATION CONVERSION ############
    #################################################

    if lin2circ and circ2lin:
        sys.exit('Both polarization conversions are given, please choose 1.')

    elif lin2circ or circ2lin:

        if lin2circ:
            h5_polchange = h5_out[0:-3] + '_circ.h5'
            print('\nPolarization will be converted from linear to circular')
        else:
            h5_polchange = h5_out[0:-3] + '_lin.h5'
            print('\nPolarization will be converted from circular to linear')

        # Polarization conversion
        Pol = PolChange(h5_in=h5_out, h5_out=h5_polchange)

        Pol.create_template('phase')
        if Pol.G.ndim > 1:
            Pol.create_template('amplitude')

        Pol.create_new_gain_table(lin2circ, circ2lin)
        Pol.add_antenna_source_tables()

        os.system('rm ' + h5_out + ' && cp ' + h5_polchange + ' ' + h5_out + ' && rm ' + h5_polchange)

    # Cleanup temp files
    for h5 in h5_tables:
        if '_temph5merger.h5' in h5:
            os.system('rm ' + h5)

    # Validate output
    if check_output:
        output_check(h5_out)

    print('\nSee output file --> ' + h5_out + '\n\n###################\nEND MERGE H5 TABLES \n###################\n')

    try:
        # Close all tables that are still open
        tables.file._open_files.close_all()
    except:
        pass

    # Check content of h5
    if output_summary:
        h5_check(h5_out)

    return


def parse_input():
    """
    Command line argument parser

    :return: arguments in correct format
    """

    parser = ArgumentParser(
        description='h5parm merger script for merging h5parms containing phase and/or amplitude solutions for calibrating LOFAR observations')
    parser.add_argument('-out', '--h5_out', type=str,
                        help='Solution file output name. This name cannot be in the list of input solution files.')
    parser.add_argument('-in', '--h5_tables', type=str, nargs='+',
                        help='Solution files to merge (can be both given as list with wildcard or string).',
                        required=True)
    parser.add_argument('-ms', '--ms', type=str,
                        help='Measurement Set input files (can be both given as list with wildcard or string).')
    parser.add_argument('--h5_time_freq', type=str,
                        help='Solution file to use time and frequency arrays from. This is useful if the input solution files do not have the preferred time/freq resolution. '
                             'Input can be boolean (true/fals or 1/0) to use all h5 files or input can be 1 specific h5 file (<H5.h5>).')
    parser.add_argument('--time_av', type=int, help='Time averaging factor.')
    parser.add_argument('--freq_av', type=int, help='Frequency averaging factor.')
    parser.add_argument('--keep_tec', action='store_true', help='Do not convert TEC to phase.')
    parser.add_argument('--merge_all_in_one', action='store_true', help='Merge all solutions into one direction.')
    parser.add_argument('--lin2circ', action='store_true', help='Transform linear polarization to circular.')
    parser.add_argument('--circ2lin', action='store_true', help='Transform circular polarization to linear.')
    parser.add_argument('--add_direction', default=None,
                        help='Add direction with amplitude 1 and phase 0 (example: --add_direction [0.73,0.12]).')
    parser.add_argument('--single_pol', action='store_true', default=None,
                        help='Return only a single polarization axis if both polarizations are the same.')
    parser.add_argument('--no_pol', action='store_true', default=None,
                        help='Remove polarization axis if both polarizations are the same.')
    # parser.add_argument('--combine_h5', action='store_true', default=None, help='Merge H5 files with different time axis into 1.')
    parser.add_argument('--usesolset', type=str, default='sol000',
                        help='Choose a solset to merge from your input solution files (only necessary if not sol000 is used).')
    parser.add_argument('--filter_directions', type=str, default=None,
                        help='Filter out a list of indexed directions from your output solution file. Only lists allowed (example: --filter_directions [2, 3]).')
    parser.add_argument('--add_cs', action='store_true', default=None,
                        help='Add core stations to antenna output from MS (needs --ms).')
    parser.add_argument('--add_ms_stations', action='store_true', default=None,
                        help='Use only antenna stations from measurement set (needs --ms). Note that this is different from --add_cs, as it does not keep the international stations if these are not in the MS.')
    parser.add_argument('--no_stationflag_check', action='store_true', default=None,
                        help='Do not flag complete station (for all directions) if entire station is flagged somewhere in input solution file.')
    parser.add_argument('--propagate_flags', action='store_true', default=None,
                        help='Interpolate weights and return in output file.')
    parser.add_argument('--no_antenna_crash', action='store_true', default=None,
                        help='Do not check if antennas are in h5.')
    parser.add_argument('--output_summary', action='store_true', default=None, help='Give output summary.')
    parser.add_argument('--check_output', action='store_true', default=None,
                        help='Check if the output has all the correct output information.')
    parser.add_argument('--merge_diff_freq', action='store_true', default=None,
                        help='Merging tables over different frequency bands --> same as old "freq_concat" setting')
    parser.add_argument('--time_concat', action='store_true', default=None,
                        help='Merging tables over different time slots (ensuring correct interpolation)')
    parser.add_argument('--freq_concat', action='store_true', default=None,
                        help='Merging tables over different frequency bands (ensuring correct interpolation) --> same as old "merge_diff_freq" setting')
    parser.add_argument('--min_distance', type=float,
                        help='(ONLY PYTHON 3) Minimal coordinate distance between sources for merging directions (in degrees). If smaller, directions will be merged together.',
                        default=0.)

    args = parser.parse_args()

    # check if solset name is accepted
    if 'sol' not in args.usesolset or sum([c.isdigit() for c in args.usesolset]) != 3:
        sys.exit(
            args.usesolse + ' not an accepted name. Only sol000, sol001, sol002, ... are accepted names for solsets.')

    if args.filter_directions:
        if (args.filter_directions.startswith("[") and args.filter_directions.endswith("]")):
            filtered_dir = args.filter_directions.replace(' ', '').replace('[', '').replace(']', '').split(',')
            for n, v in enumerate(filtered_dir):
                if not v.isdigit():
                    sys.exit('--filter_directions can only have integers in the list.')
                else:
                    filtered_dir[n] = int(v)
            args.filter_directions = filtered_dir
        else:
            sys.exit('--filter_directions given but no list format. Please pass a list to --filter_directions.')
    else:
        args.filter_directions = []

    # make sure h5 tables in right format
    if '[' in args.h5_tables[0]:
        args.h5_tables = args.h5_tables[0].replace('[', '').replace(']', '').replace(' ', '').split(',')
    elif ' ' in args.h5_tables:
        args.h5_tables = args.h5_tables[0].split()

    if type(args.h5_tables) == str:
        args.h5_tables = glob(args.h5_tables)
    elif type(args.h5_tables) == list and len(args.h5_tables) == 1:
        args.h5_tables = glob(args.h5_tables[0])
    elif type(args.h5_tables) == list:
        h5tablestemp = []
        for h5 in args.h5_tables:
            h5tablestemp += glob(h5)
        args.h5_tables = h5tablestemp

    if args.add_direction:
        args.add_direction = args.add_direction.replace('[', '').replace(']', '').split(',')
        args.add_direction = [float(args.add_direction[0]), float(args.add_direction[1])]
        if args.add_direction[0] > pi * 6 or args.add_direction[1] > pi * 6:
            sys.exit('ERROR: Please give --add_direction values in radian.')

    if args.h5_time_freq is not None:
        if args.h5_time_freq.lower() == 'true' or str(args.h5_time_freq) == '1':
            args.h5_time_freq = True
        elif args.h5_time_freq.lower() == 'false' or str(args.h5_time_freq) == '0':
            args.h5_time_freq = False
    else:
        args.h5_time_freq = False

    return args


def main():
    """
    Main function
    """
    if sys.version_info.major == 2:
        print('WARNING: This code is optimized for Python 3. Please switch to Python 3 if possible.')

    args = parse_input()

    if args.merge_diff_freq:
        print('WARNING: --merge_diff_freq given, please use --freq_concat.')

    merge_h5(h5_out=args.h5_out,
             h5_tables=args.h5_tables,
             ms_files=args.ms,
             h5_time_freq=args.h5_time_freq,
             convert_tec=not args.keep_tec,
             merge_all_in_one=args.merge_all_in_one,
             lin2circ=args.lin2circ,
             circ2lin=args.circ2lin,
             add_directions=args.add_direction,
             single_pol=args.single_pol,
             no_pol=args.no_pol,
             use_solset=args.usesolset,
             filtered_dir=args.filter_directions,
             add_cs=args.add_cs,
             add_ms_stations=args.add_ms_stations,
             check_output=args.check_output,
             time_av=args.time_av,
             freq_av=args.freq_av,
             check_flagged_station=not args.no_stationflag_check,
             propagate_flags=args.propagate_flags,
             freq_concat=args.merge_diff_freq or args.freq_concat,
             time_concat=args.time_concat,
             no_antenna_crash=args.no_antenna_crash,
             output_summary=args.output_summary,
             min_distance=args.min_distance)


if __name__ == '__main__':
    main()
