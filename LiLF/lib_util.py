import os, sys, glob
import socket
import tempfile

import numpy as np
import pyregion
from astropy.io import fits

# remove some annoying warnings from astropy
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)

# load here to be sure to have "Agg" at the beginning
import matplotlib as mpl
mpl.use("Agg")

from LiLF import lib_img
from LiLF.lib_log import logger

def create_extregion(ra, dec, extent, color='yellow'):
    """
    Parameters
    ----------
    ra
    dec
    extent
    color

    Returns
    -------
    DS9 region centered on ra, dec with radius = extent
    """

    regtext = ['# Region file format: DS9 version 4.1']
    regtext.append(
        f'global color={color} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')
    regtext.append('fk5')
    regtext.append('circle(' + str(ra) + ',' + str(dec) + f',{extent})')
    nline = '\n'
    target = f"{nline}{nline.join(regtext)}"

    return target


def getCalibratorProperties():
    """
    Return properties of known calibrators.
    The lists below (sorted in RA) are incomplete,
    and should be expanded to include all calibrators that could possibly be used.
    """

    calibratorRAs           = np.array([24.4220808, 85.6505746, 123.4001379, 202.784479167, 202.8569, 212.835495, 277.3824204, 299.8681525]) # in degrees
    calibratorDecs          = np.array([33.1597594, 49.8520094, 48.2173778,  30.509088,     25.1429,  52.202770,  48.7461556,  40.7339156])  # in degrees
    calibratorNames         = np.array(["3C48",     "3C147",    "3C196",     "3C286",       "3C287",  "3C295",    "3C380",     "CygA"])

    return calibratorRAs, calibratorDecs, calibratorNames


def distanceOnSphere(RAs1, Decs1, RAs2, Decs2, rad=False):
    """
    Return the distances on the sphere from the set of points '(RAs1, Decs1)' to the
    set of points '(RAs2, Decs2)' using the spherical law of cosines.

    Using 'numpy.clip(..., -1, 1)' is necessary to counteract the effect of numerical errors, that can sometimes
    incorrectly cause '...' to be slightly larger than 1 or slightly smaller than -1. This leads to NaNs in the arccosine.
    """
    if rad: # rad in rad out
        return np.radians(np.arccos(np.clip(
            np.sin(Decs1) * np.sin(Decs2) +
            np.cos(Decs1) * np.cos(Decs2) *
            np.cos(RAs1 - RAs2), -1, 1)))
    else: # deg in deg out
        return np.degrees(np.arccos(np.clip(
               np.sin(np.radians(Decs1)) * np.sin(np.radians(Decs2)) +
               np.cos(np.radians(Decs1)) * np.cos(np.radians(Decs2)) *
               np.cos(np.radians(RAs1 - RAs2)), -1, 1)))


def check_rm(regexp):
    """
    Check if file exists and remove it
    Handle reg exp of glob and spaces
    """
    filenames = regexp.split(' ')
    for filename in filenames:
        # glob is used to check if file exists
        for f in glob.glob(filename):
            os.system("rm -r " + f)


class Sol_iterator(object):
    """
    Iterator on a list that keeps on returing
    the last element when the list is over
    """

    def __init__(self, vals=[]):
        self.vals = vals
        self.pos = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.pos < len(self.vals):
            val = self.vals[self.pos]
            self.pos += 1
            return val
        else:
            return self.vals[-1]


def lofar_nu2num(nu):
    """
    Get LOFAR SB number from the freq
    """
    nu_clk = 200. # 160 or 200 MHz, clock freq
    # nyquist zone (1 for LBA, 2 for HBA low, 3 for HBA mid-high)
    if nu < 100:
        n = 1
    elif nu < 190:
        n = 2
    else:
        n = 3

    if nu_clk == 200:
        SBband = 195312.5/1e6
    elif nu_clk == 160:
        SBband = 156250.0/1e6

    return int(np.floor((1024./nu_clk) * (nu - (n-1) * nu_clk/2.)))

def run_losoto(s, c, h5s, parsets, plots_dir=None, h5_dir=None) -> object:
    """
    s : scheduler
    c : cycle name, e.g. "final" (or h5 output in a format filename.h5)
    h5s : lists of H5parm files or string of 1 h5parm
    parsets : lists of parsets to execute
    plots_dir : rename the "plots" dir to this name at the end
    h5_dir : dir where to move the new h5parm
    """

    logger.info("Running LoSoTo...")
    if c[-3:] == '.h5':
        h5out = c
    else:
        h5out = 'cal-'+c+'.h5'

    if type(h5s) is str: h5s = [h5s]

    # convert from killMS
    for i, h5 in enumerate(h5s):
        if h5[-3:] == 'npz':
            newh5 = h5.replace('.npz','.h5')
            s.add('killMS2H5parm.py -V --nofulljones %s %s ' % (newh5, h5), log='losoto-'+c+'.log', commandType="python")
            s.run(check = True)
            h5s[i] = newh5

    # concat/move
    if len(h5s) > 1:
        check_rm(h5out)
        s.add('H5parm_collector.py -V -s sol000 -o '+h5out+' '+' '.join(h5s), log='losoto-'+c+'.log', commandType="python")
        s.run(check = True)
    elif h5s[0] != h5out:
        os.system('cp -r %s %s' % (h5s[0], h5out) )

    if h5_dir:
        os.system(f'mv {h5out} {h5_dir}/{h5out}')
        h5out = f'{h5_dir}/{h5out}'

    check_rm('plots')
    os.makedirs('plots')

    for parset in parsets:
        logger.debug('-- executing '+parset+'...')
        s.add('losoto -V '+h5out+' '+parset, log='losoto-'+c+'.log', logAppend=True, commandType="python")
        s.run(check = True)

    if plots_dir is None:
        check_rm('plots-' + c)
        os.system('mv plots plots-' + c)
    else:
        if not os.path.exists(plots_dir): os.system('mkdir '+plots_dir)
        os.system('mv plots/* '+plots_dir)
        check_rm('plots')


def run_wsclean(s, logfile, MSs_files, do_predict=False, concat_mss=False, keep_concat=False, reuse_concat=False, use_shm=False, **kwargs):
    """
    s : scheduler
    concat_mss : try to concatenate mss files to speed up wsclean
    keep_concat : keep the concat MSs for re-use either by -cont or by reuse_concat=True
    reuse_concat : reuse concatenated MS previously kept with keep_concat=True
    temp_dir : if true try to store temp file in /dev/shm to speed up wsclean
    args : parameters for wsclean, "_" are replaced with "-", any parms=None is ignored.
           To pass a parameter with no values use e.g. " no_update_model_required='' "
    """

    # Check whether we can combine MS files in time, if some (or all) of them have the same antennas.
    # This can speed up WSClean significantly (or slow it down, depending on the number and size of MSs and the clean call).
    if concat_mss:
        if not 'cont' in kwargs.keys() and not reuse_concat:
            from LiLF import lib_ms
            from itertools import groupby

            keyfunct = lambda x: ' '.join(sorted(lib_ms.MS(x).getAntennas()))
            MSs_list = sorted(MSs_files.split(), key=keyfunct) # needs to be sorted
            groups = []
            for k, g in groupby(MSs_list, keyfunct):
                g = list(g)
                # reorder in time to prevent wsclean bug
                times = [lib_ms.MS(MS).getTimeRange()[0] for MS in g]
                g = [MS for _, MS in sorted(zip(times, g))]
                groups.append(g)
            logger.info(f"Found {len(groups)} groups of datasets with same antennas.")
            for i, group in enumerate(groups, start=1):
                antennas = ', '.join(lib_ms.MS(group[0]).getAntennas())
                logger.info(f"WSClean MS group {i}: {group}")
                logger.debug(f"List of antennas: {antennas}")

            MSs_files_clean = []
            for g, group in enumerate(groups):
                check_rm(f'wsclean_concat_{g}.MS')
                # simply make a symlink for groups of 1, faster
                if len(group) == 1:
                    os.system(f'ln -s {group[0]} wsclean_concat_{g}.MS') # TEST - symlink should be the quickest
                    # os.system(f'cp -r {group[0]} wsclean_concat_{g}.MS')
                else:
                    if 'data_column' in kwargs.keys():
                        data_column = kwargs['data_column']
                    else:
                        data_column = 'CORRECTED_DATA'
                    s.add(f'taql select UVW, FLAG_CATEGORY, WEIGHT, SIGMA, ANTENNA1, ANTENNA2, ARRAY_ID, DATA_DESC_ID, EXPOSURE, FEED1, FEED2, FIELD_ID, FLAG_ROW, INTERVAL, OBSERVATION_ID, PROCESSOR_ID, SCAN_NUMBER, STATE_ID, TIME, TIME_CENTROID, {data_column}, FLAG, WEIGHT_SPECTRUM from {group} giving wsclean_concat_{g}.MS as plain', log=logfile, commandType='general')
                    s.run(check=True)
                MSs_files_clean.append(f'wsclean_concat_{g}.MS')
        else:
            # continue clean
            MSs_files_clean = glob.glob('wsclean_concat_*.MS')
            logger.info(f'Continue clean on concat MSs {MSs_files_clean}')
        MSs_files_clean = ' '.join(MSs_files_clean)
    else:
        MSs_files_clean = MSs_files

    wsc_parms = []
    #reordering_processors = np.min([len(MSs_files_clean),s.max_proc])

    # basic parms
    wsc_parms.append( '-j '+str(s.max_proc)+' -reorder -parallel-reordering 4 -wgridder-accuracy 0.0001 ' )
    if 'use_idg' in kwargs.keys():
        if s.cluster == 'Hamburg_fat' and socket.gethostname() in ['node31', 'node32', 'node33', 'node34', 'node35']:
            wsc_parms.append( '-idg-mode hybrid' )
            wsc_parms.append( '-mem 10' )
        else:
            wsc_parms.append( '-idg-mode cpu' )
            
    # limit parallel gridding to max_proc
    if 'parallel_gridding' in kwargs.keys() and kwargs['parallel_gridding'] > s.max_proc:
            kwargs['parallel_gridding'] = s.max_proc

    tmp_dir = None
    if use_shm and os.access('/dev/shm/', os.W_OK) and 'temp_dir' not in kwargs:
        tmp_dir = tempfile.mkdtemp(dir='/dev/shm')
        wsc_parms.append(f'-temp-dir {tmp_dir}')
        wsc_parms.append('-mem 90')  # use 90% of memory
    elif s.cluster == 'Spider':
        wsc_parms.append('-temp-dir /tmp/')
    #elif s.cluster == 'Hamburg_fat' and not 'temp_dir' in list(kwargs.keys()):
    #    wsc_parms.append( '-temp-dir /localwork.ssd' )

    try:
        # user defined parms
        for parm, value in list(kwargs.items()):
            if value is None: continue
            if parm == 'baseline_averaging' and value == '':
                scale = float(kwargs['scale'].replace('arcsec','')) # arcsec
                value = 1.87e3*60000.*2.*np.pi/(24.*60.*60*np.max(kwargs['size'])) # the np.max() is OK with both float and arrays
                if value > 10: value=10
                if value < 1: continue
            if parm == 'cont': 
                parm = 'continue'
                value = ''
                # if continue, remove nans from previous models
                lib_img.Image(kwargs['name']).nantozeroModel()
            if parm == 'size' and type(value) is int: value = '%i %i' % (value, value)
            if parm == 'size' and type(value) is list: value = '%i %i' % (value[0], value[1])
            wsc_parms.append( '-%s %s' % (parm.replace('_','-'), str(value)) )
    
        # files
        wsc_parms.append( MSs_files_clean )
    
        # create command string
        command_string = 'wsclean '+' '.join(wsc_parms)
        s.add(command_string, log=logfile, commandType='wsclean')
        logger.info('Running WSClean...')
        s.run(check=True)
    
        # Predict in case update_model_required cannot be used
        if do_predict == True:
            if 'apply_facet_solutions' in kwargs.keys():
                raise NotImplementedError('do_predict in combination with apply_facet_solutions is not implemented.')
            wsc_parms = []
            # keep imagename and channel number
            for parm, value in list(kwargs.items()):
                if value is None: continue
                #if 'min' in parm or 'max' in parm or parm == 'name' or parm == 'channels_out':
                if parm == 'name' or parm == 'channels_out' or parm == 'wgridder_accuracy' or parm == 'shift':
                    wsc_parms.append( '-%s %s' % (parm.replace('_','-'), str(value)) )
    
            # files (the original, not the concatenated)
            wsc_parms.append( MSs_files )
            lib_img.Image(kwargs['name']).nantozeroModel() # If we have fully flagged channel, set to zero so we don't get error
    
            # Test without reorder as it apperas to be faster
            # wsc_parms.insert(0, ' -reorder -parallel-reordering 4 ')
            command_string = 'wsclean -predict -padding 1.8 ' \
                             '-j '+str(s.max_proc)+' '+' '.join(wsc_parms)
            s.add(command_string, log=logfile, commandType='wsclean')
            s.run(check=True)
    finally:
        if tmp_dir and os.path.exists(tmp_dir):
            logger.info(f"Deleting temporary directory {tmp_dir}")
            check_rm(tmp_dir)
    if not keep_concat:
        check_rm('wsclean_concat_*.MS')


class Region_helper():
    """
    Simple class to get the extent of a ds9 region file containing one or more circles or polygons.
    All properties are returned in degrees.

    Parameters
    ----------
    filename: str
        Path to ds9 region file.
    """
    def __init__(self, filename):
        self.filename = filename
        self.reg_list = pyregion.open(filename)
        min_ra, max_ra, min_dec, max_dec = [], [], [], []
        for r in self.reg_list:
            # TODO: if necessary, box, ellipse and polygon can be added.
            if r.name == 'circle':
                c = r.coord_list # c_ra, c_dec, radius
                # how much RA does the radius correspond to
                radius_ra = np.rad2deg(2*np.arcsin(np.sin(np.deg2rad(c[2])/2)/np.cos(np.deg2rad(c[1]))))
                min_ra.append(c[0] - radius_ra)
                max_ra.append(c[0] + radius_ra)
                min_dec.append(c[1] - c[2])
                max_dec.append(c[1] + c[2])
            elif r.name == 'polygon':
                c = np.array(r.coord_list) # ra_i, dec_i, ra_i+1, dec_i+1
                ra_mask = np.zeros(len(c), dtype=bool)
                ra_mask[::2] = True
                p_ra  = c[ra_mask]
                p_dec = c[~ra_mask]
                min_ra.append(np.min(p_ra))
                max_ra.append(np.max(p_ra))
                min_dec.append(np.min(p_dec))
                max_dec.append(np.max(p_dec))
            else:
                logger.error('Region type {} not supported.'.format(r.name))
                sys.exit(1)
        self.min_ra = np.min(min_ra)
        self.max_ra = np.max(max_ra)
        self.min_dec = np.min(min_dec)
        self.max_dec = np.max(max_dec)

    def get_center(self):
        """ Return center point [ra, dec] """
        return 0.5 * np.array([self.min_ra + self.max_ra, self.min_dec + self.max_dec])

    def get_width(self):
        """ Return RA width in degree (at center declination)"""
        delta_ra = self.max_ra - self.min_ra
        width = 2*np.arcsin(np.cos(np.deg2rad(self.get_center()[1]))*np.sin(np.deg2rad(delta_ra/2)))
        width = np.rad2deg(width)
        return width

    def get_height(self):
        """ Return height in degree"""
        return self.max_dec - self.min_dec

    def __len__(self):
        return len(self.reg_list)

def get_template_image(reference_ra_deg, reference_dec_deg, ximsize=512, yimsize=512, cellsize_deg=0.000417, fill_val=0):
    """
    Make a blank image and return
    adapted from https://github.com/darafferty/LSMTool/blob/master/lsmtool/operations_lib.py#L619

    Parameters
    ----------
    reference_ra_deg : float
        RA for center of output image
    reference_dec_deg : float
        Dec for center of output image
    ximsize : int, optional
        Size of output image
    yimsize : int, optional
        Size of output image
    cellsize_deg : float, optional
        Size of a pixel in degrees
    fill_val : int, optional
        Value with which to fill the image
    """

    # Make fits hdu
    # Axis order is [STOKES, FREQ, DEC, RA]
    shape_out = [yimsize, ximsize]
    hdu = fits.PrimaryHDU(np.ones(shape_out, dtype=np.float32)*fill_val)
    hdulist = fits.HDUList([hdu])
    header = hdulist[0].header

    # Add RA, Dec info
    i = 1
    header['CRVAL{}'.format(i)] = reference_ra_deg
    header['CDELT{}'.format(i)] = -cellsize_deg
    header['CRPIX{}'.format(i)] = ximsize / 2.0
    header['CUNIT{}'.format(i)] = 'deg'
    header['CTYPE{}'.format(i)] = 'RA---SIN'
    i += 1
    header['CRVAL{}'.format(i)] = reference_dec_deg
    header['CDELT{}'.format(i)] = cellsize_deg
    header['CRPIX{}'.format(i)] = yimsize / 2.0
    header['CUNIT{}'.format(i)] = 'deg'
    header['CTYPE{}'.format(i)] = 'DEC--SIN'

    # Add equinox
    header['EQUINOX'] = 2000.0

    # Add telescope
    header['TELESCOP'] = 'LOFAR'

    hdulist[0].header = header
    return hdulist
