import os, sys, glob
import socket
import datetime

from casacore import tables
import lsmtool.skymodel
import numpy as np
import multiprocessing, subprocess
from threading import Thread
from queue import Queue
import pyregion
from astropy.io import fits
import gc

# remove some annoying warnings from astropy
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)

if (sys.version_info > (3, 0)):
    from configparser import ConfigParser
else:
    from ConfigParser import ConfigParser

# load here to be sure to have "Agg" at the beginning
import matplotlib as mpl
mpl.use("Agg")

from LiLF import lib_img
from LiLF.lib_log import logger

def getParset(parsetFile=''):
    """
    Get parset file and return dict of values
    """
    def add_default(section, option, val):
        if not config.has_option(section, option): config.set(section, option, val)
    
    if parsetFile == '' and os.path.exists('lilf.config'): parsetFile='lilf.config'
    if parsetFile == '' and os.path.exists('../lilf.config'): parsetFile='../lilf.config'

    config = ConfigParser(defaults=None)
    config.read(parsetFile)
    
    # add pipeline sections and defaul parset dir:
    for pipeline in glob.glob(os.path.dirname(__file__)+'/../parsets/*'):
        pipeline = os.path.basename(pipeline)
        if not config.has_section(pipeline): config.add_section(pipeline)
        if not config.has_option(pipeline, 'parset_dir'):
                config.set(pipeline, 'parset_dir', os.path.dirname(__file__)+'/../parsets/'+pipeline)
    # add other sections
    if not config.has_section('flag'): config.add_section('flag')
    if not config.has_section('model'): config.add_section('model')
    if not config.has_section('PiLL'): config.add_section('PiLL')

    ### LOFAR ###

    # PiLL
    add_default('PiLL', 'working_dir', os.getcwd())
    add_default('PiLL', 'redo_cal', 'False') # re-do the calibrator although it is in the archive
    add_default('PiLL', 'download_file', '') # html.txt file to use instead of staging
    add_default('PiLL', 'project', '')
    add_default('PiLL', 'target', '')
    add_default('PiLL', 'obsid', '') # unique ID
    # preprocess
    add_default('LOFAR_preprocess', 'fix_table', 'True') # fix bug in some old observations
    add_default('LOFAR_preprocess', 'renameavg', 'True')
    add_default('LOFAR_preprocess', 'keep_IS', 'False')
    add_default('LOFAR_preprocess', 'backup_full_res', 'False')
    add_default('LOFAR_preprocess', 'demix_sources', '')  # Demix  sources in these patches (e.g. [VirA,TauA], default: No demix
    add_default('LOFAR_preprocess', 'demix_skymodel', '')  # Use non-default demix skymodel.
    add_default('LOFAR_preprocess', 'demix_field_skymodel', 'gsm')  # Provide a custom target skymodel instead of online gsm model. Set to '' to ignore target.
    add_default('LOFAR_preprocess', 'run_aoflagger', 'False')  # run aoflagger on individual sub-bands, only in cases where this was not one by the observatory!
    add_default('LOFAR_preprocess', 'tar', 'False')  # Tar MS files at the end 
    # cal
    add_default('LOFAR_cal', 'data_dir', 'data-bkp/')
    add_default('LOFAR_cal', 'skymodel', '') # by default use calib-simple.skydb for LBA and calib-hba.skydb for HBA
    add_default('LOFAR_cal', 'imaging', 'False')
    add_default('LOFAR_cal', 'fillmissingedges', 'True')
    add_default('LOFAR_cal', 'sparse_sb', 'False') # change flagging so that we can handle data with alternating SBs only
    add_default('LOFAR_cal', 'develop', 'False') # if true prevents the deletion of files
    # timesplit
    add_default('LOFAR_timesplit', 'data_dir', 'data-bkp/')
    add_default('LOFAR_timesplit', 'cal_dir', '') # by default the repository is tested, otherwise ../obsid_3[c|C]*
    add_default('LOFAR_timesplit', 'ngroups', '1')
    add_default('LOFAR_timesplit', 'initc', '0')
    add_default('LOFAR_timesplit', 'no_aoflagger', 'False') # TEST: Skip aoflagger (e.g. for observations of A-Team sources)
    # self
    add_default('LOFAR_ddparallel', 'maxIter', '2')
    add_default('LOFAR_ddparallel', 'subfield', '') # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
    add_default('LOFAR_ddparallel', 'subfield_min_flux', '20') # min flux within calibration subfield
    add_default('LOFAR_ddparallel', 'ph_sol_mode', 'phase') # phase or tecandphase
    add_default('LOFAR_ddparallel', 'intrinsic', 'True')
    # dd
    add_default('LOFAR_ddserial', 'maxIter', '2')
    add_default('LOFAR_ddserial', 'minCalFlux60', '0.7')
    add_default('LOFAR_ddserial', 'solve_amp', 'True') # to disable amp sols
    # add_default('LOFAR_ddserial', 'removeExtendedCutoff', '0.0005')
    add_default('LOFAR_ddserial', 'target_dir', '') # ra,dec
    add_default('LOFAR_ddserial', 'manual_ddserial_cal', '')
    # add_default('LOFAR_ddserial', 'solve_tec', 'False') # per default, solve each dd for scalarphase. if solve_tec==True, solve for TEC instead.
    # extract
    add_default('LOFAR_extract', 'max_niter', '10')
    add_default('LOFAR_extract', 'subtract_region', '') # Sources inside extract-reg that should still be subtracted! Use this e.g. for individual problematic sources in a large extractReg
    add_default('LOFAR_extract', 'ph_sol_mode', 'phase') # tecandphase, phase
    add_default('LOFAR_extract', 'amp_sol_mode', 'diagonal') # diagonal, fulljones
    add_default('LOFAR_extract', 'beam_cut', '0.3') # up to which distance a pointing will be considered
    add_default('LOFAR_extract', 'no_selfcal', 'False') # just extract the data, do not perform selfcal - use this if u want to use e.g. Reinout van Weeren's facet_seflcal script
    add_default('LOFAR_extract', 'ampcal', 'auto')
    add_default('LOFAR_extract', 'extractRegion', 'target.reg')
    # quality
    add_default('LOFAR_quality', 'self_dir', 'self')
    add_default('LOFAR_quality', 'ddcal_dir', 'ddcal')
    # virgo
    add_default('LOFAR_virgo', 'cal_dir', '')
    add_default('LOFAR_virgo', 'data_dir', './')
    # m87
    add_default('LOFAR_m87', 'data_dir', './')
    add_default('LOFAR_m87', 'updateweights', 'False')
    add_default('LOFAR_m87', 'skipmodel', 'False')
    add_default('LOFAR_m87', 'model_dir', '')
    # peel
    ### uGMRT ###
    # init - deprecated
    #add_default('uGMRT_init', 'data_dir', './datadir')
    # cal - deprecated
    #add_default('uGMRT_cal', 'skymodel', os.path.dirname(__file__)+'/../models/calib-simple.skydb')

    ### General ###

    # flag
    add_default('flag', 'stations', '') # LOFAR
    add_default('flag', 'antennas', '') # uGMRT
    # model
    add_default('model', 'sourcedb', '')
    add_default('model', 'fits_model', '')
    add_default('model', 'apparent', 'False')
    add_default('model', 'userReg', '')


    return config

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


def columnAddSimilar(pathMS, columnNameNew, columnNameSimilar, dataManagerInfoNameNew, overwrite = False, fillWithOnes = True, comment = "", verbose = False):
    # more to lib_ms
    """
    Add a column to a MS that is similar to a pre-existing column (in shape, but not in values).
    pathMS:                 path of the MS
    columnNameNew:          name of the column to be added
    columnNameSimilar:      name of the column from which properties are copied (e.g. "DATA")
    dataManagerInfoNameNew: string value for the data manager info (DMI) keyword "NAME" (should be unique in the MS)
    overwrite:              whether or not to overwrite column 'columnNameNew' if it already exists
    fillWithOnes:           whether or not to fill the newly-made column with ones
    verbose:                whether or not to produce abundant output
    """
    t = tables.table(pathMS, readonly = False)

    if (columnExists(t, columnNameNew) and not overwrite):
        logger.warning("Attempt to add column '" + columnNameNew + "' aborted, as it already exists and 'overwrite = False' in columnAddSimilar(...).")
    else: # Either the column does not exist yet, or it does but overwriting is allowed.

        # Remove column if necessary.
        if (columnExists(t, columnNameNew)):
            logger.info("Removing column '" + columnNameNew + "'...")
            t.removecols(columnNameNew)

        # Add column.
        columnDescription       = t.getcoldesc(columnNameSimilar)
        dataManagerInfo         = t.getdminfo(columnNameSimilar)

        if (verbose):
            logger.debug("columnDescription:")
            logger.debug(columnDescription)
            logger.debug("dataManagerInfo:")
            logger.debug(dataManagerInfo)

        columnDescription["comment"] = ""
        # What about adding something here like:
        #columnDescription["dataManagerGroup"] = ...?
        dataManagerInfo["NAME"]      = dataManagerInfoNameNew

        if (verbose):
            logger.debug("columnDescription (updated):")
            logger.debug(columnDescription)
            logger.debug("dataManagerInfo (updated):")
            logger.debug(dataManagerInfo)

        logger.info("Adding column '" + columnNameNew + "'...")
        t.addcols(tables.makecoldesc(columnNameNew, columnDescription), dataManagerInfo)

        # Fill with ones if desired.
        if (fillWithOnes):
            logger.info("Filling column '" + columnNameNew + "' with ones...")
            columnDataSimilar = t.getcol(columnNameSimilar)
            t.putcol(columnNameNew, np.ones_like(columnDataSimilar))

    # Close the table to avoid that it is locked for further use.
    t.close()


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
    if nu < 90:
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
            s.add('killMS2H5parm.py -V --nofulljones %s %s ' % (newh5, h5), log='losoto-'+c+'.log', commandType="python", processors='max')
            s.run(check = True)
            h5s[i] = newh5

    # concat/move
    if len(h5s) > 1:
        check_rm(h5out)
        s.add('H5parm_collector.py -V -s sol000 -o '+h5out+' '+' '.join(h5s), log='losoto-'+c+'.log', commandType="python", processors='max')
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
        s.add('losoto -V '+h5out+' '+parset, log='losoto-'+c+'.log', logAppend=True, commandType="python", processors='max')
        s.run(check = True)

    if plots_dir is None:
        check_rm('plots-' + c)
        os.system('mv plots plots-' + c)
    else:
        if not os.path.exists(plots_dir): os.system('mkdir '+plots_dir)
        os.system('mv plots/* '+plots_dir)
        check_rm('plots')


def run_wsclean(s, logfile, MSs_files, do_predict=False, concat_mss=False, keep_concat=False, **kwargs):
    """
    s : scheduler
    concat_mss : try to concatenate mss files to speed up wsclean
    args : parameters for wsclean, "_" are replaced with "-", any parms=None is ignored.
           To pass a parameter with no values use e.g. " no_update_model_required='' "
    """

    # Check whether we can combine MS files in time, if some (or all) of them have the same antennas.
    # This speeds up WSClean significantly.
    if concat_mss:
        if not 'cont' in kwargs.keys():
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
                # cp is faster than taql, ok for groups of 1
                if len(group) == 1:
                    os.system(f'cp -r {group[0]} wsclean_concat_{g}.MS')
                else:
                    s.add(f'taql select from {group} giving wsclean_concat_{g}.MS as plain', log=logfile, commandType='general')
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
    reordering_processors = np.min([len(MSs_files_clean),s.max_processors])

    # basic parms
    wsc_parms.append( '-j '+str(s.max_processors)+' -reorder -parallel-reordering 4 ' )
    if 'use_idg' in kwargs.keys():
        if s.get_cluster() == 'Hamburg_fat' and socket.gethostname() in ['node31', 'node32', 'node33', 'node34', 'node35']:
            wsc_parms.append( '-idg-mode hybrid' )
            wsc_parms.append( '-mem 10' )
        else:
            wsc_parms.append( '-idg-mode cpu' )

    # temp dir
    #if s.get_cluster() == 'Hamburg_fat' and not 'temp_dir' in list(kwargs.keys()):
    #    wsc_parms.append( '-temp-dir /localwork.ssd' )
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
    s.add(command_string, log=logfile, commandType='wsclean', processors='max')
    logger.info('Running WSClean...')
    s.run(check=True)

    # Predict in case update_model_required cannot be used
    if do_predict == True:
        wsc_parms = []
        # keep imagename and channel number
        for parm, value in list(kwargs.items()):
            if value is None: continue
            #if 'min' in parm or 'max' in parm or parm == 'name' or parm == 'channels_out':
            if parm == 'name' or parm == 'channels_out' or parm == 'use_wgridder' or parm == 'wgridder_accuracy':
                wsc_parms.append( '-%s %s' % (parm.replace('_','-'), str(value)) )

        # files (the original, not the concatenated)
        wsc_parms.append( MSs_files )
        lib_img.Image(kwargs['name']).nantozeroModel() # If we have fully flagged channel, set to zero so we don't get error

        # Test without reorder as it apperas to be faster
        # wsc_parms.insert(0, ' -reorder -parallel-reordering 4 ')
        command_string = 'wsclean -predict -padding 1.8 ' \
                         '-j '+str(s.max_processors)+' '+' '.join(wsc_parms)
        s.add(command_string, log=logfile, commandType='wsclean', processors='max')
        s.run(check=True)
    if not keep_concat:
        check_rm('wsclean_concat_*.MS')

def run_DDF(s, logfile, **kwargs):
    """
    s : scheduler
    args : parameters for ddfacet, "_" are replaced with "-", any parms=None is ignored.
           To pass a parameter with no values use e.g. " no_update_model_required='' "
    """
    
    ddf_parms = []

    # basic parms
    ddf_parms.append( '--Log-Boring 1 --Debug-Pdb never --Parallel-NCPU %i --Misc-IgnoreDeprecationMarking=1 ' % (s.max_processors) )

    # cache dir
    if not 'Cache_Dir' in list(kwargs.keys()):
        ddf_parms.append( '--Cache-Dir .' )

    # user defined parms
    for parm, value in list(kwargs.items()):
        if value is None: continue
        if isinstance(value, str):
            if '$' in value: # escape dollar signs (e.g. of BeamFits)
                value = "'" + value + "'"
        ddf_parms.append( '--%s=%s' % (parm.replace('_','-'), str(value)) )

    # files
    #wsc_parms.append( MSs_files )

    # create command string
    command_string = 'DDF.py '+' '.join(ddf_parms)
    s.add(command_string, log=logfile, commandType='DDFacet', processors='max')
    s.run(check=True)


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

class Skip(Exception):
    pass

class Exit(Exception):
    pass

class Walker():
    """
    An object of this class may be used to re-run a pipeline without repeating steps that were completed previously.
    Use like:
    w = Walker("filename.walker")
    with w.if_todo("stepname"):
        Do whatever...

    Adopted from https://stackoverflow.com/questions/12594148/skipping-execution-of-with-block
    """
    def __init__(self, filename):
        open(filename, 'a').close() # create the file if doesn't exists
        self.filename = os.path.abspath(filename)
        self.__skip__ = False
        self.__step__ = None
        self.__inittime__ = None
        self.__globaltimeinit__ = datetime.datetime.now()

    def if_todo(self, stepname):
        """
        This is basically a way to get a context manager to accept an argument. Will return "self" as context manager
        if called as context manager.
        """
        self.__skip__ = False
        self.__step__ = stepname
        with open(self.filename, "r") as f:
            for stepname_done in f:
                if stepname == stepname_done.split('#')[0].rstrip():
                    self.__skip__ = True
        return self

    def __enter__(self):
        """
        Skips body of with-statement if __skip__.
        This uses some kind of dirty hack that might only work in CPython.
        """
        if self.__skip__:
            sys.settrace(lambda *args, **keys: None)
            frame = sys._getframe(1)
            frame.f_trace = self.trace
        else:
            logger.log(20, '>> start >> {}'.format(self.__step__))
            self.__timeinit__ = datetime.datetime.now()

    def trace(self, frame, event, arg):
        raise Skip()

    def __exit__(self, type, value, traceback):
        """
        Catch "Skip" errors, if not skipped, write to file after exited without exceptions.
        """
        if type is None:
            with open(self.filename, "a") as f:
                delta = 'h '.join(str(datetime.datetime.now() - self.__timeinit__).split(':')[:-1])+'m'
                f.write(self.__step__ + ' # '+delta+' ' +'\n')
            logger.info('<< done << {}'.format(self.__step__))
            return  # No exception
        if issubclass(type, Skip):
            logger.warning('>> skip << {}'.format(self.__step__))
            return True  # Suppress special SkipWithBlock exception
        if issubclass(type, Exit):
            logger.error('<< exit << {}'.format(self.__step__))
            return True
        
    def alldone(self):
        delta = 'h '.join(str(datetime.datetime.now() - self.__globaltimeinit__).split(':')[:-1])+'m'
        logger.info('Done. Total time: '+delta)

class Scheduler():
    def __init__(self, qsub = None, maxThreads = None, max_processors = None, log_dir = 'logs', dry = False):
        """
        qsub:           if true call a shell script which call qsub and then wait
                        for the process to finish before returning
        maxThreads:    max number of parallel processes
        dry:            don't schedule job
        max_processors: max number of processors in a node (ignored if qsub=False)
        """
        self.hostname = socket.gethostname()
        self.cluster = self.get_cluster()
        self.log_dir = log_dir
        self.qsub    = qsub
        # if qsub/max_thread/max_processors not set, guess from the cluster
        # if they are set, double check number are reasonable
        if (self.qsub == None):
            if (self.cluster == "Hamburg"):
                self.qsub = True
            else:
                self.qsub = False
        else:
            if ((self.qsub is False and self.cluster == "Hamburg") or
               (self.qsub is True and (self.cluster == "Leiden" or self.cluster == "CEP3" or
                                       self.cluster == "Hamburg_fat" or self.cluster == "Pleiadi" or self.cluster == "Herts"))):
                logger.critical('Qsub set to %s and cluster is %s.' % (str(qsub), self.cluster))
                sys.exit(1)

        if (maxThreads is None):
            if (self.cluster == "Hamburg"):
                self.maxThreads = 32
            else:
                self.maxThreads = multiprocessing.cpu_count()
        else:
            self.maxThreads = maxThreads

        if (max_processors == None):
            if   (self.cluster == "Hamburg"):
                self.max_processors = 6
            else:
                self.max_processors = multiprocessing.cpu_count()
        else:
            self.max_processors = max_processors

        self.dry = dry

        logger.info("Scheduler initialised for cluster " + self.cluster + ": " + self.hostname + " (maxThreads: " + str(self.maxThreads) + ", qsub (multinode): " +
                     str(self.qsub) + ", max_processors: " + str(self.max_processors) + ").")


        self.action_list = []
        self.log_list    = []  # list of 2-tuples of the type: (log filename, type of action)


    def get_cluster(self):
        """
        Find in which computing cluster the pipeline is running
        """
        hostname = self.hostname
        if (hostname == 'lgc1' or hostname == 'lgc2'):
            return "Hamburg"
        elif ('r' == hostname[0] and 'c' == hostname[3] and 's' == hostname[6]):
            return "Pleiadi"
        elif ('node3' in hostname):
            return "Hamburg_fat"
        elif ('node' in hostname):
            return "Herts"
        elif ('leidenuniv' in hostname):
            return "Leiden"
        elif (hostname[0 : 3] == 'lof'):
            return "CEP3"
        else:
            logger.warning('Hostname %s unknown.' % hostname)
            return "Unknown"


    def add(self, cmd = '', log = '', logAppend = True, commandType = '', processors = None):
        """
        Add a command to the scheduler list
        cmd:         the command to run
        log:         log file name that can be checked at the end
        logAppend:  if True append, otherwise replace
        commandType: can be a list of known command types as "wsclean", "DP3", ...
        processors:  number of processors to use, can be "max" to automatically use max number of processors per node
        """

        if (log != ''):
            log = self.log_dir + '/' + log

            if (logAppend):
                cmd += " >> "
            else:
                cmd += " > "
            cmd += log + " 2>&1"

        # if running wsclean add the string
        if commandType == 'wsclean':
            logger.debug('Running wsclean: %s' % cmd)
        elif commandType == 'DP3':
            logger.debug('Running DP3: %s' % cmd)
        #elif commandType == 'singularity':
        #    cmd = 'SINGULARITY_TMPDIR=/dev/shm singularity exec -B /tmp,/dev/shm,/localwork,/localwork.ssd,/home /home/fdg/node31/opt/src/lofar_sksp_ddf.simg ' + cmd
        #    logger.debug('Running singularity: %s' % cmd)
        elif (commandType.lower() == "ddfacet" or commandType.lower() == 'ddf'):
            logger.debug('Running DDFacet: %s' % cmd)
        elif commandType == 'python':
            logger.debug('Running python: %s' % cmd)
        else:
            logger.debug('Running general: %s' % cmd)


        if (processors != None and processors == 'max'):
            processors = self.max_processors

        if self.qsub:
            # if number of processors not specified, try to find automatically
            if (processors == None):
                processors = 1 # default use single CPU
                if ("DP3" == cmd[ : 4]):
                    processors = 1
                if ("wsclean" == cmd[ : 7]):
                    processors = self.max_processors
            if (processors > self.max_processors):
                processors = self.max_processors

            self.action_list.append([str(processors), '\'' + cmd + '\''])
        else:
            self.action_list.append(cmd)

        if (log != ""):
            self.log_list.append((log, commandType))


    def run(self, check = False, maxThreads = None):
        """
        If 'check' is True, a check is done on every log in 'self.log_list'.
        If max_thread != None, then it overrides the global values, useful for special commands that need a lower number of threads.
        """

        def worker(queue):
            for cmd in iter(queue.get, None):
                if self.qsub and self.cluster == "Hamburg":
                    cmd = 'salloc --job-name LBApipe --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                            ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                gc.collect()
                subprocess.call(cmd, shell = True)

        # limit threads only when qsub doesn't do it
        if (maxThreads == None):
            maxThreads_run = self.maxThreads
        else:
            maxThreads_run = min(maxThreads, self.maxThreads)

        q       = Queue()
        threads = [Thread(target = worker, args=(q,)) for _ in range(maxThreads_run)]

        for i, t in enumerate(threads): # start workers
            t.daemon = True
            t.start()

        for action in self.action_list:
            if (self.dry):
                continue # don't schedule if dry run
            q.put_nowait(action)
        for _ in threads:
            q.put(None) # signal no more commands
        for t in threads:
            t.join()

        # check outcomes on logs
        if (check):
            for log, commandType in self.log_list:
                self.check_run(log, commandType)

        # reset list of commands
        self.action_list = []
        self.log_list    = []


    def check_run(self, log = "", commandType = ""):
        """
        Produce a warning if a command didn't close the log properly i.e. it crashed
        NOTE: grep, -L inverse match, -l return only filename
        """

        if (not os.path.exists(log)):
            logger.warning("No log file found to check results: " + log)
            return 1

        if (commandType == "DP3"):
            out = subprocess.check_output(r'grep -L "Finishing processing" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Segmentation fault\|Killed" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Aborted (core dumped)" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -i -l "Exception" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # this interferes with the missingantennabehaviour=error option...
            # out += subprocess.check_output('grep -l "error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "misspelled" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "CASA"):
            out = subprocess.check_output(r'grep -l "[a-z]Error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "An error occurred running" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "\*\*\* Error \*\*\*" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "wsclean"):
            out = subprocess.check_output(r'grep -l "exception occur" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Segmentation fault\|Killed" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Aborted" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # out += subprocess.check_output('grep -L "Cleaning up temporary files..." '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType.lower() == "ddfacet" or commandType.lower() == 'ddf'):
            out = subprocess.check_output(r'grep -l "Traceback (most recent call last):" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "exception occur" ' + log + ' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "raise Exception" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Segmentation fault\|Killed" ' + log + ' ; exit 0', shell=True,
                                           stderr=subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "killed by signal" ' + log + ' ; exit 0', shell=True,
                                           stderr=subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Aborted" ' + log + ' ; exit 0', shell=True, stderr=subprocess.STDOUT)

        elif (commandType == "python"):
            out = subprocess.check_output(r'grep -l "Traceback (most recent call last):" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Segmentation fault\|Killed" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # out += subprocess.check_output(r'grep -i -l \'(?=^((?!error000).)*$).*Error.*\' '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -i -l "Critical" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "ERROR" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Traceback (most recent call last)" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

#        elif (commandType == "singularity"):
#            out = subprocess.check_output('grep -l "Traceback (most recent call last):" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
#            out += subprocess.check_output('grep -i -l \'(?=^((?!error000).)*$).*Error.*\' '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
#            out += subprocess.check_output('grep -i -l "Critical" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "general"):
            out = subprocess.check_output('grep -l -i "error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        else:
            logger.warning("Unknown command type for log checking: '" + commandType + "'")
            return 1

        if out != b'':
            out = out.split(b'\n')[0].decode()
            logger.error(commandType+' run problem on:\n'+out)
            raise RuntimeError(commandType+' run problem on:\n'+out)

        return 0

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
