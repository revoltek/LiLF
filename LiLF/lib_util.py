import os, sys, glob
import socket
import datetime

import numpy as np
import multiprocessing, subprocess
from threading import Thread
from queue import Queue
import pyregion
from astropy.io import fits
import gc
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, as_completed
import subprocess, os, time, shlex, signal, gc

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
    #reordering_processors = np.min([len(MSs_files_clean),s.maxProcs])

    # basic parms
    wsc_parms.append( '-j ' + str(s.maxWorkers) + ' -reorder -parallel-reordering 4 ')
    if 'use_idg' in kwargs.keys():
        if s.cluster == 'Hamburg_fat' and socket.gethostname() in ['node31', 'node32', 'node33', 'node34', 'node35']:
            wsc_parms.append( '-idg-mode hybrid' )
            wsc_parms.append( '-mem 10' )
        else:
            wsc_parms.append( '-idg-mode cpu' )
            
    # limit parallel gridding to maxProcs
    if 'parallel_gridding' in kwargs.keys() and kwargs['parallel_gridding'] > s.maxWorkers:
            kwargs['parallel_gridding'] = s.maxWorkers

    # set the tmp dir to speed up
    if use_shm and os.access('/dev/shm/', os.W_OK) and not 'temp_dir' in list(kwargs.keys()):
        check_rm('/dev/shm/*') # remove possible leftovers
        wsc_parms.append( '-temp-dir /dev/shm/' )
        wsc_parms.append( '-mem 90' ) # use 90% of memory
    elif s.cluster == 'Spider':
        wsc_parms.append( '-temp-dir /tmp/' )
    #elif s.cluster == 'Hamburg_fat' and not 'temp_dir' in list(kwargs.keys()):
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
                         '-j ' + str(s.maxWorkers) + ' ' + ' '.join(wsc_parms)
        s.add(command_string, log=logfile, commandType='wsclean')
        s.run(check=True)
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


def _run_cmd(cmd):
    try:
        os.system(cmd)
        return 0, 0.0
    except Exception as e:
        logger.error(f"Error occurred while running command: {cmd}\n{e}")
        return 1, 0.0

def _run_cmd_gpt(cmd, log_path=None, timeout=None):
    # ChatGPT
    """
    Run a shell command, tee to log file, return (returncode, walltime_s).
    """
    t0 = time.time()
    # Ensure parent dirs for logs exist
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)

    # Open log and stream stdout/stderr
    with open(log_path, "a" if log_path else os.devnull) as logf:
        # Start in a process group so we can kill children if timeout
        with subprocess.Popen(cmd, shell=True, stdout=logf, stderr=logf,
                              preexec_fn=os.setsid) as p:
            try:
                p.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                # Kill the whole process group
                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                return (124, time.time() - t0)   # 124 like GNU timeout
    return (p.returncode, time.time() - t0)

def get_slurm_max_walltime():
    # ChatGPT
    import os, re, subprocess
    def _parse_slurm_time(tstr):
        if not tstr or tstr.upper() == "UNLIMITED":
            return None
        if '-' in tstr:
            days, hms = tstr.split('-', 1)
            days = int(days)
        else:
            days = 0
            hms = tstr
        hh, mm, ss = (list(map(int, hms.split(':'))) + [0, 0, 0])[:3]
        return days*86400 + hh*3600 + mm*60 + ss
    
    partition = os.getenv('SLURM_PARTITION')
    try:
        if partition:
            outp = subprocess.check_output(['scontrol', 'show', 'partition', partition], text=True, stderr=subprocess.DEVNULL)
            m = re.search(r'MaxTime=([^\s]+)', outp)
            if m:
                timestr = m.group(1)
                return _parse_slurm_time(timestr), timestr
        # try sinfo: list partitions and limits, take first non-empty limit
        out = subprocess.check_output(['sinfo', '-h', '-o', '%P %l'], text=True, stderr=subprocess.DEVNULL)
        for line in out.splitlines():
            part, limit = line.split(None, 1)
            if limit and limit != "UNKNOWN":
                timestr = limit.strip()
                return _parse_slurm_time(timestr), timestr
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    return None, None


class Scheduler():
    def __init__(self, backend='slurm', slurm_max_jobs=244, max_cpus_per_node=None, slurm_max_walltime=None, slurm_mem_per_cpu='8GB',
                 log_dir = 'logs', dry = False, container_path=None):
        """
        TODO max walltime
        TODO max_jobs autoset?
        backend: string, backend used to launch jobs. 'local' and 'slurm' are supported
        max_jobs:       max number of parallel processes (either on local node or on slurm cluster)
        max_cpucores:   max number of cpu cores usable in a node
        dry:            don't schedule job
        """
        self.backend = backend.lower()
        self.hostname = socket.gethostname()
        self.cluster = self.get_cluster()
        self.log_dir = log_dir

        # automatically set max cpucores if not set manually
        if max_cpus_per_node:
            self.max_cpus_per_node = int(max_cpus_per_node)
        elif max_cpus_per_node is None:
            if backend == 'local':
                self.max_cpus_per_node = multiprocessing.cpu_count()
            elif backend == 'slurm':
                slurm_cpus = os.getenv('SLURM_CPUS_ON_NODE')
                if slurm_cpus:
                    self.max_cpus_per_node = int(slurm_cpus)
                else:
                    logger.warning('Neither max_cpucores_per_node nor $SLURM_CPUS_ON_NODE defined - guessing cpus per node.')
                    self.max_cpus_per_node = multiprocessing.cpu_count()
             
        # HARDCODED LIMIT FOR NOW            
        self.max_cpucores_per_node = 16
        
        # automatically set maxJobs if not manually set
        if slurm_max_jobs is None:
            logger.warn(f'max_jobs not set - what to do in this case?')
            self.slurm_max_jobs = 1
        else:
            self.slurm_max_jobs = int(slurm_max_jobs)

        self.dry = dry


        self.action_list = []
        self.log_list    = []  # list of 2-tuples of the type: (log filename, type of action)
        self.futures     = []

        if self.backend == "slurm":
            # sensible defaults; override with slurm_opts
            LILFDIR = os.path.realpath(__file__).split('LiLF')[0] + 'LiLF'
            # We mount only the parent directory of the current working directory
            singularity_command = f"singularity exec --cleanenv --pwd {os.getcwd()} \
                --env PYTHONPATH=\$PYTHONPATH:{LILFDIR},PATH=\$PATH:{LILFDIR}/scripts/ --pid \
                --writable-tmpfs -B{os.path.dirname(os.getcwd())} {container_path}"
            self.slurm_max_walltime = slurm_max_walltime if slurm_max_walltime else get_slurm_max_walltime()[1],  # auto-find max walltime if not set

            so = {
                    'cores': min(self.max_cpus_per_node, 32),
                    'memory': f'{8*self.max_cpus_per_node}GB',
                    'walltime': self.slurm_max_walltime,
                    'python': 'python',
                    'log_directory': self.log_dir,
                    'job_script_prologue': [
                        "unset PYTHONPATH",
                        "unset PYTHONHOME",
                        "unset LD_LIBRARY_PATH",
                        f"{singularity_command} \\"
                    ]
                }
            
            if self.cluster == "Herts":
                so.update({'shebang': '#!/bin/tcsh'})
                so.update({'queue': 'core32'})
            else: 
                logger.warning(f'Slurm cluster {self.cluster} not specifically supported, trying generic settings.')
                
            self._cluster = SLURMCluster(**so)
            # Test if we want adaptive scaling if you like
            self._cluster.adapt(minimum=1, maximum=slurm_max_jobs)
            self._client = Client(self._cluster)
            
            logger.debug(f"Dask SLURM cluster script:\n{self._cluster.job_script()}")

        if backend == 'slurm':
            logger.info(f'SLURM scheduler initialised  for cluster {self.cluster}:{self.hostname} \
                         (slurm_max_jobs: {self.slurm_max_jobs}, max_cpus_per_node: {self.max_cpus_per_node}, \
                         slurm_max_walltime: {self.slurm_max_walltime}, slurm_mem_per_cpu: {slurm_mem_per_cpu})')
        else:
            logger.info(f'Local scheduler initialised  for cluster {self.cluster}:{self.hostname} \
                         (max_cpucores_per_node: {self.max_cpus_per_node})')

    def add(self, cmd='', log='', commandType='general', threads=1, mem=None, time='00:30:00', timeout=None):
        """
        Add a command with optional resources (mapped to SLURM/Dask).
        """
        
        if mem is None:
            mem =threads
        log_path = os.path.join(self.log_dir, log) if log else ''
        if log:
            # Truncate the log on first write
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            open(log_path, 'w').close()
            self.log_list.append((log_path, commandType))

        if commandType in ['wsclean', 'DP3', 'python']:
            logger.debug(f'Running {commandType}: {cmd}')
        else:
            logger.debug(f'Running general: {cmd}')
            

        self.action_list.append(dict(
            cmd=cmd, log=log_path, commandType=commandType,
            threads=threads, mem=mem, time=time, timeout=timeout
        ))

        if self.backend == "slurm":
            action = dict(
                cmd=cmd, log=log_path, commandType=commandType,
                threads=threads, mem=mem, time=time, timeout=timeout
            )
            #fut = self._client.submit(_run_cmd, cmd)
            fut = self._client.submit(_run_cmd_gpt, cmd, log_path, timeout, resources=None, pure=False)
            self.futures.append(fut, action)

    def run(self, check=False, maxProcs=None):
        if self.dry:
            return
        maxProcs_run = maxProcs

        if self.backend == "slurm":
            # Gather and raise on failure
            for fut, action in as_completed([f for f, _ in self.futures]):
                rc, wall = fut.result()
                if rc != 0:
                    tail = ''
                    if action['log'] and os.path.exists(action['log']):
                        tail = subprocess.check_output(f'tail -n 40 {shlex.quote(action["log"])}', shell=True).decode()
                    raise RuntimeError(f"Command failed (rc={rc}): {action['cmd']}\nLog: {action['log']}\n{tail}")
                
            self.futures.clear()

        else:
            # local thread pool (your existing behavior)
            from queue import Queue
            from threading import Thread
            q = Queue()

            def worker():
                for item in iter(q.get, None):
                    cmd, log, timeout, env = item['cmd'], item['log'], item['timeout'], item['env']
                    gc.collect()
                    rc, wall = _run_cmd(cmd, log, timeout, env)
                    if rc != 0:
                        tail = ''
                        if log and os.path.exists(log):
                            tail = subprocess.check_output(f'tail -n 40 {shlex.quote(log)}', shell=True).decode()
                        raise RuntimeError(f"Command failed (rc={rc}): {cmd}\nLog: {log}\n{tail}")
                    q.task_done()

            threads = [Thread(target=worker, daemon=True) for _ in range(maxProcs_run)]
            for t in threads: t.start()
            for a in self.action_list:
                q.put_nowait(a)
            q.join()
            for _ in threads: q.put(None)
            for t in threads: t.join()

        if check:
            for log, ctype in self.log_list:
                self.check_run(log, ctype)

        self.action_list.clear()
        self.log_list.clear()
        

        #if self.qsub:
        #    if qsub_cpucores == 'max':
        #        qsub_cpucores = self.max_cpucores
        #    # if number of cores not specified, try to find automatically
        #    elif qsub_cpucores == None:
        #        qsub_cpucores = 1 # default use single CPU
        #        if ("DP3" == cmd[ : 4]):
        #            qsub_cpucores = 1
        #        if ("wsclean" == cmd[ : 7]):
        #            qsub_cpucores = self.max_cpucores
        #    if (qsub_cpucores > self.max_cpucores):
        #        qsub_cpucores = self.max_cpucores
        #    self.action_list.append([str(qsub_cpucores), '\'' + cmd + '\''])
        #else:
        self.action_list.append(cmd)

        if (log != ""):
            self.log_list.append((log, commandType))


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
            #out += subprocess.check_output(r'grep -i -l "already a beam correction applied" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
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
            out += subprocess.check_output(r'grep -l "Bus error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "(core dumped)" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # out += subprocess.check_output('grep -L "Cleaning up temporary files..." '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "python"):
            out = subprocess.check_output(r'grep -l "Traceback (most recent call last):" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Segmentation fault\|Killed" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # out += subprocess.check_output(r'grep -i -l \'(?=^((?!error000).)*$).*Error.*\' '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -i -l "Critical" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "ERROR" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "ImportError" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Traceback (most recent call last)" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Permission denied" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)


        elif (commandType == "general"):
            out = subprocess.check_output('grep -l -i "error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

        else:
            logger.warning("Unknown command type for log checking: '" + commandType + "'")
            return 1

        if out != b'':
            out = out.split(b'\n')[0].decode()
            logger.error(commandType+' run problem on:\n'+out)
            errlines = subprocess.check_output('tail -n 10 '+log, shell = True, stderr = subprocess.STDOUT).decode()
            raise RuntimeError(commandType+' run problem on:\n'+out+'\n'+errlines)

        return 0


    def get_cluster(self):
        """
        TODO required?
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
        elif ('spider' in hostname):
            return "Spider"
        else:
            logger.debug('Hostname %s unknown.' % hostname)
            return "Unknown"
        
    def close(self):
        if self.backend == "slurm":
            self._client.close()
            self._cluster.close()
        else:
            pass

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
