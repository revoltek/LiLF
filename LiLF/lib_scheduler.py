import os
import subprocess
import multiprocessing
from threading import Thread
from queue import Queue

from LiLF.lib_log import logger
from LiLF.lib_util import check_rm

class Scheduler():
    def __init__(self, max_proc = None, max_cpucores = None, log_dir = '.', dry_run = False):
        """
        max_proc:       max number of parallel processes
        max_cpucores:   max number of cpu cores usable in a node
        dry_run:        don't schedule job
        """
        self.log_dir = log_dir
        self.dry_run = dry_run

        if max_cpucores is None:
            # check if running in a slurm environment with a limited number of CPUs (less than cpu_count())
            slurm_cpus = os.getenv('SLURM_CPUS_ON_NODE', False)
            if slurm_cpus:
                self.max_cpucores = int(slurm_cpus)
            else:
                self.max_cpucores = multiprocessing.cpu_count()
        else:
            self.max_cpucores = max_cpucores

        if (max_proc is None) or (max_proc > self.max_cpucores):
            self.max_proc = self.max_cpucores
        else:
            self.max_proc = max_proc

        logger.info(f"Scheduler initialised: max_proc={self.max_proc}, max_cpucores={self.max_cpucores}.")

        self.action_list = []
        self.log_list    = []  # list of 2-tuples of the type: (log filename, type of action)


    def add(self, cmd = '', log = '', logAppend = True, commandType = ''):
        """
        Add a command to the scheduler list
        cmd:         the command to run
        log:         log file name that can be checked at the end
        logAppend:   if True append, otherwise replace
        commandType: can be a list of known command types as "wsclean", "DP3", ...
        """

        if log:
            log = os.path.join(self.log_dir, log)
            redirect = '>>' if logAppend else '>'
            cmd += f' {redirect} {log} 2>&1'

        if commandType == 'wsclean':
            logger.debug(f'Running wsclean: {cmd}')
        elif commandType == 'DP3':
            logger.debug(f'Running DP3: {cmd}')
        elif (commandType.lower() == "ddfacet" or commandType.lower() == 'ddf'):
            logger.debug(f'Running DDFacet: {cmd}')
        elif commandType == 'python':
            logger.debug(f'Running python: {cmd}')
        else:
            logger.debug(f'Running general: {cmd}')

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
        if log:
            self.log_list.append((log, commandType))


    def run(self, check = False, max_proc = None):
        """
        If 'check' is True, a check is done on every log in 'self.log_list'.
        If max_proc != None, then it overrides the global values, useful for special commands that need a lower number of threads.
        """

        def worker(queue):
            for cmd in iter(queue.get, None):
                #if self.qsub:
                #    cmd = 'salloc --job-name LBApipe --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                #            ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                subprocess.call(cmd, shell = True)

        # limit number of processes
        if max_proc is None:
            max_proc_run = self.max_proc
        else:
            max_proc_run = min(max_proc, self.max_proc)

        q       = Queue()
        threads = [Thread(target = worker, args=(q,)) for _ in range(max_proc_run)]

        for t in threads: # start workers
            t.daemon = True
            t.start()

        if not self.dry_run:
            for action in self.action_list:
                q.put_nowait(action)
        for _ in threads:
            q.put(None) # signal no more commands
        for t in threads:
            t.join()

        # check outcomes on logs
        if check:
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
            logger.warning(f'No log file found to check results: {log}')
            return 1

        if (commandType == "DP3"):
            out = subprocess.check_output(f'grep -L "Finishing processing" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(f'grep -El "Segmentation fault|Killed|Aborted \(core dumped\)|misspelled" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(f'grep -il "Exception" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            #out += subprocess.check_output(f'grep -i -l "already a beam correction applied" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # this interferes with the missingantennabehaviour=error option...
            # out += subprocess.check_output(f'grep -l "error" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "CASA"):
            out = subprocess.check_output(f'grep -El "[a-z]Error|An error occurred running|\*\*\* Error \*\*\*" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "wsclean"):
            out = subprocess.check_output(f'grep -El "exception occur|Segmentation fault|Killed|Aborted|Bus error|\(core dumped\)" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # out += subprocess.check_output(f'grep -L "Cleaning up temporary files..." {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType.lower() == "ddfacet" or commandType.lower() == 'ddf'):
            out = subprocess.check_output(f'grep -El "Traceback \(most recent call last\):|exception occur|raise Exception|Segmentation fault|Killed|killed by signal|Aborted" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "python"):
            out = subprocess.check_output(f'grep -El "Traceback \(most recent call last\):|Segmentation fault|Killed|ImportError|Permission denied|ERROR" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(f'grep -il "Critical" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)
            # out += subprocess.check_output(f'grep -i -l \'(?=^((?!error000).)*$).*Error.*\' {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)

        elif (commandType == "general"):
            out = subprocess.check_output(f'grep -l -i "error" {log} ; exit 0', shell = True, stderr = subprocess.STDOUT)

        else:
            logger.warning(f"Unknown command type for log checking: '{commandType}'")
            return 1

        if out != b'':
            out = out.split(b'\n')[0].decode()
            logger.error(f'{commandType} run problem on:\n{out}')
            errlines = subprocess.check_output(f'tail -n 10 {log}', shell = True, stderr = subprocess.STDOUT).decode()
            raise RuntimeError(f'{commandType} run problem on:\n{out}\n{errlines}')

        return 0
    

def run_losoto(s, h5s, parsets, logname='losoto.log', plots_dir=None, h5_out=None) -> None:
    """
    s         : scheduler
    h5s       : list of input H5parm files; if more than one they are concatenated into h5_out
    parsets   : list of LoSoTo parset files to execute sequentially
    logname   : log file name
    plots_dir : directory to collect plots into; if None, plots are kept in ./plots
    h5_out    : output H5parm name; if None, h5s[0] is used in place (only valid when len(h5s) == 1)
    """

    logger.info("Running LoSoTo...")
    assert isinstance(h5s, list)

    # Determine working h5 path
    if h5_out is None:
        if len(h5s) > 1:
            raise ValueError("h5_out must be specified when multiple input h5 files are provided")
        else:
            h5_out = h5s[0]
    else:
        h5_out_dir = os.path.dirname(h5_out)
        if h5_out_dir:
            os.makedirs(h5_out_dir, exist_ok=True)

    # concat/move
    if len(h5s) > 1:
        check_rm(h5_out)
        s.add(f'H5parm_collector.py -V -s sol000 -o {h5_out} {" ".join(h5s)}', log=f'{logname}', commandType="python")
        s.run(check=True)
    elif h5s[0] != h5_out:
        os.system(f'cp -r {h5s[0]} {h5_out}')

    check_rm('plots')
    os.makedirs('plots')

    for parset in parsets:
        logger.debug(f'-- executing {parset}...')
        s.add(f'losoto -V {h5_out} {parset}', log=f'{logname}', logAppend=True, commandType="python")
        s.run(check=True)

    if plots_dir is not None:
        os.makedirs(plots_dir, exist_ok=True)
        os.system(f'mv plots/* {plots_dir}')
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