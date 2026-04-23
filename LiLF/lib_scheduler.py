import os
import subprocess
import multiprocessing
from threading import Thread
from queue import Queue

from LiLF.lib_log import logger

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