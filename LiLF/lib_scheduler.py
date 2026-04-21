
class Scheduler():
    def __init__(self, maxProcs = None, max_cpucores = None, log_dir = 'logs', dry = False):
        """
        maxProcs:       max number of parallel processes
        max_cpucores:   max number of cpu cores usable in a node
        dry:            don't schedule job
        """
        self.hostname = socket.gethostname()
        self.cluster = self.get_cluster()
        self.log_dir = log_dir
        #self.qsub    = qsub
        # if qsub/max_thread/max_cpucores not set, guess from the cluster
        # if they are set, double check number are reasonable
        #if (self.qsub == None):
        #    self.qsub = False
        #else:
        #    if ((self.qsub is False and self.cluster == "Hamburg") or
        #       (self.qsub is True and (self.cluster == "Leiden" or self.cluster == "CEP3" or
        #                               self.cluster == "Hamburg_fat" or self.cluster == "Pleiadi" or self.cluster == "Herts"))):
        #        logger.critical('Qsub set to %s and cluster is %s.' % (str(qsub), self.cluster))
        #        sys.exit(1)

        if (max_cpucores == None):
            # check if running in a slurm environment with a limited number of CPUs (less than cpu_count())
            slurm_cpus = os.getenv('SLURM_CPUS_ON_NODE', False)
            if slurm_cpus:
                self.max_cpucores = int(slurm_cpus)
            else:
                self.max_cpucores = multiprocessing.cpu_count()
        else:
            self.max_cpucores = max_cpucores

        if (maxProcs is None) or (maxProcs > self.max_cpucores):
            self.maxProcs = self.max_cpucores
        else:
            self.maxProcs = maxProcs

        self.dry = dry

        logger.info("Scheduler initialised for cluster " + self.cluster + ": " + self.hostname +
                    " (maxProcs: " + str(self.maxProcs) + ", max_cpucores: " + str(self.max_cpucores) + ").")

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
        elif ('spider' in hostname):
            return "Spider"
        else:
            logger.debug('Hostname %s unknown.' % hostname)
            return "Unknown"


    def add(self, cmd = '', log = '', logAppend = True, commandType = ''):
        """
        Add a command to the scheduler list
        cmd:         the command to run
        log:         log file name that can be checked at the end
        logAppend:   if True append, otherwise replace
        commandType: can be a list of known command types as "wsclean", "DP3", ...
        """

        if (log != ''):
            log = self.log_dir + '/' + log

            if (logAppend):
                cmd += " >> "
            else:
                cmd += " > "
            cmd += log + " 2>&1"

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


    def run(self, check = False, maxProcs = None):
        """
        If 'check' is True, a check is done on every log in 'self.log_list'.
        If max_thread != None, then it overrides the global values, useful for special commands that need a lower number of threads.
        """

        def worker(queue):
            for cmd in iter(queue.get, None):
                #if self.qsub and self.cluster == "Hamburg":
                #    cmd = 'salloc --job-name LBApipe --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                #            ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                gc.collect()
                subprocess.call(cmd, shell = True)

        # limit number of processes
        if (maxProcs == None):
            maxProcs_run = self.maxProcs
        else:
            maxProcs_run = min(maxProcs, self.maxProcs)

        q       = Queue()
        threads = [Thread(target = worker, args=(q,)) for _ in range(maxProcs_run)]

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
            out += subprocess.check_output(r'grep -l "ImportError" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Traceback (most recent call last)" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output(r'grep -l "Permission denied" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)

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
            errlines = subprocess.check_output('tail -n 10 '+log, shell = True, stderr = subprocess.STDOUT).decode()
            raise RuntimeError(commandType+' run problem on:\n'+out+'\n'+errlines)

        return 0