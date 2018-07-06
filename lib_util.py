import os, sys, re, pickle, random, shutil, glob

from casacore import tables
import numpy as np

if (sys.version_info > (3, 0)):
    from configparser import ConfigParser
else:
    from ConfigParser import ConfigParser

# load here to be sure to have "Agg" at the beginning
import matplotlib as mpl
mpl.use("Agg")

from lib_log import logger

def getParset(parsetFile='../lilf.config'):
    """
    Get parset file and return dict of values
    """
    config = ConfigParser()
    config.read(parsetFile)
    
    # populate defaults sections
    if not config.has_section('flag'): config.add_section('flag')
    if not config.has_section('timesplit'): config.add_section('timesplit')
    if not config.has_section('model'): config.add_section('model')

    # flag
    if not config.has_option('flag', 'stations'): config.set('flag', 'stations', 'DE*;FR*;SE*;UK*')
    # timesplit
    if not config.has_option('timesplit', 'ngroups'): config.set('timesplit', 'ngroups', 2)
    # model
    if not config.has_option('model', 'sourcedb'): config.set('model', 'sourcedb', '')
    if not config.has_option('model', 'apparent'): config.set('model', 'apparent', False)
    if not config.has_option('model', 'useReg'): config.set('model', 'useReg', None)

    return config


def columnExists(tableObject, columnName):
    '''
    Check whether a column with name 'columnName' exists in  table 'tableObject'.
    '''
    columnNames = tableObject.colnames()
    return (columnName in columnNames)


def columnAddSimilar(pathMS, columnNameNew, columnNameSimilar, dataManagerInfoNameNew, overwrite = False, fillWithOnes = True, comment = "", verbose = False):
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


def distanceOnSphere(RAs1, Decs1, RAs2, Decs2):
    """
    Return the distances on the sphere from the set of points '(RAs1, Decs1)' to the
    set of points '(RAs2, Decs2)' using the spherical law of cosines.

    It assumes that all inputs are given in degrees, and gives the output in degrees, too.

    Using 'numpy.clip(..., -1, 1)' is necessary to counteract the effect of numerical errors, that can sometimes
    incorrectly cause '...' to be slightly larger than 1 or slightly smaller than -1. This leads to NaNs in the arccosine.
    """

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


def run_losoto(s, c, h5s, parsets):
    """
    s : scheduler
    c : cycle name, e.g. "final"
    h5s : lists of H5parm files
    parsets : lists of parsets to execute
    """

    logger.info("Running LoSoTo...")

    # concat
    if len(h5s) > 1:
        h5 = 'cal-'+c+'.h5'
        check_rm("cal-" + c + ".h5")
        s.add('H5parm_collector.py -V -s sol000 -o '+h5+' '+' '.join(h5s), log='losoto-'+c+'.log', commandType="python", processors='max')
        s.run(check = True)
    else:
        h5 = h5s[0]

    check_rm('plots')
    os.makedirs('plots')

    for parset in parsets:
        logger.debug('-- executing '+parset+'...')
        s.add('losoto -V '+h5+' '+parset, log='losoto-'+c+'.log', logAppend=True, commandType="python", processors='max')
        s.run(check = True)

    check_rm('plots-' + c)
    os.system('mv plots plots-' + c)


class Scheduler():
    def __init__(self, qsub = None, maxThreads = None, max_processors = None, log_dir = 'logs', dry = False):
        """
        qsub:           if true call a shell script which call qsub and then wait
                        for the process to finish before returning
        maxThreads:    max number of parallel processes
        dry:            don't schedule job
        max_processors: max number of processors in a node (ignored if qsub=False)
        """
        self.cluster = self.get_cluster()
        self.qsub    = qsub
        # if qsub/max_thread/max_processors not set, guess from the cluster
        # if they are set, double check number are reasonable
        if (self.qsub == None):
            if (self.cluster == "Hamburg"):
                self.qsub = True
            else:
                self.qsub = False
        else:
            if ((self.qsub == False and self.cluster == "Hamburg") or \
               (self.qsub == True and (self.cluster == "Leiden" or self.cluster == "CEP3"))):
                logger.critical('Qsub set to %s and cluster is %s.' % (str(qsub), self.cluster))
                sys.exit(1)

        if (maxThreads == None):
            if   (self.cluster == "Hamburg"):
                self.maxThreads = 32
            elif (self.cluster == "Leiden"):
                self.maxThreads = 64
            elif (self.cluster == "CEP3"):
                self.maxThreads = 40
            else:
                self.maxThreads = 12
        else:
            self.maxThreads = maxThreads

        if (max_processors == None):
            if   (self.cluster == "Hamburg"):
                self.max_processors = 6
            elif (self.cluster == "Leiden"):
                self.max_processors = 64
            elif (self.cluster == "CEP3"):
                self.max_processors = 40
            else:
                self.max_processors = 12
        else:
            self.max_processors = max_processors

        self.dry = dry
        logger.info("Scheduler initialised for cluster " + self.cluster + " (maxThreads: " + str(self.maxThreads) + ", qsub (multinode): " +
                     str(self.qsub) + ", max_processors: " + str(self.max_processors) + ").")

        self.action_list = []
        self.log_list    = [] # list of 2-tuples of the type: (log filename, type of action)

        if (not os.path.isdir(log_dir)):
            logger.info("Creating log dir '" + log_dir + "'.")
            os.makedirs(log_dir)
        self.log_dir = log_dir


    def get_cluster(self):
        """
        Find in which computing cluster the pipeline is running
        """
        import socket
        hostname = socket.gethostname()
        if (hostname == 'lgc1' or hostname == 'lgc2'):
            return "Hamburg"
        elif ('leidenuniv' in hostname):
            return "Leiden"
        elif (hostname[0 : 3] == 'lof'):
            return "CEP3"
        else:
            logger.error('Hostname %s unknown.' % hostname)
            return "Unknown"


    def add(self, cmd = '', log = '', logAppend = True, commandType = '', processors = None):
        """
        Add a command to the scheduler list
        cmd:         the command to run
        log:         log file name that can be checked at the end
        logAppend:  if True append, otherwise replace
        commandType: can be a list of known command types as "BBS", "DPPP", ...
        processors:  number of processors to use, can be "max" to automatically use max number of processors per node
        """
        if (log != ''):
            log = self.log_dir + '/' + log

            if (logAppend):
                cmd += " >> "
            else:
                cmd += " > "
            cmd += log + " 2>&1"

        if (processors != None and processors == 'max'):
            processors = self.max_processors

        if self.qsub:
            # if number of processors not specified, try to find automatically
            if (processors == None):
                processors = 1 # default use single CPU
                if ("calibrate-stand-alone" == cmd[ : 21]):
                    processors = 1
                if ("DPPP" == cmd[ : 5]):
                    processors = 1
                if ("wsclean" == cmd[ : 7]):
                    processors = self.max_processors
                if ("awimager" == cmd[ : 8]):
                    processors = self.max_processors
            if (processors > self.max_processors):
                processors = self.max_processors

            self.action_list.append([str(processors), '\'' + cmd + '\''])
        else:
            self.action_list.append(cmd)

        if (log != ""):
            self.log_list.append((log, commandType))

#    def add_casa(self, cmd = '', params = {}, wkd = None, log = '', logAppend = False, processors = None):
#        """
#        Run a casa command pickling the parameters passed in params
#        NOTE: running casa commands in parallel is a problem for the log file, better avoid
#        alternatively all used MS and CASA must be in a separate working dir
#
#        wkd = working dir (logs and pickle are in the pipeline dir)
#        """
#
#        if processors != None and processors == 'max': processors = self.max_processors
#        if processors == None: processors=self.max_processors # default use entire node
#
#        # since CASA can run in another dir, be sure log and pickle are in the pipeline working dir
#        if log != '': log = os.getcwd()+'/'+self.log_dir+'/'+log
#        pfile = os.getcwd()+'/casaparams_'+str(random.randint(0, 1e9))+'.pickle'
#        pickle.dump( params, open( pfile, "wb" ) )
#
#        # exec in the script dir?
#        if wkd == None: casacmd = 'casa --nogui --log2term --nologger -c '+cmd+' '+pfile
#        elif os.path.isdir(wkd):
#            casacmd = 'cd '+wkd+'; casa --nogui --log2term --nologger -c '+cmd+' '+pfile
#        else:
#            logger.error('Cannot find CASA working dir: '+wkd)
#            sys.exit(1)
#
#        if self.qsub:
#            if log != '' and not logAppend: casacmd = str(processors)+' \''+casacmd+' > '+log+' 2>&1'
#            elif log != '' and logAppend: casacmd = str(processors)+' \''+casacmd+' >> '+log+' 2>&1'
#            else: casacmd = str(processors)+' \''+casacmd
#
#            # clean up casa remnants in Hamburg cluster
#            if self.cluster == "Hamburg":
#                self.action_list.append(casacmd+'; killall -9 -r dbus-daemon Xvfb python casa\*\'')
#                if processors != self.max_processors:
#                    logger.error('To clean annoying CASA remnants no more than 1 CASA per node is allowed.')
#                    sys.exit(1)
#            else:
#                self.action_list.append(casacmd+'\'')
#        else:
#            if (log != '' and not logAppend):
#                self.action_list.append(casacmd + ' > ' + log + ' 2>&1')
#            elif (log != '' and logAppend):
#                self.action_list.append(casacmd + ' >> ' + log + ' 2>&1')
#            else:
#                self.action_list.append(casacmd)
#
#        if log != '':
#            self.log_list.append((log, "CASA"))


    def run(self, check = False, maxThreads = None):
        """
        If 'check' is True, a check is done on every log in 'self.log_list'.
        If max_thread != None, then it overrides the global values, useful for special commands that need a lower number of threads.
        """
        from threading import Thread
        from Queue import Queue
        import subprocess

        def worker(queue):
            for cmd in iter(queue.get, None):
                if self.qsub and self.cluster == "Hamburg":
                    # run in priority nodes
                    #cmd = 'salloc --job-name LBApipe --reservation=important_science --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                    #        ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                    # run on all cluster
                    cmd = 'salloc --job-name LBApipe --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                            ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
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
        import subprocess

        if (not os.path.exists(log)):
            logger.warning("No log file found to check results: " + log)
            return 1

        if (commandType == "BBS"):
            out = subprocess.check_output('grep -L success '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            if out != '':
                logger.error('BBS run problem on:\n'+out.split("\n")[0])
                return 1

        elif (commandType == "DPPP"):
            out = subprocess.check_output('grep -L "Finishing processing" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -l "Exception" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -l "**** uncaught exception ****" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            #out += subprocess.check_output('grep -l "misspelled" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            if out != '':
                logger.error('DPPP run problem on:\n'+out.split("\n")[0])
                return 1

        elif (commandType == "CASA"):
            out = subprocess.check_output('grep -l "[a-z]Error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -l "An error occurred running" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -l "\*\*\* Error \*\*\*" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            if out != '':
                logger.error('CASA run problem on:\n'+out.split("\n")[0])
                return 1

        elif (commandType == "wsclean"):
            out = subprocess.check_output('grep -l "exception occurred" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -L "Cleaning up temporary files..." '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            if out != '':
                logger.error('WSClean run problem on:\n'+out.split("\n")[0])
                return 1

        elif (commandType == "python"):
            out = subprocess.check_output('grep -l "Traceback (most recent call last):" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -i -l "Error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            out += subprocess.check_output('grep -i -l "Critical" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            if out != '':
                logger.error('Python run problem on:\n'+out.split("\n")[0])
                return 1

        elif (commandType == "general"):
            out = subprocess.check_output('grep -l -i "error" '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
            if out != '':
                logger.error('Run problem on:\n'+out.split("\n")[0])
                return 1

        else:
            logger.warning("Unknown command type for log checking: '" + commandType + "'")
            return 1

        return 0
