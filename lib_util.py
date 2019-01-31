import os, sys, re, time, pickle, random, shutil, glob

from casacore import tables
import numpy as np
import multiprocessing

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
    def add_default(section, option, val):
        if not config.has_option(section, option): config.set(section, option, val)

    config = ConfigParser()
    config.read(parsetFile)
    
    # add pipeline sections and defaul parset dir:
    for pipeline in ['download','demix','cal','timesplit','self','dd', 'ateam']:
        if not config.has_section(pipeline): config.add_section(pipeline)
        if not config.has_option(pipeline, 'parset_dir'): config.set(pipeline, 'parset_dir', os.path.dirname(__file__)+'/parsets/LOFAR_'+pipeline)
    # add other sections
    if not config.has_section('flag'): config.add_section('flag')
    if not config.has_section('model'): config.add_section('model')

    # download
    add_default('download', 'fix_table', 'True') # fix bug in some old observations
    add_default('download', 'renameavg', 'True')
    add_default('download', 'flag_elev', 'True')
    # demix
    add_default('demix', 'data_dir', '../cals-bkp/')
    add_default('demix', 'demix_model', '/home/fdg/scripts/model/demix_all.skydb')
    # cal
    add_default('cal', 'imaging', 'False')
    add_default('cal', 'skymodel', os.path.dirname(__file__)+'/models/calib-simple.skydb')
    add_default('cal', 'data_dir', '../cals-bkp/')
    # timesplit
    add_default('timesplit', 'data_dir', '../tgts-bkp/')
    add_default('timesplit', 'cal_dir', '../cals/')
    add_default('timesplit', 'ngroups', '1')
    add_default('timesplit', 'initc', '0')
    # self
    # dd
    add_default('dd', 'maxniter', '10')

    # flag
    add_default('flag', 'stations', 'DE*;FR*;SE*;UK*;IR*;PL*')
    # model
    add_default('model', 'sourcedb', None)
    add_default('model', 'apparent', 'False')
    add_default('model', 'userreg', None)

    return config


def columnExists(tableObject, columnName):
    # more to lib_ms
    '''
    Check whether a column with name 'columnName' exists in  table 'tableObject'.
    '''
    columnNames = tableObject.colnames()
    return (columnName in columnNames)


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


def run_wsclean(s, logfile, MSs_files, **kwargs):
    """
    Use only for imaging - not for predict
    s : scheduler
    args : parameters for wsclean
    """
    
    wsc_parms = []

    # basic parms
    wsc_parms.append( '-reorder -j '+str(s.max_processors) )
    # other stanrdard parms
    wsc_parms.append( '-clean-border 1' )
    # temp dir
    if s.get_cluster() == 'Hamburg_fat' and not 'temp_dir' in kwargs.keys():
        wsc_parms.append( '-temp-dir /localwork.ssd' )
    # user defined parms
    for parm, value in kwargs.items():
        if parm == 'cont': 
            parm = 'continue'
            value = ''
        if parm == 'size': value = '%i %i' % (value, value)
        wsc_parms.append( '-%s %s' % (parm.replace('_','-'), str(value)) )

    # files
    wsc_parms.append( MSs_files )

    # create command string
    command_string = 'wsclean '+' '.join(wsc_parms)
    s.add(command_string, log=logfile, commandType='wsclean', processors='max')
    logger.debug('Running wsclean: %s' % command_string)
    s.run(check=True)


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
            if ((self.qsub == False and self.cluster == "Hamburg") or \
               (self.qsub == True and (self.cluster == "Leiden" or self.cluster == "CEP3" or self.cluster == "Hamburg_fat"))):
                logger.critical('Qsub set to %s and cluster is %s.' % (str(qsub), self.cluster))
                sys.exit(1)

        if (maxThreads == None):
            if   (self.cluster == "Hamburg"):
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
        logger.info("Scheduler initialised for cluster " + self.cluster + " (maxThreads: " + str(self.maxThreads) + ", qsub (multinode): " +
                     str(self.qsub) + ", max_processors: " + str(self.max_processors) + ").")

        self.action_list = []
        self.log_list    = [] # list of 2-tuples of the type: (log filename, type of action)


    def get_cluster(self):
        """
        Find in which computing cluster the pipeline is running
        """
        import socket
        hostname = socket.gethostname()
        if (hostname == 'lgc1' or hostname == 'lgc2'):
            return "Hamburg"
        elif ('node3' in hostname):
            return "Hamburg_fat"
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
        commandType: can be a list of known command types as "BBS", "DPPP", ...
        processors:  number of processors to use, can be "max" to automatically use max number of processors per node
        """

        #cmd = cmd.replace(';','\;') # escape ;

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
            out += subprocess.check_output('grep -i -l \'(?=^((?!error000).)*$).*Error.*\' '+log+' ; exit 0', shell = True, stderr = subprocess.STDOUT)
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
