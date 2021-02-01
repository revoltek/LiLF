#!/usr/bin/env python3
# To do a second run, rename the uris.pickle

# Idea of the code:
# The code crates 4 suprocesses:
# - 1 stager
# - 1 checker
# - 2 downlaoders

import os, sys, time, glob, pickle, argparse, re
import subprocess, multiprocessing
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
import stager_access as stager
from casacore import tables
from download_file import download_file

#project = 'LC9_017' # 3c first part
#project = 'LC10_020' # 3c second part
#project = 'LC13_011' # cluster

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--projects', '-p', dest='projects', help='Comma separated list of project names')
parser.add_argument('--obsID', '-o', dest='obsID', help='Comma separated list of project ids')
parser.add_argument('--target', '-t', dest='target', help='')
parser.add_argument('--calonly', '-c', dest='calonly', action='store_true', help='')
parser.add_argument('--nocal', '-n', dest='nocal', action='store_true', help='')
args = parser.parse_args()

if args.projects is None:
    print('ERROR: --project needs to be specified.')
    sys.exit()
projects = args.projects.split(',')

if args.obsID is not None:
    obsIDs = [int(obsID) for obsID in args.obsID.split(',')]

target = args.target
calonly = args.calonly
nocal = args.nocal

# Login/Passwd for LTA
login = None
password = None
file_rc = os.path.expanduser('~/.stagingrc')
if os.path.exists(file_rc):
    with open(file_rc, 'r') as f:
        for line in f:
            if 'user' in line:
                login = line.split('=')[-1].strip(' \t\n\r')
            elif 'password' in line:
                password = line.split('=')[-1].strip(' \t\n\r')
            else:
                print('Unknown content of .wgetrc: %s' % line)
                sys.exit()

# The class of data to query
cls = CorrelatedDataProduct
re_cal = re.compile('.*3[c|C](196|295|380).*')

# First: collect all uris
# This part of the code simply selects the uris to stage starting from the project names and the target name
if not os.path.exists('uris.pickle'):
    uris = set() # All URIS to stage
    for project in projects:
        print("Quering project: %s" % project)
        query_observations = Observation.select_all().project_only(project)
        for observation in query_observations :
            obsID = int(observation.observationId)
            if args.obsID is not None:
                if obsID not in obsIDs:
                    continue
            print("Querying ObservationID %i" % obsID, end='')
            # Instead of querying on the Observations of the DataProduct, all DataProducts could have been queried
            dataproduct_query = cls.observations.contains(observation)
            # isValid = 1 means there should be an associated URI
            dataproduct_query &= cls.isValid == 1
            #if target is not None: dataproduct_query &= CorrelatedDataProduct.subArrayPointing.targetName == target

            for i, dataproduct in enumerate(dataproduct_query):
                # apply selections
                name = dataproduct.subArrayPointing.targetName
                if re_cal.match(name) and nocal: continue
                if not re_cal.match(name) and calonly: continue
                if target is not None and not target in name: continue

                # This DataProduct should have an associated URL
                fileobject = ((FileObject.data_object == dataproduct) & (FileObject.isValid > 0)).max('creation_date')
                if fileobject:
                    uris.add(fileobject.URI)
                    if i%10 == 0:
                        print(".", end='')
                        sys.stdout.flush()
                else :
                    print("No URI found for %s with dataProductIdentifier %d" % (dataproduct.__class__.__name__, dataproduct.dataProductIdentifier))
            print("")
                
                #if len(uris) == 1: break # TEST
            #break # TEST
     
        pickle.dump(uris, open('uris.pickle', 'wb'))
else:
    uris = pickle.load(open('uris.pickle','rb'))

# remove files already downloaded/renamed
downloaded_mss = glob.glob('*MS')
if os.path.exists('renamed.txt'):
    with open('renamed.txt','r') as flog:
        for line in flog:
            downloaded_mss.append(line[:-1])
            
len_all_uris = len(uris)
uris = [uri for uri in uris if uri.split('/')[-1][:-13] not in downloaded_mss]
print(("Total URI's: %i (after removal of already downloaded: %i)" % (len_all_uris,len(uris))))
if len(uris) == 0:
    print("Done.")
    sys.exit()
 
# Queue of data to stage
L_toStage = multiprocessing.Manager().list()  # list of surls to download
L_inStage = multiprocessing.Manager().list()  # list of sids of active staging processes

# Queue of data to download
L_toDownload = multiprocessing.Manager().list()  # list of surls ready to download
L_inDownload = multiprocessing.Manager().list()  # list of surls being downloaded
L_Downloaded = multiprocessing.Manager().list()  # list of surls downlaoded

class Worker(multiprocessing.Process):
    """
    This is a global worker class to be inherited by the 3 specialized workers
    """
    def __init__(self, stager, L_toStage, L_inStage, L_toDownload, L_inDownload, L_Downloaded):
        multiprocessing.Process.__init__(self)
        self.exit = multiprocessing.Event()
        self.stager = stager
        self.L_toStage = L_toStage
        self.L_inStage = L_inStage
        self.L_toDownload = L_toDownload
        self.L_inDownload = L_inDownload
        self.L_Downloaded = L_Downloaded

    def run(self):
        while not self.exit.is_set():
            raise "To be set."

    def terminate(self):
        self.exit.set()

class Worker_stager(Worker):
    """
    This worker is specialized in staging the data
    """
    def run(self):
        import time
        while not self.exit.is_set():
            # if there's space add a block of 200
            if len(self.L_inStage) < 5 and len(self.L_toStage) > 0:
                uris = self.L_toStage[:200]
                #uris = [self.L_toStage[0]] # debug to stage 1 uri at a time
                print("Stager -- Staging %i uris" % len(uris))
                try:
                    sids = self.stager.stage(uris)
                    self.L_inStage.append( sids )
                    for uri in uris:
                        self.L_toStage.remove(uri)
                except Exception as e:
                    print("Error at staging...", e)
    
            time.sleep(600)
    
    
class Worker_checker(Worker):
    """
    This worker check if the data are staged and in that case it passes the uris to the downloaders
    """
    def run(self):
        import time
        while not self.exit.is_set():
            for sid in self.L_inStage:

                try:
                    status = self.stager.get_status(sid)
                    surls = self.stager.get_surls_online(sid)
                except:
                    status = ""
                    surls = []
                    print("Checker -- Failed to get status for sid %i. Continue." % sid)

                for surl in surls:
                    if not surl in self.L_toDownload and not surl in self.L_inDownload and not surl in self.L_Downloaded:
                        # this should always be the case, but if the process is re-started it might have collected old staging processes,
                        # this if prevents us from downloading useless files
                        if surl in uris:
                            self.L_toDownload.append(surl)

                # pass to download
                if status == 'success' or status == 'partial success':
                    print("Checker -- Sid %i completed." % sid)
                    self.L_inStage.remove(sid)

                elif status == 'in progress' or status == 'new' or status == 'scheduled' or status == '':
                    print("Checker -- WARNING: Sid %i status is: '%s'" % (sid, status) )
                    continue
    
                # reschedule
                else:
                    print("Checker -- ERROR: Sid %i status is: '%s' (reschedule submitted)" % (sid, status) )
                    self.stager.reschedule(sid)
    
            time.sleep(300)
    
    
class Worker_downloader(Worker):
    """
    This worker download the data
    """
    def run(self):
        import os
        while not self.exit.is_set():
            if len(self.L_toDownload) > 0:
                surl = self.L_toDownload.pop()
                self.L_inDownload.append(surl)

                tar_file = surl.split('/')[-1]  # e.g. .../L769079_SB020_uv.MS_daf24388.tar
                ms_file = surl.split('/')[-1].split('.MS')[0]+'.MS'  # e.g. .../L769079_SB020_uv.MS

                if 'psnc.pl' in surl:
                    url = 'https://lta-download.lofar.psnc.pl/lofigrid/SRMFifoGet.py?surl=%s' % surl
                    LTA_site = 'PL'
                elif 'sara.nl' in surl:
                    url = 'https://lofar-download.grid.surfsara.nl/lofigrid/SRMFifoGet.py?surl=%s' % surl
                    LTA_site = 'NL'
                elif 'juelich.de' in surl:
                    url = 'https://lofar-download.fz-juelich.de/webserver-lofar/SRMFifoGet.py?surl=%s' % surl
                    LTA_site = 'DE'
                else:
                    print('ERROR: unknown archive for %s...' % surl)
                    sys.exit()

                print("Downloader -- Download: %s (from: %s) " % (tar_file, LTA_site))

                # loop until the sanity check on the downloaded MS is ok
                while True:
                    download_file(url, tar_file, login, password)
                    os.system('tar xf %s' % tar_file)
                    print(tar_file)
                    try:
                        t = tables.table(ms_file, ack=False)
                        break
                    except:
                        print('ERROR opening %s, probably corrupted - redownload it' % ms_file)
                        #os.system('rm -r %s %s' % (tar_file, ms_file))
                os.system('rm -r %s' % tar_file)

                self.L_inDownload.remove(surl)
                self.L_Downloaded.append(surl)

            time.sleep(2)

# start processes
w_stager = Worker_stager(stager, L_toStage, L_inStage, L_toDownload, L_inDownload, L_Downloaded)
w_checker = Worker_checker(stager, L_toStage, L_inStage, L_toDownload, L_inDownload, L_Downloaded)
w_downloader1 = Worker_downloader(stager, L_toStage, L_inStage, L_toDownload, L_inDownload, L_Downloaded)
w_downloader2 = Worker_downloader(stager, L_toStage, L_inStage, L_toDownload, L_inDownload, L_Downloaded)

# fill the queue with uris
[L_toStage.append(uri) for uri in uris]

# add things already staged (here we do the get_progress())
i=0
try:
    for sid, _ in stager.get_progress().items():
        sid = int(sid)
        print("Found an active staging process: %i" % sid)
        for surl in stager.get_surls_online(sid):
            if surl in L_toStage:
                L_toStage.remove(surl)
                if sid not in L_inStage:
                    L_inStage.append(sid)  # the worker will take care of starting downloads
                i += 1
    print("Removed %i already staged surls." % i)
except Exception as e:
    print("Error recovering staged surls...", e)
    pass

w_stager.start()
w_checker.start()
w_downloader1.start()
w_downloader2.start()

# this part creates some output to monitor the progress
while True:
    sys.stdout.write("\r%s: To stage: %i -- In staging: %i (blocks) -- To download: %i -- In downloading: %i || " % \
            ( time.ctime(), len(L_toStage), len(L_inStage), len(L_toDownload), len(L_inDownload) ) )
    #print(L_toStage,L_inStage,L_toDownload,L_inDownload)
    sys.stdout.flush()
    time.sleep(2)
    
    # if all queues are empty, kill children and exit
    if len(L_toStage) + len(L_toDownload) + len(L_inStage) + len(L_inDownload) == 0 :
        print("Done.")
        break

w_stager.terminate()
w_checker.terminate()
w_downloader1.terminate()
w_downloader2.terminate()