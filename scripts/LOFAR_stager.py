#!/usr/bin/env python3
# To do a second run, rename the uris.pickle

import os, sys, time, glob, pickle, argparse
from datetime import datetime
import subprocess, multiprocessing
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
import stager_access as stager

#project = 'LC9_017' # 3c first part
#project = 'LC10_020' # 3c second part
#project = 'LC13_011' # cluster

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--project', '-p', dest='project', help='')
parser.add_argument('--sasID', '-o', dest='sasID', help='')
args = parser.parse_args()

project = args.project
sasID = args.sasID

if project is None and sasID is None:
    print('ERROR: at least one between --project and --sasID needs to be specified.')
    sys.exit()

# The class of data to query
cls = CorrelatedDataProduct

downloaded_mss = glob.glob('*MS')
if os.path.exists('renamed.txt'):
    with open('renamed.txt','r') as flog:
        for line in flog:
            downloaded_mss.append(line[:-1]+'.MS')

# first collect all uris
if not os.path.exists('uris.pickle'):
    query_observations = Observation.select_all().project_only(project)
    uris = set() # All URIS to stage
    for observation in query_observations :
        print("Querying ObservationID %s" % observation.observationId, end='')
        # Instead of querying on the Observations of the DataProduct, all DataProducts could have been queried
        dataproduct_query = cls.observations.contains(observation)
        # isValid = 1 means there should be an associated URI
        dataproduct_query &= cls.isValid == 1
        #dataproduct_query &= cls.pipelineName == 1
        for i, dataproduct in enumerate(dataproduct_query):
            # This DataProduct should have an associated URL
            fileobject = ((FileObject.data_object == dataproduct) & (FileObject.isValid > 0)).max('creation_date')
            #print (fileobject.URI)
            if fileobject:
                if i%10 == 0:
                    print(".", end='')
                    sys.stdout.flush()
                skip = False
                for ms in downloaded_mss:
                    if ms in fileobject.URI:
                        print("%s: already downloaded in %s." % (fileobject.URI, ms) )
                        skip = True
                if not skip: uris.add(fileobject.URI)
            else :
                print("No URI found for %s with dataProductIdentifier %d" % (dataproduct.__class__.__name__, dataproduct.dataProductIdentifier))
        
        print("")
            
            #if len(uris) == 1: break # TEST
        #break # TEST
 
    pickle.dump(uris, open('uris.pickle', 'wb'))
else:
    uris = pickle.load(open('uris.pickle','rb'))

print(("Total URI's found %d" % len(uris)))
 
# Queue of data to stage
L_toStage = multiprocessing.Manager().list() # list of surls to download
L_inStage = multiprocessing.Manager().list() # list of sids of active staging processes

# Queue of data to download
L_toDownload = multiprocessing.Manager().list() # list of surls ready to download
L_inDownload = multiprocessing.Manager().list() # list of surls being downloaded

class Worker(multiprocessing.Process):

    def __init__(self, stager, L_toStage, L_inStage, L_toDownload, L_inDownload):
        multiprocessing.Process.__init__(self)
        self.exit = multiprocessing.Event()
        self.stager = stager
        self.L_toStage = L_toStage
        self.L_inStage = L_inStage
        self.L_toDownload = L_toDownload
        self.L_inDownload = L_inDownload

    def run(self):
        while not self.exit.is_set():
            raise "To be set."

    def terminate(self):
        self.exit.set()

class Worker_stager(Worker):
    def run(self):
        import time
        while not self.exit.is_set():
            # if there's space add a block of 500
            if len(self.L_inStage) < 10 and len(self.L_toStage) > 0:
                uris = self.L_toStage[:500]
                print("Stager -- Staging %i uris" % len(uris))
                sids = self.stager.stage(uris)
                self.L_inStage.append( sids )
                for uri in uris:
                    self.L_toStage.remove(uri)
    
            time.sleep(60)
    
    
class Worker_checker(Worker):
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
                    if not surl in self.L_toDownload and not surl in self.L_inDownload:
                        self.L_toDownload.append(surl)

                # pass to download
                if status == 'success' or status == 'partial success':
                    print("Checker -- Sid %i completed." % sid)
                    self.L_inStage.remove(sid)

                elif status == 'in progress' or status == 'new' or status == 'scheduled':
                    continue
    
                # reschedule
                else:
                    print("Checker -- ERROR: Sid %i status is %s!" % (sid, status) )
                    self.stager.reschedule(sid)
    
            time.sleep(60)
    
    
class Worker_downloader(Worker):
    def run(self):
        import os
        while not self.exit.is_set():
            if len(self.L_toDownload) > 0:
                surl = self.L_toDownload.pop()
                print("Downloader -- Download: "+surl.split('/')[-1])
                self.L_inDownload.append(surl)
                with open("wgetout.txt","wb") as out, open("wgeterr.txt","wb") as err:
                    if 'psnc.pl' in surl: p = subprocess.Popen('wget -nv https://lta-download.lofar.psnc.pl/lofigrid/SRMFifoGet.py?surl=%s -O - | tar -x' % surl, shell=True,stdout=out,stderr=err)
                    if 'sara.nl' in surl: p = subprocess.Popen('wget -nv https://lofar-download.grid.surfsara.nl/lofigrid/SRMFifoGet.py?surl=%s -O - | tar -x' % surl, shell=True,stdout=out,stderr=err)
                    os.waitpid(p.pid, 0)
                self.L_inDownload.remove(surl)

            time.sleep(2)

# start processes
w_stager = Worker_stager(stager, L_toStage, L_inStage, L_toDownload, L_inDownload)
w_checker = Worker_checker(stager, L_toStage, L_inStage, L_toDownload, L_inDownload)
w_downloader1 = Worker_downloader(stager, L_toStage, L_inStage, L_toDownload, L_inDownload)
w_downloader2 = Worker_downloader(stager, L_toStage, L_inStage, L_toDownload, L_inDownload)

# fill the queue with uris
[L_toStage.append(uri) for uri in uris]

# add things already staged
i=0
try:
    for sid, _ in stager.get_progress():
        sid = int(sid)
        L_inStage.append(sid) # the worker will take care of starting downloads
        for surl in stager.get_surls_online(sid):
            if surl in L_toStage:
                L_toStage.remove(surl)
                i+=1
    print("Removed %i already staged surls." % i)
except:
    pass

w_stager.start()
w_checker.start()
w_downloader1.start()
w_downloader2.start()

while True:
    sys.stdout.write("\r%s: To stage: %i -- In staging: %i -- To download: %i -- In downloading: %i || " % \
            ( time.ctime(), len(L_toStage), len(L_inStage)*500, len(L_toDownload), len(L_inDownload) ) )
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
