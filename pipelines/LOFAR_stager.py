#!/usr/bin/env python3

import os, sys, time, glob, pickle
from datetime import datetime
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
import subprocess, multiprocessing
import stager_access as stager

# Should the found files be staged ?
do_stage = False
# The project to query, LC2_035 has public data
project = 'LC9_017'
# The class of data to query
cls = CorrelatedDataProduct

downloaded_mss = glob.glob('*MS')
 
query_observations = Observation.select_all().project_only(project)
uris = set() # All URIS to stage
for observation in query_observations :
    print("Querying ObservationID %s" % observation.observationId)
    # Instead of querying on the Observations of the DataProduct, all DataProducts could have been queried
    dataproduct_query = cls.observations.contains(observation)
    # isValid = 1 means there should be an associated URI
    dataproduct_query &= cls.isValid == 1
    for i, dataproduct in enumerate(dataproduct_query):
        # This DataProduct should have an associated URL
        fileobject = ((FileObject.data_object == dataproduct) & (FileObject.isValid > 0)).max('creation_date')
        if fileobject :
            #print("URI found %s" % fileobject.URI)
            if i%10 == 0: print(".")
            skip = False
            for ms in downloaded_mss:
                if ms in fileobject.URI:
                    print("%s: already downloaded in %s." % (fileobject.URI, ms) )
                    skip = True
            if not skip: uris.add(fileobject.URI)
        else :
            print("No URI found for %s with dataProductIdentifier %d" % (dataproduct.__class__.__name__, dataproduct.dataProductIdentifier))
        
        #if len(uris) == 1: break # TEST
    #break # TEST
 
pickle.dump(uris, open('uris.pickle', 'wb'))
print("Total URI's found %d" % len(uris))
 
# Queue of data to stage
Q_toStage = multiprocessing.Manager().Queue()
L_inStage = multiprocessing.Manager().list()

# Queue of data to download
Q_toDownload = multiprocessing.Manager().Queue()
L_inDownload = multiprocessing.Manager().list()

class Worker(multiprocessing.Process):

    def __init__(self, stager, Q_toStage, L_inStage, Q_toDownload, L_inDownload):
        multiprocessing.Process.__init__(self)
        self.exit = multiprocessing.Event()
        self.stager = stager
        self.Q_toStage = Q_toStage
        self.L_inStage = L_inStage
        self.Q_toDownload = Q_toDownload
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
            # if there's space add a block of 100
            if ( len(self.L_inStage) + self.Q_toDownload.qsize() + len(self.L_inDownload) ) < 4800 and self.Q_toStage.qsize() > 0:
                uris = []
                for i in range(100):
                    try:
                        uris.append( self.Q_toStage.get(block=False) )
                    except:
                        break
    
                print ("Stager -- Staging %i uris" % len(uris))
                sids = self.stager.stage(uris)
                self.L_inStage.append( sids )
    
            time.sleep(5)
    
    
class Worker_checker(Worker):
    def run(self):
        import time
        while not self.exit.is_set():
            for sid in self.L_inStage:
    
                # pass to download
                if self.stager.get_status(sid) == 'success':
                    print ("Checker -- Sid %i ready." % sid)
                    surls = self.stager.get_surls_online(sid)
                    for surl in surls:
                        self.Q_toDownload.put(surl)
                    self.L_inStage.remove(sid)
    
                # reschedule
                if self.stager.get_status(sid) == 'fail':
                    pass
    
            time.sleep(2)
    
    
class Worker_downloader(Worker):
    def run(self):
        import os
        while not self.exit.is_set():
            if not Q_toDownload.empty():
                surl = self.Q_toDownload.get(block=False)
                print ("Downloader -- Download: "+surl)
                self.L_inDownload.append(surl)
                with open("wgetout.txt","wb") as out, open("wgeterr.txt","wb") as err:
                    p = subprocess.Popen('wget -nv https://lta-download.lofar.psnc.pl/lofigrid/SRMFifoGet.py?surl=%s -O - | tar -x' % surl, shell=True,stdout=out,stderr=err)
                    os.waitpid(p.pid, 0)
                self.L_inDownload.remove(surl)

            time.sleep(2)

# start processes
w_stager = Worker_stager(stager, Q_toStage, L_inStage, Q_toDownload, L_inDownload)
w_checker = Worker_checker(stager, Q_toStage, L_inStage, Q_toDownload, L_inDownload)
w_downloader1 = Worker_downloader(stager, Q_toStage, L_inStage, Q_toDownload, L_inDownload)
w_downloader2 = Worker_downloader(stager, Q_toStage, L_inStage, Q_toDownload, L_inDownload)

# fill the queue with uris
[Q_toStage.put(uri) for uri in uris]

w_stager.start()
w_checker.start()
w_downloader1.start()
w_downloader2.start()

while True:
    sys.stdout.write("\r%s: To stage: %i -- In staging: %i -- To download: %i -- In downloading: %i || " % ( time.ctime(), Q_toStage.qsize(), len(L_inStage), Q_toDownload.qsize(), len(L_inDownload) ) )
    sys.stdout.flush()
    time.sleep(2)
    
    # if all queues are empty, kill children and exit
    if Q_toStage.empty() and Q_toDownload.empty() and (len(L_inStage) + len(L_inDownload) == 0 ):
        print("Done.")
        break

w_stager.terminate()
w_checker.terminate()
w_downloader1.terminate()
w_downloader2.terminate()
