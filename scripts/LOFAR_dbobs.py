#!/usr/bin/env python3

import os, sys, argparse, re
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation, Pointing, SubArrayPointing
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
from surveys_db import SurveysDB

survey_projects = 'LT14_002,LC12_017,LC9_016,LC8_031' # list of projects related with the LBA survey

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--projects', '-p', dest='project', help='Comma separated list of project names', default=survey_projects)
args = parser.parse_args()

projects = args.project.split(',')
re_cal = re.compile('3[c|C](196|295|380)')

if projects is None:
    print('ERROR: --project needs to be specified.')
    sys.exit()

# get obsid already done
with SurveysDB(survey='lba',readonly=True) as sdb:
    sdb.execute('select obs_id from field_obs')
    obs_to_skip = [ x['obs_id'] for x in sdb.cur.fetchall() ]

uris = set() # All URIS to stage
for project in projects:
    print('project %s' % project)
    query_observations = Observation.select_all().project_only(project)
    for observation in query_observations:
        obsID = int(observation.observationId)
        if obsID in obs_to_skip: continue
        print('check %i' % obsID)
        dataproduct_query = CorrelatedDataProduct.observations.contains(observation)
        # isValid = 1 means there should be an associated URI
        dataproduct_query &= CorrelatedDataProduct.isValid == 1
        dataproduct_query &= CorrelatedDataProduct.minimumFrequency >= 59
        dataproduct_query &= CorrelatedDataProduct.maximumFrequency <= 59.3
        for i, dataproduct in enumerate(dataproduct_query):
            # apply selections
            name = dataproduct.subArrayPointing.targetName.split('_')[-1]
            if re_cal.match(name): continue
            print('Add to the db: %i -> %s' % (obsID, name))
            with SurveysDB(survey='lba',readonly=False) as sdb:
                sdb.execute('INSERT INTO field_obs (obs_id,field_id) VALUES (%i,"%s")' % (obsID, name))

 
len_all_uris = len(uris)
uris = [uri for uri in uris if uri.split('/')[-1][:-13] not in downloaded_mss]
print(("Total URI's: %i (after removal of already downloaded: %i)" % (len_all_uris,len(uris))))
 
