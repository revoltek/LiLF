#!/usr/bin/env python3
# This script query the LTA and populate the field_obs table
# At the same time set to "Observed" all fields that have at least 3 observed hours

# To reset everything:
#with SurveysDB(survey='lba',readonly=False) as sdb: 
#   sdb.execute('DELETE FROM field_obs') 
#   sdb.execute('UPDATE fields SET status="Not started"') 

# To reset all runs:
#with SurveysDB(survey='lba',readonly=False) as sdb: 
#   sdb.execute('UPDATE fields SET status="Observed" where status!="Not started"') 

import os, sys, argparse, re
from awlofar.database.Context import context
from awlofar.main.aweimports import CorrelatedDataProduct, \
    FileObject, \
    Observation, Pointing, SubArrayPointing
from awlofar.toolbox.LtaStager import LtaStager, LtaStagerError
from surveys_db import SurveysDB

survey_projects = 'LT14_002,LC12_017,LC9_016,LC8_031' # list of projects related with the LBA survey

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--projects', '-p', dest='project', help='Comma separated list of project names.', default=survey_projects)
parser.add_argument('--skip', '-s', action="store_true", help='Skip observations already present in field_obs, \
        this is faster but might miss some target to update as "Observed" in the field table.')
args = parser.parse_args()

projects = args.project.split(',')
skip_obs = args.skip
re_cal = re.compile('3[c|C](196|295|380)')

if projects is None:
    print('ERROR: --project needs to be specified.')
    sys.exit()

# get obs_id already done
if skip_obs:
    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('select obs_id from field_obs')
        obs_to_skip = [ x['obs_id'] for x in sdb.cur.fetchall() ]
    print('Skip the following obs:', obs_to_skip)

id_all={}
with SurveysDB(survey='lba',readonly=False) as sdb:
    for project in projects:
        print('Checking project: %s' % project)
        query_observations = Observation.select_all().project_only(project)
        for observation in query_observations:
            obs_id = int(observation.observationId)
            id_all[obs_id]=[]
            if skip_obs and obs_id in obs_to_skip: continue
            print('Checking obs_id: %i' % obs_id)
            dataproduct_query = CorrelatedDataProduct.observations.contains(observation)
            # isValid = 1 means there should be an associated URI
            dataproduct_query &= CorrelatedDataProduct.isValid == 1
            dataproduct_query &= CorrelatedDataProduct.minimumFrequency >= 59
            dataproduct_query &= CorrelatedDataProduct.maximumFrequency <= 59.3
            for i, dataproduct in enumerate(dataproduct_query):
                # apply selections
                field_id = dataproduct.subArrayPointing.targetName.split('_')[-1]
                if re_cal.match(field_id): continue
                if field_id in id_all[obs_id]: continue # prevent multiple entries
                print('Add to the db: %i -> %s' % (obs_id, field_id))
                sdb.execute('INSERT INTO field_obs (obs_id,field_id) VALUES (%i,"%s")' % (obs_id, field_id))
                id_all[obs_id].append(field_id)

#print (id_all)
field_id_all = []
for obs_id, field_id in id_all.items():
    field_id_all += field_id

#print (field_id_all)
with SurveysDB(survey='lba',readonly=False) as sdb:
    for field_id in set(field_id_all):
        nobs = field_id_all.count(field_id)
        if nobs >= 3:
            print("Set %s as observed (%i)" % (field_id, nobs))
            sdb.execute('UPDATE fields SET status="Observed" WHERE id="%s"' % (field_id))
