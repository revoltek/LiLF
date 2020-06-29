#!/usr/bin/env python

import os, sys, glob
from surveys_db import SurveysDB
from LiLF import lib_util
parset = lib_util.getParset(parsetFile='lilf.config')
LiLF_dir = os.path.dirname(lib_util.__file__)

# get parameters
# use lilf.config (this is also used by all other scripits)
working_dir = parset.get('PiLL','working_dir')
redo_cal = parset.getboolean('PiLL','redo_cal')
download_file = parset.get('PiLL','download_file')
project = parset.get('PiLL','project')
target = parset.get('PiLL','target')

# fix non absolute paths
if working_dir[0] != '/':
    working_dir = os.get_cwd()+'/'+working_dir

def calibrator_tables_available(obsid):
    """
    check if calibrator data exist in the database
    """
    with SurveysDB(survey='lba',readonly=True) as sdb:
        r = sdb.get_field(obsid)
        return r['']

def local_calibrator_dirs(searchdir='.', obsid=None):
    """
    Return the dirname of the calibrators
    """
    if obsid is None:
        calibrators = glob.glob(searchdir+'/id*_3[C|c]196') + \
                  glob.glob(searchdir+'/id*_3[C|c]295') + \
                  glob.glob(searchdir+'/id*_3[C|c]380')
    else:
        calibrators = glob.glob(searchdir+'/id%i_3[C|c]196' % obsid) + \
                  glob.glob(searchdir+'/id%i_3[C|c]295' % obsid) + \
                  glob.glob(searchdir+'/id%i_3[C|c]380' % obsid)

    if len(calibrators) == 0: return None
    else: return calibrators


# query the database for data to process
if project == '' and target == '':
    pass
    target = 'xxx'

#######
# setup
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
if os.path.exists('lilf.config') and os.get_cwd() != working_dir: 
    os.system('cp lilf.config '+working_dir)

os.chdir(working_dir)
if not os.path.exists('download'): os.makedirs('download')

if download_file != '':
    os.system('cp %s download/html.txt' % download_file)

##########
# data download
os.chdir(working_dir+'/download')

if download_file == '':
    cmd = LiLF_dir+'/scripts/LOFAR_stager.py'
    if not redo_cal:
        cmd += ' --nocal'
    if project != '':
        cmd += ' --project %s' % project
    if target != '':
        cmd += ' --target %s' % target
    os.system(cmd)

# TODO: how to be sure all MS were downloaded?
os.system(LiLF_dir+'/pipelines/LOFAR_download.py')

os.chdir(working_dir)
os.system('mv download/mss/* ./')

calibrators = local_calibrator_dirs()
targets = [t for t in glob.glob('id*') if t not in calibrators]
print ('CALIBRATORS:', calibrators)
print ('TARGET:', targets)

for target in targets:

    ##########
    # calibrator
    obsid = int(target.split('_')[0][2:])
    if redo_cal or not calibrator_tables_available(obsid):
        # if calibrator not downaloaded, do it
        cal_dir = local_calibrator_dirs(working_dir, obsid)
        if cal_dir is None:
            os.chdir(working_dir+'/download')
            os.system(LiLF_dir+'/scripts/LOFAR_stager.py --cal --obsid %i' % obsid)
            os.system(LiLF_dir+'/pipelines/LOFAR_download.py')
    
            calibrator = local_calibrator_dirs('./mss/', obsid)[0]
            os.system('mv '+calibrator+' '+working_dir)

        os.chdir(local_calibrator_dirs(working_dir, obsid)[0])
        if not os.path.exists('data-bkp'):
            os.makedirs('data-bkp')
            os.system('mv *MS data-bkp')
        os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')
    
    ##########
    # timesplit
    os.chdir(working_dir+'/'+target)
    if not os.path.exists('data-bkp'):
        os.makedirs('data-bkp')
        os.system('mv *MS data-bkp')

    os.system(LiLF_dir+'/pipelines/LOFAR_timesplit.py')

# group targets with same name, assuming they are different pointings of the same dir
grouped_targets = set([t.split('_')[1] for t in targets])

for grouped_target in grouped_targets:
    if not os.path.exists(working_dir+'/'+grouped_target):
        os.makedirs(working_dir+'/'+grouped_target)
    os.chdir(working_dir+'/'+grouped_target)
    
    # collet mss
    if not os.path.exists('mss'):
        os.makedirs('mss')
        for i, tc in enumerate(glob.glob('../id*_'+target+'/mss/TC*MS')):
            tc_ren = 'TC%02i.MS' % i
            print('mv -r %s mss/%s' % (tc,tc_ren))
            os.system('mv -r %s mss/%s' % (tc,tc_ren))

    ##########
    # selfcal
    os.system(LiLF_dir+'/pipelines/LOFAR_self.py')

    ##########
    # DD-cal
    os.system(LiLF_dir+'/pipelines/LOFAR_dd-serial.py')

    # TODO: update the database



