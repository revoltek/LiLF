#!/usr/bin/env python

import os, sys, glob
from LiLF import lib_util
parset = lib_util.getParset(parsetFile='lilf.config')
LiLF_dir = os.path.dirname(lib_util.__file__)

# get parameters
# use lilf.config (this is also used by all other scripits)
working_dir = parset.get('PiLL','working_dir')

def calibrator_tables_available(target):
    """
    check if calibrator data exist in the database
    """
    # TODO
    return False

# query the database for data to process

#######
# setup
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
os.chdir(working_dir)
if os.path.exists('lilf.config'): 
    os.system('cp lilf.config '+working_dir)
if not os.path.exists('download'): os.makedirs('download')

##########
# download
os.chdir(working_dir+'/download')
#os.system(LiLF_dir+'/scripts/LOFAR_stager.py')
# TODO: how to be sure all MS were downloaded?
os.system(LiLF_dir+'/pipelines/LOFAR_download.py')

calibrators = glob.glob('mss/*3[C|c][196|295|380]')
targets = [t for t in glob.glob('mss/id*') if t not in calibrators]

os.system('mv '+working_dir+'/download/mss/* '+working_dir)

for target in targets:

    ##########
    # calibrator
    if not calibrator_tables_available(targets):
        os.chdir(working_dir+'/download')
        os.system(LiLF_dir+'/scripts/LOFAR_stager.py --cal --obsid ...')
        os.system(LiLF_dir+'/pipelines/LOFAR_download.py')

        calibrator = glob.glob('mss/id*')[0]
        os.system('mv mss/'+calibrator+' '+working_dir)

        os.chdir(working_dir+'/'+calibrator)
        os.makedirs('data-bkp')
        os.system('mv *MS data-kp')
        os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')
    
    ##########
    # timesplit
    os.chdir(working_dir+'/'+target)
    os.makedirs('data-bkp')
    os.system('mv *MS data-kp')
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



