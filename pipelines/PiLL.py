#!/usr/bin/env python

import os, sys
from LiLF import lib_util
parset = lib_util.getParset(parsetFile='lilf.config')
LiLF_dir = os.path.dirname(lib_util.__file__)

# get parameters
# use lilf.config (this is also used by all other scripits)
working_dir = parset.get('PiLL','working_dir')

# query the database for data to process
# TODO

# setup directiories
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
os.chdir(working_dir)
os.system('cp lilf.config '+working_dir)
os.makedirs('download')
os.makedirs('cals-bkp')
os.makedirs('cals')
os.makedirs('tgts-bkp')
os.makedirs('tgts')

##########
# download
os.chdir(working_dir+'/download')
os.system(LiLF_dir+'/../scripts/LOFAR_stager.py')
# TODO: how to be sure all MS were downloaded?
os.system(LiLF_dir+'/pipelines/LOFAR_download.py')

##########
# calibrator
if calibrator_tables_available:
    # TODO: Download calibration tables

else:
    # TODO
    # how to download only calibrator and then only target?
    os.system('mv '+working_dir+'/download/*3[C|c][196|295|380]/*MS '+working_dir+'/cals-bkp/')
    os.chdir(working_dir+'/cals')
    os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')

    # TODO: Upload new calibration tables

##########
# timesplit
os.chdir(working_dir+'/tgts')
os.system(LiLF_dir+'/pipelines/LOFAR_timesplit.py')

##########
# selfcal
os.chdir(working_dir+'/tgts')
os.system(LiLF_dir+'/pipelines/LOFAR_self.py')

##########
# DD-cal
os.chdir(working_dir+'/tgts')
os.system(LiLF_dir+'/pipelines/LOFAR_dd-serial.py')

# update the database
# TODO
