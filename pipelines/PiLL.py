#!/usr/bin/env python

import os, sys

# get parameters
working_dir = '.'
LiLF_dir = '/home/baq1889/scripts/LiLF/'

# query the database for data to process


# setup directiories
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
os.chdir(working_dir)
os.makedirs('cals-bkp')
os.makedirs('cals')
os.makedirs('tgts-bkp')
os.makedirs('tgts')

# download
os.system(LiLF_dir+'/pipelines/LOFAR_download.py')

# calibrator
os.chdir(working_dir+'/cals')
os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')

# timesplit
os.chdir(working_dir+'/tgts')
os.system(LiLF_dir+'/pipelines/LOFAR_timesplit.py')

# selfcal
os.chdir(working_dir+'/tgts')
os.system(LiLF_dir+'/pipelines/LOFAR_self.py')

# DD-cal
os.chdir(working_dir+'/tgts')
os.system(LiLF_dir+'/pipelines/LOFAR_dd-serial.py')

# update the database
