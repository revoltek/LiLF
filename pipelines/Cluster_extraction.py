# Pipeline for extraction of target region after LOFAR_dd-serial.
# This pipeline will subtract sources outside of the region and
# perform subsequent self-calibration.
# It can work with multiple pointings.
# An extraction region can be provided. If not,
# an automatic region will be produced based on the amount of
# flux around the source of interest.
# Amplitude calibration can be forced, excluded or set to auto.
# A userReg may be specified as clean mask.
# phSolMode can be used to solve either using phases or phaseandtec.

import os, argparse, sys
from LiLF import lib_util, lib_log
import astropy.io.fits as pyfits
import csv
import numpy as np

logger_obj = lib_log.Logger('cluster_extraction.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('cluster_extraction.walker')

print("""

  _      ____     _       ____          _            _____        _                      _    _               
 | |    | __ )   / \     |  _ \   __ _ | |_  __ _   | ____|__  __| |_  _ __  __ _   ___ | |_ (_)  ___   _ __  
 | |    |  _ \  / _ \    | | | | / _` || __|/ _` |  |  _|  \ \/ /| __|| '__|/ _` | / __|| __|| | / _ \ | '_ \ 
 | |___ | |_) |/ ___ \   | |_| || (_| || |_| (_| |  | |___  >  < | |_ | |  | (_| || (__ | |_ | || (_) || | | |
 |_____||____//_/   \_\  |____/  \__,_| \__|\__,_|  |_____|/_/\_\ \__||_|   \__,_| \___| \__||_| \___/ |_| |_|
                                                                                                            
            
""")

def get_data(fname,colname):
    data=pyfits.open(fname)
    data=data[1].data
    return data[colname]

parser = argparse.ArgumentParser(description='Extraction of HETDEX LBA Clusters')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for observations. It must lead to a place where there are directories containing /ddcal and /mss-avg subdirectories.')
parser.add_argument('-l', '--list', dest='cllist', action='store', default='', type=str, help='Name of .fits file which lists Name, RA, DEC and z. Optionally an extraction region and a mask region can be added.')

args = parser.parse_args()
cluster_list = args.cllist
pathdir = args.path

cl_name = get_data(cluster_list,'Name')
cl_ra = get_data(cluster_list,'RA')
cl_dec = get_data(cluster_list,'DEC')
cl_z = get_data(cluster_list, 'z')

try:
    cl_extractreg =  get_data(cluster_list,'EXTREG')
    ext=1
except:
    ext=0
    pass

try:
    cl_maskreg =  get_data(cluster_list,'MASKREG')
    mreg=1
except:
    mreg=0
    pass

for n, cluster in enumerate(cl_name):
    print('')
    with w.if_todo('Extraction object ' + str(cluster)):
        if not os.path.exists(str(cluster)):
            os.system('mkdir '+str(cluster))
        os.system('cd '+str(cluster))
        with open(str(cluster)+'/redshift_temp.txt', 'w') as f:
            writer = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
            if ext==1:
                if mreg==0:
                    writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n]])
                else:
                    writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n], cl_maskreg[n]])
            else:
                if mreg==0:
                    writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n]])
                else:
                    #writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n]], cl_maskreg[n])
                    logger.error('To specify a mask region, an extraction region must also be provided.')
                    break
        os.chdir(str(cluster))
        cmd = '/beegfs/bax8338/LiLF/pipelines/LOFAR_extract.py -p ' + str(pathdir)
        os.system(cmd)
        error_check = np.loadtxt('error_checker.txt')
        if error_check == 0:
            logger.error(f'Error: no observation covers object {str(cluster)} in path {str(pathdir)}.')
        else:
            logger.info('Cluster '+str(cluster)+' has been extracted.')
        os.chdir('../')
        size = float(os.path.getsize(str(cluster))/1e+9) #get directory size in GB, delete if too small
        # if size < 0.05:
        #     os.system(f'rm -r {str(cluster)}')

logger.info('Concluded. Extracted datasets are in the mss-extract directory of each target.')
logger.info('A new column DIFFUSE_SUB has been added to each extracted ms file.')
logger.info('It contains the compact-source-subtracted visibilities.')
logger.info('Nominal, high and source-subtracted low resolution images are in the /img directory of each target.')

