#!/usr/bin/env python3

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

import os, argparse, sys, glob
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
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for observations. It must lead to a directory where there are subdirectories containing /ddcal and /mss-avg derived from calibration.')
parser.add_argument('-l', '--list', dest='cllist', action='store', default=None, type=str, help='Name of .fits file which lists Name, RA, DEC and z. Optionally an extraction region and a mask region can be added.')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, default=None, help='RA/DEC where to center the extraction in deg. Use if you wish to extract only one target.')
parser.add_argument('--z', dest='redshift', type=float, default=-99, help='Redshift of the target.')
parser.add_argument('--name', dest='name', type=str, default='target_extracted', help='Name of the target. Will be used to create the directory containing the extracted data.')

args = parser.parse_args()

cluster_list = args.cllist
coords = args.radec

if cluster_list is not None and coords is not None:
    logger.error('Provide either a fits file (-l) with Name, RA, DEC, z, or RA and DEC (--radec) in deg.')
    sys.exit()

if cluster_list is None and coords is None:
    logger.error('Provide either a fits file (-l) with Name, RA, DEC, z, or RA and DEC (--radec) in deg.')
    sys.exit()

pathdir = args.path
ztarget = args.redshift
targetname = args.name

if not pathdir:
    logger.error('Provide a path (-p) where to look for LBA observations.')
    sys.exit()

if cluster_list:
    cl_name = get_data(cluster_list,'Name')
    cl_ra = get_data(cluster_list,'RA')
    cl_dec = get_data(cluster_list,'DEC')
    cl_z = get_data(cluster_list, 'z')
else:
    cl_name = [targetname]
    cl_ra = [coords[0]]
    cl_dec = [coords[1]]
    cl_z = [ztarget]
cl_failed = [] # this is a list of targets where extraction failed to print at the end

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
    with w.if_todo(f'Extraction object {cluster}'):
        if not os.path.exists(cluster):
            os.system(f'mkdir {cluster}')
        os.system(f'cd {cluster}')
        with open(f'{cluster}/redshift_temp.txt', 'w') as f:
            writer = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
            if len(cl_name) == 1: # case if only one target to extract:
                writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n]])
            else: # case if multiple targets
                # TODO: see if we can simplify the following lines here
                if ext==1:
                    os.system(f'cp {cl_extractreg[n]} {cluster}')
                    cl_extractreg[n] = os.path.basename(cl_extractreg[n])
                    if mreg==0:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n]])
                    else:
                        os.system(f'cp {cl_maskreg[n]} {cluster}')
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n], cl_maskreg[n]])
                else:
                    if mreg==0:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n]])
                    else:
                        #writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n]], cl_maskreg[n])
                        logger.error('To specify a mask region, an extraction region must also be provided.')
                        break
        os.chdir(str(cluster))
        os.system(f'LOFAR_extract.py -p {pathdir}')
        # check LOFAR_extract log if pipeline finished as expected
        logfile = sorted(glob.glob('pipeline-extract_*.logger'))[-1]
        with open(logfile, 'r') as f:
            last_line = f.readlines()[-1]
            if not "Done" in last_line:
                logger.error(f'Something went wrong in the extraction of {cluster} - check the logfile {cluster}/{logfile}.')
                cl_failed.append(cluster)
                if (len(cl_name) > 1) and (n < len(cl_name) -1):
                    logger.warning(f'Continuing with the next extraction target.')
                os.chdir('../')
                raise lib_util.Skip
            else:
                logger.info(f'Cluster {cluster} has been extracted.')
                os.chdir('../')
        size = float(os.path.getsize(str(cluster))/1e+9) #get directory size in GB, delete if too small


if len(cl_failed) < len(cl_name): # if all or some are successfull
    logger.info('Concluded. Extracted datasets are in the mss-extract directory of each target.')
    logger.info('A new column DIFFUSE_SUB has been added to each successfully extracted .ms file.')
    logger.info('It contains the compact-source-subtracted visibilities.')
    logger.info('Nominal, high and source-subtracted low resolution images are in the /img directory of each target.')
else: # none worked :(
    logger.error(f'Extraction of all target(s): {",".join(cl_failed)} failed')
if (len(cl_failed) < len(cl_name)) and len(cl_failed): # case some but not all failed
    logger.warning(f'Extraction of target(s): {",".join(cl_failed)} failed')

os.system('rm -r logs_cluster_extraction.logger*') # directory not used in target_extract
#os.system('rm cluster_extraction.logger*')
