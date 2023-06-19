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
parser.add_argument('-lp', '--lilfpath', dest='lilfpath', action='store', default='', type=str, help='Path where to look for LiLF.')
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

pathlilf = args.lilfpath
pathdir = args.path
ztarget = args.redshift
targetname = args.name

if not pathlilf:
    logger.error('Provide a path (-lp) where to look for LiLF.')
    sys.exit()

if not pathdir:
    logger.error('Provide a path (-p) where to look for LBA observations.')
    sys.exit()

if cluster_list:
    multiple = True
    cl_name = get_data(cluster_list,'Name')
    cl_ra = get_data(cluster_list,'RA')
    cl_dec = get_data(cluster_list,'DEC')
    cl_z = get_data(cluster_list, 'z')
else:
    multiple = False
    cl_name = targetname
    cl_ra = coords[0]
    cl_dec = coords[1]
    cl_z = ztarget

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

if multiple == True:
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
            cmd = f'{pathlilf}/LiLF/pipelines/LOFAR_extract.py -p ' + str(pathdir)
            os.system(cmd)
            error_check = np.loadtxt('error_checker.txt')
            if error_check == 0:
                logger.error(f'Error: no observation covers object {str(cluster)} in path {str(pathdir)}.')
                logger.error(f'If this is somehow unexpected, check the path (-p) and the coordinates and try again.')
                logger.error(f'If you wish to force the extraction, you can lower the beam sensitivity threshold (default = 0.3) in lilf.config.')
                os.system(f'rm cluster_extraction.walker')
                sys.exit()
            else:
                logger.info('Cluster '+str(cluster)+' has been extracted.')
            os.chdir('../')
            size = float(os.path.getsize(str(cluster))/1e+9) #get directory size in GB, delete if too small

else:
    print('')
    cluster = cl_name
    with w.if_todo('Extraction object ' + str(cluster)):
        if not os.path.exists(str(cluster)):
            os.system('mkdir '+str(cluster))
        os.system('cd '+str(cluster))
        with open(str(cluster)+'/redshift_temp.txt', 'w') as f:
            writer = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
            writer.writerow([cl_z, cl_ra, cl_dec, cl_name])
        os.chdir(str(cluster))
        cmd = f'{pathlilf}/LiLF/pipelines/LOFAR_extract.py -p ' + str(pathdir)
        os.system(cmd)
        error_check = np.loadtxt('error_checker.txt')
        if error_check == 0:
            logger.error(f'Error: no observation covers object {str(cluster)} in path {str(pathdir)}, or the beam sensitivity at its position is too low.')
            logger.error(f'If this is somehow unexpected, check the path (-p) and the coordinates and try again.')
            logger.error(f'If you wish to force the extraction, you can lower the beam sensitivity threshold (default = 0.3) in lilf.config.')
            os.system(f'rm cluster_extraction.walker')
            sys.exit()
        else:
            logger.info('Target '+str(cluster)+' has been extracted.')
        os.chdir('../')
        size = float(os.path.getsize(str(cluster))/1e+9) #get directory size in GB, delete if too small


logger.info('Concluded. Extracted datasets are in the mss-extract directory of each target.')
logger.info('A new column DIFFUSE_SUB has been added to each extracted .ms file.')
logger.info('It contains the compact-source-subtracted visibilities.')
logger.info('Nominal, high and source-subtracted low resolution images are in the /img directory of each target.')

