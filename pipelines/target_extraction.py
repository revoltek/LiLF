#!/usr/bin/env python3

# Pipeline for extraction of target region after LOFAR_dd.
# This pipeline will subtract sources outside of the region and
# perform subsequent self-calibration.
# It can work with multiple pointings.
# Check README file for additional detail and parameters.

import os, argparse, sys, glob
from LiLF import lib_util, lib_log
import astropy.io.fits as pyfits
import csv

logger_obj = lib_log.Logger('target_extraction.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('target_extraction.walker')

print("""

  _      ____     _       ____          _            _____        _                      _    _               
 | |    | __ )   / \     |  _ \   __ _ | |_  __ _   | ____|__  __| |_  _ __  __ _   ___ | |_ (_)  ___   _ __  
 | |    |  _ \  / _ \    | | | | / _` || __|/ _` |  |  _|  \ \/ /| __|| '__|/ _` | / __|| __|| | / _ \ | '_ \ 
 | |___ | |_) |/ ___ \   | |_| || (_| || |_| (_| |  | |___  >  < | |_ | |  | (_| || (__ | |_ | || (_) || | | |
 |_____||____//_/   \_\  |____/  \__,_| \__|\__,_|  |_____|/_/\_\ \__||_|   \__,_| \___| \__||_| \___/ |_| |_|
                                                                                                            
            
""")

def get_data(fname,colname):
    data=pyfits.open(fname).

    data=data[1].data
    return data[colname]

parser = argparse.ArgumentParser(description='Extraction of HETDEX LBA Clusters')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for observations. It must lead to a directory where subdirectories contain /ddcal and /mss-avg derived from calibration.')
parser.add_argument('-l', '--list', dest='cllist', action='store', default=None, type=str, help='Name of .fits file which lists Name, RA, DEC and z. Optionally an extraction region and a mask region can be added.')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, default=None, help='RA/DEC where to center the extraction in deg. Use if you wish to extract only one target.')
parser.add_argument('--z', dest='redshift', type=float, default=-99, help='Redshift of the target. Not necessary unless one wants to perform compact source subtraction.')
parser.add_argument('--name', dest='name', type=str, default='target_extracted', help='Name of the target. Will be used to create the directory containing the extracted data.')
parser.add_argument('--beamcut', dest='beamcut', type=float, default=0.3, help='Beam sensitivity threshold.')
parser.add_argument('--noselfcal', dest='noselfcal', help='Do not perform selfcalibration.', action='store_true')
parser.add_argument('--extreg', dest='extreg', action='store', default=None, type=str, help='Provide an optional extraction region. If not, one will be created automatically.')
parser.add_argument('--maskreg', dest='maskreg', action='store', default=None, type=str, help='Provide an optional user mask for cleaning.')
parser.add_argument('--ampcal', dest='ampcal', action='store', default='auto', type=str, help='Perform amplitude calibration. Can be set to True, False or auto.')
parser.add_argument('--ampsol', dest='ampsol', action='store', default='diagonal', type=str, help='How to solve for amplitudes. Can be set to diagonal or fulljones.')
parser.add_argument('--phsol', dest='phsol', action='store', default='tecandphase', type=str, help='How to solve for phases. Can be set to tecandphase or phase.')
parser.add_argument('--maxniter', dest='maxniter', type=float, default=10, help='Maximum number of selfcalibration cycles to perform.')
parser.add_argument('--subreg', dest='subreg', action='store', default=None, type=str, help='Provide an optional mask for sources that need to be removed.')

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
beam_cut = args.beamcut
no_selfcal = args.noselfcal
cl_extractreg = args.extreg
userReg = args.maskreg
ampcal = args.ampcal
ampsol = args.ampsol
phsol = args.phsol
maxniter = args.maxniter
subtract_reg_file = args.subreg

if not pathdir:
    logger.error('Provide a path (-p) where to look for LBA observations.')
    sys.exit()

if cluster_list:
    cl_name = get_data(cluster_list,'Name')
    cl_ra = get_data(cluster_list,'RA')
    cl_dec = get_data(cluster_list,'DEC')
    cl_z = get_data(cluster_list, 'z')

    try:
        cl_extractreg = get_data(cluster_list, 'EXTREG')
        ext = 1
    except:
        ext = 0
        pass

    try:
        cl_maskreg = get_data(cluster_list, 'MASKREG')
        mreg = 1
    except:
        mreg = 0
        pass

    try:
        cl_subreg = get_data(cluster_list, 'SUBREG')
        subreg = 1
    except:
        subreg = 0
        pass

else:
    cl_name = [targetname]
    cl_ra = [coords[0]]
    cl_dec = [coords[1]]
    cl_z = [ztarget]

    if cl_extractreg is not None:
        ext = 1
        if userReg is not None:
            mreg = 1
            cl_maskreg = userReg
        else:
            mreg = 0
        if subtract_reg_file is not None:
            subreg = 1
        else:
            subreg = 0
    else:
        if userReg is not None:
            logger.error('To specify a mask region, an extraction region must also be provided.')
            sys.exit()
        elif subtract_reg_file is not None:
            logger.error('To specify a subtraction region, an extraction region must also be provided.')
            sys.exit()
        else:
            ext = 0

cl_failed = [] # this is a list of targets where extraction failed to print at the end

for n, cluster in enumerate(cl_name):
    print('')
    with w.if_todo(f'Extraction target {cluster}'):
        if not os.path.exists(cluster):
            os.system(f'mkdir {cluster}')
        os.system(f'cd {cluster}')
        with open(f'{cluster}/redshift_temp.txt', 'w') as f:
            writer = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
            #if len(cl_name) == 1: # case if only one target to extract:
            if ext == 0:
                    writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n]])
            else:
                # TODO: see if we can simplify the following lines here
                if len(cl_name) == 1:
                    os.system(f'cp {cl_extractreg} {cluster}')
                    cl_extractreg = os.path.basename(cl_extractreg)
                else:
                    os.system(f'cp {cl_extractreg[n]} {cluster}')
                    cl_extractreg[n] = os.path.basename(cl_extractreg[n])

                if mreg==0 and subreg==0:
                    if len(cl_name) == 1:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg])
                    else:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n]])
                elif mreg==1 and subreg==0:
                    os.system(f'cp {cl_maskreg} {cluster}')
                    if len(cl_name) == 1:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg, cl_maskreg])
                    else:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n], cl_maskreg[n]])
                elif mreg==0 and subreg==1:
                    os.system(f'cp {subtract_reg_file} {cluster}')
                    if len(cl_name) == 1:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg, 'None', subtract_reg_file])
                    else:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n], 'None', cl_subreg[n]])
                else:
                    os.system(f'cp {subtract_reg_file} {cluster}')
                    os.system(f'cp {cl_maskreg} {cluster}')
                    if len(cl_name) == 1:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg, cl_maskreg, subtract_reg_file])
                    else:
                        writer.writerow([cl_z[n], cl_ra[n], cl_dec[n], cl_name[n], cl_extractreg[n], cl_maskreg[n], cl_subreg[n]])

        os.chdir(str(cluster))
        if no_selfcal:
            os.system(f'LOFAR_extract.py -p {pathdir} --beamcut {beam_cut} --extreg {cl_extractreg} --maskreg {userReg} --ampcal {ampcal} --ampsol {ampsol} --phsol {phsol} --maxniter {maxniter} --subreg {subtract_reg_file} --no_selfcal')
        else:
            os.system(f'LOFAR_extract.py -p {pathdir} --beamcut {beam_cut} --extreg {cl_extractreg} --maskreg {userReg} --ampcal {ampcal} --ampsol {ampsol} --phsol {phsol} --maxniter {int(maxniter)} --subreg {subtract_reg_file}')
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
                raise lib_util.Exit
            else:
                logger.info(f'Target {cluster} has been extracted.')
                os.chdir('../')
        size = float(os.path.getsize(str(cluster))/1e+9) #get directory size in GB, delete if too small?


if len(cl_failed) < len(cl_name): # if all or some are successfull
    logger.info('Concluded. Extracted datasets are in the mss-extract directory of each target.')
    logger.info('A new column DIFFUSE_SUB has been added to each successfully extracted .ms file.')
    logger.info('It contains the compact-source-subtracted visibilities.')
    logger.info('Nominal, high and source-subtracted low resolution images are in the /img directory of each target.')
else: # none worked :(
    logger.error(f'Extraction of all target(s): {",".join(cl_failed)} failed')
if (len(cl_failed) < len(cl_name)) and len(cl_failed): # case some but not all failed
    logger.warning(f'Extraction of target(s): {",".join(cl_failed)} failed')

os.system('rm -r logs_target_extraction.logger*') # directory not used in target_extract
#os.system('rm cluster_extraction.logger*')
