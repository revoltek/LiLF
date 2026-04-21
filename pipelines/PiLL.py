#!/usr/bin/env python3

import os, sys, glob, getpass, socket, pickle
from LiLF import lib_util, lib_log
logger_obj = lib_log.Logger('PiLL')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('PiLL.walker')

LiLF_dir = os.path.dirname(os.path.dirname(lib_util.__file__))
parset = lib_util.getParset(parsetFile='lilf.config')

# get parameters
# use lilf.config (this is also used by all other scripits)
working_dir = os.path.abspath(parset.get('PiLL','working_dir'))
redo_cal = parset.getboolean('PiLL','redo_cal')
project = parset.get('PiLL','project')
target = parset.get('PiLL','target')
obsid = parset.get('PiLL','obsid')
download_file = parset.get('PiLL','download_file')
if download_file != '': download_file = os.path.abspath(download_file)

def local_calibrator_dirs(searchdir='', obsid=None):
    """
    Return the dirname of the calibrators
    """
    if searchdir != '': searchdir += '/'
    if obsid is None:
        calibrators = glob.glob(searchdir+'id*_-_*3[C|c]196*') + \
                  glob.glob(searchdir+'id*_-_*3[C|c]295*') + \
                  glob.glob(searchdir+'id*_-_*3[C|c]380*')
    else:
        calibrators = glob.glob(searchdir+'/id%i_-_*3[C|c]196*' % obsid) + \
                  glob.glob(searchdir+'id%i_-_*3[C|c]295*' % obsid) + \
                  glob.glob(searchdir+'id%i_-_*3[C|c]380*' % obsid)
    if len(calibrators) == 0: return []
    else: return calibrators


def check_done(pipename):
    """
    check if "Done" is written in the last line of the log file, otherwise quit with error.
    """
    logfile = sorted(glob.glob(pipename+'_*.logger'))[-1]
    with open(logfile, 'r') as f:
        last_line = f.readlines()[-1]
    if not "Done" in last_line:
        logger.error('Something went wrong in the last pipeline call.')
        sys.exit()


def fix_dirs(target):
    """
    fix dir names in the working_dir
    """
    calibrators = local_calibrator_dirs(working_dir)
    targets = [t for t in glob.glob(working_dir+'/id*') if t not in calibrators]
    # rename targets to consistent format
    for t in targets:
        idx = t.split('_')[0]
        if t != '%s_-_%s' % (idx,target):
            os.system('mv %s %s_-_%s' % (t,idx,target))
    # rename calibrators to consistent format
    for c in calibrators:
        idx = c.split('_')[0]
        cal = c.split('_')[-1]
        if c != '%s_-_%s' % (idx,cal):
            os.system('mv %s %s_-_%s' % (c,idx,cal))


#######
# setup
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
if os.path.exists('lilf.config') and os.getcwd() != working_dir: 
    os.system('cp lilf.config '+working_dir)

os.chdir(working_dir)
if not os.path.exists(working_dir+'/download'):
    os.makedirs(working_dir+'/download')

if download_file != '':
    os.system('cp %s download/html.txt' % download_file)

##########
# data download
# here the pipeline downloads only the target, not the calibrator
with w.if_todo('download'):
    logger.info('### Starting download... #####################################')
    os.chdir(working_dir+'/download')

    if download_file == '':
        cmd = LiLF_dir+'/scripts/LOFAR_stager.py --quiet --nocal'
        if target != '':
            cmd += ' --target %s' % target
        if obsid != '':
            cmd += ' --obsID %s' % obsid
        logger.debug("Exec: %s" % cmd)
        os.system(cmd+' >> stager.log')

    logger.info('Downloaded %i files' % len(glob.glob('*MS')))
    os.system(LiLF_dir+'/pipelines/LOFAR_preprocess.py')
    check_done('pipeline-preprocess')
    os.system('mv mss/* ../')
### DONE

os.chdir(working_dir)
if target != '': fix_dirs(target) # try to harmonize target dirs
target_dirs = [t for t in glob.glob('id*') if t not in local_calibrator_dirs()]
logger.debug('TARGET DIRS: %s' % (','.join(target_dirs) ) )

for target_dir in target_dirs:
    ##########
    # calibrator
    # here the pipeline checks if the calibrator is available online, otherwise it downloads it
    # then it also runs the calibrator pipeline
    obsid = int(target_dir.split('_-_')[0][2:])
    with w.if_todo('cal_id%i' % obsid):
        if redo_cal:
            logger.info('### %s: Starting calibrator... #####################################' % target_dir)
            cal_dir = local_calibrator_dirs(working_dir, obsid)
            # if calibrator not downloaded, do it
            if len(cal_dir) == 0:
                if not os.path.exists(working_dir+'/download-cal_id%i' % obsid):
                    os.makedirs(working_dir+'/download-cal_id%i' % obsid)
                os.chdir(working_dir+'/download-cal_id%i' % obsid)
                os.system(LiLF_dir+'/scripts/LOFAR_stager.py --quiet --cal --obsID %i >> stager.log' % (obsid))
                os.system(LiLF_dir+'/pipelines/LOFAR_preprocess.py')
                check_done('pipeline-preprocess')
                os.system('mv mss/* ../')

            os.chdir(local_calibrator_dirs(working_dir, obsid)[0])
            if not os.path.exists('data-bkp'):
                os.makedirs('data-bkp')
                os.system('mv *MS data-bkp')
            os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')
            check_done('pipeline-cal')
    ### DONE

    ##########
    # timesplit
    # each target of each observation is then timesplit
    with w.if_todo('timesplit_%s' % target_dir):
        logger.info('### %s: Starting timesplit... #####################################' % target_dir)
        os.chdir(working_dir+'/'+target_dir)
        if not os.path.exists('data-bkp'):
            os.makedirs('data-bkp')
            os.system('mv *MS data-bkp')

        os.system(LiLF_dir+'/pipelines/LOFAR_timesplit.py')
        check_done('pipeline-timesplit')
    ### DONE

# group targets with same name, assuming they are different pointings in the same direction
os.chdir(working_dir)
target_dirs = [t for t in glob.glob('id*') if t not in local_calibrator_dirs()]
grouped_targets = set([t.split('_-_')[-1] for t in target_dirs])

for grouped_target in grouped_targets:
    if not os.path.exists(working_dir+'/'+grouped_target):
        os.makedirs(working_dir+'/'+grouped_target)
    os.chdir(working_dir+'/'+grouped_target)
    
    # collet mss for this grouped_target
    if not os.path.exists('mss'):
        os.makedirs('mss')
        for i, tc in enumerate(sorted(glob.glob('../id*_-_'+grouped_target+'/mss/TC*MS'))):
            tc_ren = 'TC%02i.MS' % i
            logger.debug('mv %s mss/%s' % (tc,tc_ren))
            os.system('mv %s mss/%s' % (tc,tc_ren))

    # selfcal
    with w.if_todo('self_%s' % grouped_target):
        logger.info('### %s: Starting selfcal #####################################' % grouped_target)
        os.system(LiLF_dir+'/pipelines/LOFAR_self.py')
        check_done('pipeline-self')
    ### DONE

    # DD-cal
    with w.if_todo('dd_%s' % grouped_target):
        logger.info('### %s: Starting ddcal #####################################' % grouped_target)
        os.system(LiLF_dir+'/pipelines/LOFAR_dd.py')
        check_done('pipeline-dd')
    ### DONE

    # Quality check
    with w.if_todo('quality_%s' % grouped_target):
        logger.info('### %s: Starting quality check #####################################' % grouped_target)
        os.system(LiLF_dir+'/pipelines/LOFAR_quality.py')
        check_done('pipeline-quality')

        with open('quality/quality.pickle', 'rb') as f:
            qdict = pickle.load(f)
        logger.info('Self residual rms noise (cycle 0): %.1f mJy/b' % (qdict["self_c0_rms"] * 1e3))
        logger.info('Self residual rms noise (cycle 1): %.1f mJy/b' % (qdict["self_c1_rms"] * 1e3))
        logger.info('DDcal residual rms noise (cycle 0): %.1f mJy/b' % (qdict['ddcal_c0_rms'] * 1e3))
        logger.info('DDcal residual rms noise (cycle 1): %.1f mJy/b' % (qdict['ddcal_c1_rms'] * 1e3))
        logger.info('DDcal NVSS ratio (cycle 1): %.1f with %i matches' % (qdict['nvss_ratio'], qdict['nvss_match']))
        logger.info('DDcal total flags: %.1f%%' % (qdict['flag_frac']*100))
     ### DONE

    logger.info('### %s: Done. #####################################' % grouped_target)
