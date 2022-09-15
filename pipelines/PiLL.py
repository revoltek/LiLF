#!/usr/bin/env python3

import os, sys, glob, getpass, socket, re, pickle
from LiLF.surveys_db import SurveysDB
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('PiLL.logger')
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

def calibrator_tables_available(obsid):
    """
    check if calibrator data exist in the database
    """
    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('SELECT * FROM observations WHERE id=%f' % obsid)
        r = sdb.cur.fetchall()
        if len(r) != 0 and r[0]['location'] != '': return True
        else: return False


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


def update_status_db(field, status):
    with SurveysDB(survey='lba',readonly=False) as sdb:
        r = sdb.execute('UPDATE fields SET status="%s" WHERE id="%s"' % (status,field))


def check_done(logfile):
    """
    check if "Done" is written in the last line of the log file, otherwise quit with error.
    """
    with open(logfile, 'r') as f:
        last_line = f.readlines()[-1]
    if not "Done" in last_line:
        if survey: update_status_db(target, 'Error') 
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


####################################################################################

# If no input given, assume survey data and query the database for data to process
survey = False
if download_file == '' and project == '' and target == '' and obsid == '':
    survey = True
    if os.path.exists('target.txt'):
        with open('target.txt', 'r') as file:
                target = file.read().replace('\n', '')
    else:
        logger.info('### Quering database...')
        with SurveysDB(survey='lba',readonly=True) as sdb:
            sdb.execute('SELECT * FROM fields WHERE status="Observed" order by priority asc')
            r = sdb.cur.fetchall()
            target = r[0]['id'] # here we set $target
        # save target name
        with open("target.txt", "w") as file:
                print(target, file=file)

    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('SELECT * FROM field_obs WHERE field_id="%s"' % target)
        r = sdb.cur.fetchall()
        obsid = ','.join([str(x['obs_id']) for x in r]) # here we set $obsid

    logger.info("### Working on target: %s (obsid: %s)" % (target, obsid))
    # add other info, like cluster, node, user...
    username = getpass.getuser()
    clustername = s.cluster
    nodename = socket.gethostname()
    with SurveysDB(survey='lba',readonly=False) as sdb:
        r = sdb.execute('UPDATE fields SET username="%s" WHERE id="%s"' % (username, target))
        r = sdb.execute('UPDATE fields SET clustername="%s" WHERE id="%s"' % (clustername, target))
        r = sdb.execute('UPDATE fields SET nodename="%s" WHERE id="%s"' % (nodename, target))

if survey: update_status_db(target, 'Download') 

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
    check_done('pipeline-preprocess.logger')
    os.system('mv mss/* ../')
### DONE

os.chdir(working_dir)
if target != '': fix_dirs(target) # try to harmonize target dirs
target_dirs = [t for t in glob.glob('id*') if t not in local_calibrator_dirs()]
logger.debug('TARGET DIRS: %s' % (','.join(target_dirs) ) )

if survey: update_status_db(target, 'Calibrator')

for target_dir in target_dirs:
    ##########
    # calibrator
    # here the pipeline checks if the calibrator is available online, otherwise it downloads it
    # then it also runs the calibrator pipeline
    obsid = int(target_dir.split('_-_')[0][2:])
    with w.if_todo('cal_id%i' % obsid):
        if not survey or redo_cal or not calibrator_tables_available(obsid):
            logger.info('### %s: Starting calibrator... #####################################' % target_dir)
            cal_dir = local_calibrator_dirs(working_dir, obsid)
            # if calibrator not downloaded, do it
            if len(cal_dir) == 0:
                if not os.path.exists(working_dir+'/download-cal_id%i' % obsid):
                    os.makedirs(working_dir+'/download-cal_id%i' % obsid)
                os.chdir(working_dir+'/download-cal_id%i' % obsid)
                os.system(LiLF_dir+'/scripts/LOFAR_stager.py --quiet --cal --obsID %i >> stager.log' % (obsid))
                os.system(LiLF_dir+'/pipelines/LOFAR_preprocess.py')
                check_done('pipeline-preprocess.logger')
                os.system('mv mss/* ../')

            os.chdir(local_calibrator_dirs(working_dir, obsid)[0])
            if not os.path.exists('data-bkp'):
                os.makedirs('data-bkp')
                os.system('mv *MS data-bkp')
            os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')
            check_done('pipeline-cal.logger')

            if survey: # only backup solutions if survey
                # copy solutions in the repository
                cal_dir = os.path.basename(local_calibrator_dirs(working_dir, obsid)[0])
                logger.info('Copy: cal*h5 -> herts:/beegfs/lofar/lba/calibration_solutions/%s' % cal_dir)
                os.system('ssh herts "rm -rf /beegfs/lofar/lba/calibration_solutions/%s"' % cal_dir)
                os.system('ssh herts "mkdir /beegfs/lofar/lba/calibration_solutions/%s"' % cal_dir)
                os.system('scp -q cal-pa.h5 cal-amp.h5 cal-iono.h5 herts:/beegfs/lofar/lba/calibration_solutions/%s' % cal_dir)
                os.system('scp -q -r plots* herts:/beegfs/lofar/lba/calibration_solutions/%s' % cal_dir)

                # update the db
                with SurveysDB(survey='lba',readonly=False) as sdb:
                    sdb.execute('INSERT INTO observations (id,location,calibratordata) VALUES \
                    (%i,"herts","%s")' % (obsid, "/beegfs/lofar/lba/calibration_solutions/"+cal_dir))

        else:
            # download calibrator solutions
            with SurveysDB(survey='lba',readonly=True) as sdb:
                sdb.execute('select location,calibratordata from observations where id=%i' % obsid)
                r = sdb.cur.fetchall()[0]
            location = r['location']
            calibratordata = r['calibratordata']
            os.chdir(working_dir)
            cal_dir = calibratordata.split('/')[-1]
            logger.info('Downloading solutions: %s:%s/cal-*h5' % (location,calibratordata))
            os.system('mkdir %s' % cal_dir)
            os.system('scp -q %s:%s/cal-*h5 %s' % (location,calibratordata,cal_dir))
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
        check_done('pipeline-timesplit.logger')
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
        if survey: update_status_db(grouped_target, 'Self')
        logger.info('### %s: Starting selfcal #####################################' % grouped_target)
        os.system(LiLF_dir+'/pipelines/LOFAR_self.py')
        check_done('pipeline-self.logger')
    ### DONE

    # DD-cal
    with w.if_todo('dd_%s' % grouped_target):
        if survey: update_status_db(grouped_target, 'Ddcal')
        logger.info('### %s: Starting ddcal #####################################' % grouped_target)
        os.system(LiLF_dir+'/pipelines/LOFAR_dd-serial.py')
        check_done('pipeline-dd-serial.logger')

        if survey: # only back up solutions if survey
            logger.info('Copy: ddcal/c0*/images/img/wideDD-c*... -> lofar.herts.ac.uk:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('ssh herts "rm -rf /beegfs/lofar/lba/products/%s"' % grouped_target)
            os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s"' % grouped_target)
            os.system('scp -q self/images/wideP*.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q self/images/wideM-1-MFS-image.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q self/images/wide-largescale-MFS-image.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q ddcal/c0*/images/wideDD-c*.app.restored.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q ddcal/c0*/images/wideDD-c*.int.restored.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q '+sorted(glob.glob('ddcal/c00/images/wideDD-c00.residual*.fits'))[-1]+' herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q '+sorted(glob.glob('ddcal/c01/images/wideDD-c01.residual*.fits'))[-1]+' herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q ddcal/c01/solutions/interp.h5 herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('scp -q ddcal/c0*/skymodels/all*reg herts:/beegfs/lofar/lba/products/%s' % grouped_target)
            os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s/plots"' % grouped_target)
            os.system('scp -q ddcal/c0*/plots/* herts:/beegfs/lofar/lba/products/%s/plots' % grouped_target)
    ### DONE

    # Quality check
    with w.if_todo('quality_%s' % grouped_target):
        if survey: update_status_db(grouped_target, 'QualityCheck')
        logger.info('### %s: Starting quality check #####################################' % grouped_target)
        os.system(LiLF_dir+'/pipelines/LOFAR_quality.py')
        check_done('pipeline-quality.logger')

        with open('quality.pickle', 'rb') as f:
            qdict = pickle.load(f)
        logger.info(f'Self residual rms noise (cycle 0): %.1f mJy/b' % (qdict["self_c0_rms"] * 1e3))
        logger.info(f'Self residual rms noise (cycle 1): %.1f mJy/b' % (qdict["self_c1_rms"] * 1e3))
        logger.info('DDcal residual rms noise (cycle 0): %.1f mJy/b' % (qdict['ddcal_c0_rms'] * 1e3))
        logger.info('DDcal residual rms noise (cycle 1): %.1f mJy/b' % (qdict['ddcal_c1_rms'] * 1e3))
        logger.info('DDcal NVSS ratio (cycle 1): %.1f mJy/b' % (qdict['nvss_ratio'] * 1e3))
        if survey:
            with SurveysDB(survey='lba', readonly=False) as sdb:
                r = sdb.execute('UPDATE fields SET noise="%s", nvss_ratio="%s" WHERE id="%s"' % (qdict['ddcal_c1_rms'],qdict['nvss_ratio'], grouped_target)) # remove upper?
    ### DONE

    if survey: update_status_db(grouped_target, 'Done')
    logger.info('### %s: Done. #####################################' % grouped_target)
