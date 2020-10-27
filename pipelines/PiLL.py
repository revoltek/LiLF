#!/usr/bin/env python

import os, sys, glob, getpass, socket, re
from LiLF.surveys_db import SurveysDB
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('PiLL.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('PiLL.walker')

LiLF_dir = os.path.dirname(os.path.dirname(lib_util.__file__))
parset = lib_util.getParset(parsetFile='lilf.config')

survey_projects = 'LT14_002,LC12_017,LC9_016,LC8_031' # list of projects related with the LBA survey

# get parameters
# use lilf.config (this is also used by all other scripits)
working_dir = os.path.abspath(parset.get('PiLL','working_dir'))
redo_cal = parset.getboolean('PiLL','redo_cal')
project = parset.get('PiLL','project')
target = parset.get('PiLL','target')
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
        r = sdb.execute('UPDATE fields SET status="%s" WHERE id="%s"' % (status,field.upper()))


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

def fix_dir_format(working_dir):
    # fix for c##-o##_p##### format
    pattern = re.compile("^id[0-9]*_-_c[0-9][0-9]-o.*_.*$")
    for dir in glob.glob(working_dir+'/id*'):
        if pattern.match(dir):
            os.system('mv '+dir+' '+dir.split('_-_')[0]+'_-_'+dir.split('_')[-1])

####################################################################################

# query the database for data to process
survey = False
if download_file == '' and project == '' and target == '':
    survey = True
    project = survey_projects
    if os.path.exists('target.txt'):
        with open('target.txt', 'r') as file:
                target = file.read().replace('\n', '')
    else:
        logger.info('### Quering database...')
        with SurveysDB(survey='lba',readonly=True) as sdb:
            sdb.execute('SELECT * FROM fields WHERE status="Observed" order by priority desc')
            r = sdb.cur.fetchall()
            target = r[0]['id']
        # save target name
        with open("target.txt", "w") as file:
                print(target, file=file)

    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('SELECT * FROM field_obs WHERE field_id="%s"' % target)
        r = sdb.cur.fetchall()
        obsid = ','.join([str(x['obs_id']) for x in r])

    logger.info("### Working on target: %s (obsid: %s)" % (target, obsid))
    # add other info, like cluster, node, user...
    username = getpass.getuser()
    clustername = s.cluster
    nodename = socket.gethostname()
    with SurveysDB(survey='lba',readonly=False) as sdb:
        r = sdb.execute('UPDATE fields SET username="%s" WHERE id="%s"' % (username, target))
        r = sdb.execute('UPDATE fields SET clustername="%s" WHERE id="%s"' % (clustername, target))
        r = sdb.execute('UPDATE fields SET nodename="%s" WHERE id="%s"' % (nodename, target))

    update_status_db(target, 'Download') 

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
        cmd = LiLF_dir+'/scripts/LOFAR_stager.py --projects %s --nocal' % project
        if target != '':
            cmd += ' --target %s' % target
        if obsid != '':
            cmd += ' --obsID %s' % obsid
        logger.debug("Exec: %s" % cmd)
        os.system(cmd)

    # TODO: how to be sure all MS were downloaded?
    os.system(LiLF_dir+'/pipelines/LOFAR_download.py')
    os.system('mv mss/* ../')
    

### DONE

os.chdir(working_dir)
fix_dir_format(working_dir)
if survey: update_status_db(target, 'Calibrator')
calibrators = local_calibrator_dirs()
targets = [t for t in glob.glob('id*') if t not in calibrators]
logger.debug('CALIBRATORS: %s' % ( ','.join(calibrators) ) )
logger.debug('TARGET: %s' % (','.join(targets) ) )

for target in targets:
    
    ##########
    # calibrator
    # here the pipeline checks if the calibrator is available online, otherwise it downloads it
    # then it also runs the calibrator pipeline
    obsid = int(target.split('_-_')[0][2:])
    with w.if_todo('cal_id%i' % obsid):
        if redo_cal or not calibrator_tables_available(obsid):
            logger.info('### %s: Starting calibrator... #####################################' % target)
            # if calibrator not downaloaded, do it
            cal_dir = local_calibrator_dirs(working_dir, obsid)
            if len(cal_dir) == 0:
                if not os.path.exists(working_dir+'/download-cal_id%i' % obsid):
                    os.makedirs(working_dir+'/download-cal_id%i' % obsid)
                os.chdir(working_dir+'/download-cal_id%i' % obsid)
                os.system(LiLF_dir+'/scripts/LOFAR_stager.py --cal --projects %s --obsID %i' % (project, obsid))
                os.system(LiLF_dir+'/pipelines/LOFAR_download.py')
                check_done('pipeline-download.logger')
                os.system('mv mss/* ../')

            fix_dir_format(working_dir)
            os.chdir(local_calibrator_dirs(working_dir, obsid)[0])
            if not os.path.exists('data-bkp'):
                os.makedirs('data-bkp')
                os.system('mv *MS data-bkp')
            os.system(LiLF_dir+'/pipelines/LOFAR_cal.py')
            check_done('pipeline-cal.logger')

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
                sdb.execute('select (location, calibratordata) from observations where id=%i' % obsid)
                location = sdb.cur.fetchall()[0]['location']
                calibratordata = sdb.cur.fetchall()[0]['calibratordata']
                logger.info('Downloading solutions: %s:%s/cal-*h5' % (location,calibratordata))
                os.system('scp -q %s:%s/cal-*h5 .' % (location,calibratordata))

    ### DONE

    ##########
    # timesplit
    # each target of each observation is then timesplit
    with w.if_todo('timesplit_%s' % target):
        logger.info('### %s: Starting timesplit... #####################################' % target)
        os.chdir(working_dir+'/'+target)
        if not os.path.exists('data-bkp'):
            os.makedirs('data-bkp')
            os.system('mv *MS data-bkp')

        os.system(LiLF_dir+'/pipelines/LOFAR_timesplit.py')
        check_done('pipeline-timesplit.logger')


    ### DONE

# group targets with same name, assuming they are different pointings of the same dir
grouped_targets = set([t.split('_-_',1)[1] for t in targets])

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

        logger.info('Copy: ddcal/c0*/images/img/wideDD-c*... -> lofar.herts.ac.uk:/beegfs/lofar/lba/products/%s' % grouped_target)
        os.system('ssh herts "rm -rf /beegfs/lofar/lba/products/%s"' % grouped_target)
        os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s"' % grouped_target)
        os.system('scp -q ddcal/c0*/images/wideDD-c*.app.restored.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)
        os.system('scp -q ddcal/c0*/images/wideDD-c*.int.restored.fits herts:/beegfs/lofar/lba/products/%s' % grouped_target)

### DONE

    if survey: update_status_db(grouped_target, 'Done')
    logger.info('### %s: Done. #####################################' % grouped_target)
