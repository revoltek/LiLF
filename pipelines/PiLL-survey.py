#!/usr/bin/env python3

import os, sys, glob, getpass, socket, pickle, random
from datetime import datetime
from LiLF.surveys_db import SurveysDB
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
project = parset.get('PiLL','project')
target = parset.get('PiLL','target')
obsid = parset.get('PiLL','obsid')
download_file = parset.get('PiLL','download_file')

caldirroot = ('/iranet/groups/ulu/fdg/surveycals/done/')
tgtdirroot = ('/iranet/groups/ulu/fdg/surveytgts/download*/mss/')


def update_status_db(field, status):

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with SurveysDB(survey='lba',readonly=False) as sdb:
        r = sdb.execute('UPDATE fields SET end_date="%s", status="%s" WHERE id="%s"' % (timestamp, status, field))

def check_done(pipename):
    """
    check if "Done" is written in the last line of the log file, otherwise quit with error.
    """
    logfile = sorted(glob.glob(pipename+'_*.logger'))[-1]
    with open(logfile, 'r') as f:
        last_line = f.readlines()[-1]
    if not "Done" in last_line:
        update_status_db(target, 'Error') 

        # save everything possible for debug
        with w.if_todo('saveproduct_error_%s' % target):
            archive = '/iranet/groups/ulu/fdg/surveytgts/done/'+target+'.error'
            # copy images in herts
            logger.info(f'Copy products -> {archive}')
            lib_util.check_rm(f'{archive}')
            os.system(f'mkdir {archive}; mkdir {archive}/plots {archive}/logs')
            os.chdir(working_dir+'/'+target)
            os.system(f'cp ddparallel/images/wideM-*-MFS-image.fits {archive}')
            os.system(f'cp ddparallel/images/wide-lr-MFS-image.fits {archive}')
            os.system(f'cp ddparallel/solutions/faceets*reg {archive}')
            os.system(f'cp ddparallel/skymodel/subfield.reg {archive}')
            os.system(f'cp -r ddparallel/plots/* {archive}/plots')
            os.system(f'cp ddserial/c0*/images/*MFS-image.fits {archive}')
            os.system(f'cp ddserial/c0*/images/wideDD-*MFS-residual.fits {archive}')
            os.system(f'cp ddserial/c0*/solutions/facets-c*.reg {archive}')
            os.system(f'cp ddserial/c0*/skymodels/all*reg {archive}')
            os.system(f'cp ddserial/c0*/skymodels/mask-ddcal-c*.cat.fits {archive}')
            os.system(f'cp ddserial/c0*/skymodels/mask-ddcal-c*.reg {archive}')
            os.system(f'cp quality/quality.pickle {archive}')
            # copy logs
            logger.info(f'Copy logs -> {archive}')
            os.chdir(working_dir)
            os.system(f'cp -r PiLL_*logger PiLL*walker \
              *{target[:-1]}*/pipeline-timesplit_*logger *{target[:-1]}*/pipeline-timesplit.walker *{target[:-1]}*/logs_pipeline-timesplit_* \
              {target}/pipeline-ddparallel_*logger {target}/pipeline-ddparallel.walker {target}/logs_pipeline-ddparallel_* \
              {target}*/pipeline-ddserial_*logger {target}/pipeline-ddserial.walker {target}/logs_pipeline-ddserial_* \
              {archive}/logs')

        logger.error('Something went wrong in the last pipeline call.')
        sys.exit()

####################################################################################
# the the target from the db
logger.info('### Quering database...')
with SurveysDB(survey='lba',readonly=True) as sdb:
    if os.path.exists('target.txt'):
         with open("target.txt", "r") as file:
            target = file.readline()[:-1]
    else:
        # get all fields with max priority
        sdb.execute('SELECT * FROM fields WHERE status = "Downloaded" AND priority = (SELECT MAX(priority) FROM fields WHERE status = "Downloaded")')
        r = sdb.cur.fetchall()
        if len(r) == 0:
            logger.warning('No field left in the db...')
            sys.exit()
        
        rndidx = random.randint(0, len(r)-1) # select a random field
        target = r[rndidx]['id'] # here we set $target
        target_ra = r[rndidx]['ra']
        target_dec = r[rndidx]['decl']
        # save target name
        with open("target.txt", "w") as file:
            print(target, file=file)

    sdb.execute('SELECT * FROM field_obs WHERE field_id="%s"' % target)
    r = sdb.cur.fetchall()
    obsids = [x['obs_id'] for x in r] # here we set $obsid

# add info, like cluster, node, user to the db...
logger.info(f"### Working on target: {target} (obsids: {obsids})")
username = getpass.getuser()
clustername = s.cluster
nodename = socket.gethostname()
timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
with SurveysDB(survey='lba',readonly=False) as sdb:
    r = sdb.execute('UPDATE fields SET username="%s", clustername="%s", nodename="%s", start_date="%s" WHERE id="%s"' % \
                    (username, clustername, nodename, target, timestamp))

###################################################################################
# setup and copy

update_status_db(target, 'Copy') 

# pipelines search for lilf.config also in the ../ dir
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
if os.path.exists('lilf.config') and os.getcwd() != working_dir: 
    os.system('cp lilf.config '+working_dir)

os.chdir(working_dir)

# copy data
with w.if_todo('copy_data'):
            
    # get the targets
    for obsid in obsids:
        tgtdir = glob.glob(f'{tgtdirroot}/id{obsid}*{target[:-1]}*')
        if len(tgtdir) == 0:
            logger.error(f'Missing target data for target: {target} (id: {obsid}).')
            sys.exit()
        logger.info(f'Copy target data: {tgtdir[0]}')
        os.system(f'cp -r {tgtdir[0]} .')

    # now get the cals
    for obsid in obsids:
        caldir = glob.glob(f'{caldirroot}/id{obsid}*')
        if len(tgtdir) != 1:
            logger.error(f'Missing or too many calibrator data for target: {target} (id: {obsid}).')
            sys.exit()
        logger.info(f'Copy calibrator data: {caldir[0]}')
        os.system(f'cp -r {caldir[0]} .')
### DONE

target_dirs = glob.glob('id*_-_*'+target[:-1]+'*')
logger.debug('TARGET DIRS: %s' % (','.join(target_dirs) ) )

################################################################################
# timesplit

update_status_db(target, 'Timesplit') 

for target_dir in target_dirs:

    with w.if_todo('timesplit_%s' % target_dir):
        logger.info('### %s: Starting timesplit... #####################################' % target_dir)
        os.chdir(working_dir+'/'+target_dir)
        if not os.path.exists('data-bkp'):
            os.makedirs('data-bkp')
            os.system('mv *tar data-bkp')

        os.system(LiLF_dir+'/pipelines/LOFAR_timesplit.py')
        check_done('pipeline-timesplit')
    ### DONE

if not os.path.exists(working_dir+'/'+target):
    os.makedirs(working_dir+'/'+target)
os.chdir(working_dir+'/'+target)
    
# collet mss for this target
if not os.path.exists('mss'):
    os.makedirs('mss')
    for i, tc in enumerate(sorted(glob.glob('../id*_-_*'+target[:-1]+'*/mss/TC*MS'))):
        tc_ren = 'TC%02i.MS' % i
        logger.debug('mv %s mss/%s' % (tc,tc_ren))
        os.system('mv %s mss/%s' % (tc,tc_ren))

################################################################################
# selfcal

# DD-parallel
update_status_db(target, 'DD-parallel')

with w.if_todo('dd-parallel_%s' % target):
    logger.info('### %s: Starting dd-parallel #####################################' % target)
    os.system(LiLF_dir+'/pipelines/LOFAR_ddparallel.py')
    check_done('pipeline-ddparallel')
### DONE

# DD-serial
update_status_db(target, 'DD-serial')

with w.if_todo('dd-serial_%s' % target):
    logger.info('### %s: Starting dd-serial #####################################' % target)
    os.system(LiLF_dir+'/pipelines/LOFAR_ddserial.py')
    check_done('pipeline-ddserial')
### DONE

################################################################################
# Quality check
update_status_db(target, 'QualityCheck')

with w.if_todo('quality_%s' % target):
    logger.info('### %s: Starting quality check #####################################' % target)
    os.system(LiLF_dir+'/pipelines/LOFAR_quality.py')
    check_done('pipeline-quality')

    with open('quality/quality.pickle', 'rb') as f:
        qdict = pickle.load(f)
    logger.info('DDparallel residual rms noise (cycle 0): %.1f mJy/b' % (qdict["ddparallel_c0_rms"] * 1e3))
    logger.info('DDparallel residual rms noise (cycle 1): %.1f mJy/b' % (qdict["ddparallel_c1_rms"] * 1e3))
    logger.info('DDserial residual rms noise (cycle 0): %.1f mJy/b' % (qdict['ddserial_c0_rms'] * 1e3))
    #logger.info('DDserial residual rms noise (cycle 1): %.1f mJy/b' % (qdict['ddserial_c1_rms'] * 1e3))
    logger.info('DDserial NVSS ratio (cycle 1): %.1f with %i matches' % (qdict['nvss_ratio'], qdict['nvss_match']))
    logger.info('DDserial total flags: %.1f%%' % (qdict['flag_frac']*100))

    with SurveysDB(survey='lba', readonly=False) as sdb:
        r = sdb.execute('UPDATE fields SET noise="%s", nvss_ratio="%s", nvss_match="%s", flag_frac="%s" WHERE id="%s"' \
            % (qdict['ddserial_c0_rms'],qdict['nvss_ratio'], qdict['nvss_match'], qdict['flag_frac'],  target))
### DONE

################################################################################
# Save products
update_status_db(target, 'SaveProducts')

# on pleiadi
with w.if_todo('saveproducts_%s' % target):
    archive = '/iranet/groups/ulu/fdg/surveytgts/done/'+target
    # copy images in herts
    logger.info(f'Copy products -> {archive}')
    lib_util.check_rm(f'{archive}')
    os.system(f'mkdir {archive}; mkdir {archive}/plots {archive}/logs')
    os.system(f'cp ddparallel/images/wideM-*-MFS-image.fits {archive}')
    os.system(f'cp ddparallel/images/wide-lr-MFS-image.fits {archive}')
    os.system(f'cp ddparallel/solutions/facets*reg {archive}')
    os.system(f'cp ddparallel/skymodel/subfield.reg {archive}')
    os.system(f'cp -r ddparallel/plots/* {archive}/plots')
    os.system(f'cp ddserial/c0*/images/*image*.fits {archive}')
    os.system(f'cp ddserial/c0*/images/wideDD-*MFS-residual.fits {archive}')
    os.system(f'cp ddserial/c00/images/wideDD*model*fpb.fits {archive}')
    os.system(f'gzip ddserial/c00/solutions/interp.h5; cp ddserial/c00/solutions/interp.h5.gz {archive}')
    os.system(f'cp ddserial/c00/solutions/facets-c00.reg {archive}')
    os.system(f'cp ddserial/c0*/skymodels/all*reg {archive}')
    os.system(f'cp ddserial/c0*/skymodels/mask-ddcal-c*.cat.fits {archive}')
    os.system(f'cp ddserial/c0*/skymodels/mask-ddcal-c*.reg {archive}')
    os.system(f'cp ddserial/primarybeam.fits {archive}')
    os.system(f'cp quality/quality.pickle {archive}')
    # copy ms
    logger.info(f'Copy mss -> {archive}')
    os.system(f'tar zcf {target}.tgz mss-avg')
    os.system(f'cp {target}.tgz {archive}')
    # copy logs
    logger.info(f'Copy logs -> {archive}')
    os.chdir(working_dir)
    os.system(f'cp -r PiLL_*logger PiLL*walker \
              *{target[:-1]}*/pipeline-timesplit_*logger *{target[:-1]}*/pipeline-timesplit.walker *{target[:-1]}*/logs_pipeline-timesplit_* \
              {target}/pipeline-ddparallel_*logger {target}/pipeline-ddparallel.walker {target}/logs_pipeline-ddparallel_* \
              {target}*/pipeline-ddserial_*logger {target}/pipeline-ddserial.walker {target}/logs_pipeline-ddserial_* \
              {archive}/logs')
### DONE

update_status_db(target, 'Done')
logger.info('### %s: Done. #####################################' % target)
w.alldone()
