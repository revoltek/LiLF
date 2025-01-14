#!/usr/bin/env python3

import os, sys, glob, getpass, socket, re, pickle
from LiLF.surveys_db import SurveysDB
from LiLF import lib_ms, lib_img, lib_util, lib_log
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

caldirroot = ('iranet/groups/ulu/fdg/surveycals/done/')
tgtdirroot = ('/iranet/groups/ulu/fdg/surveytgts/download*/mss/')

#def calibrator_tables_available(obsid):
#    """
#    check if calibrator data exist in the database
#    """
#    with SurveysDB(survey='lba',readonly=True) as sdb:
#        sdb.execute('SELECT * FROM observations WHERE id=%f' % obsid)
#        r = sdb.cur.fetchall()
#        if len(r) != 0 and r[0]['location'] != '': return True
#        else: return False


#def local_calibrator_dirs(searchdir='', obsid=None):
#    """
#    Return the dirname of the calibrators
#    """
#    if searchdir != '': searchdir += '/'
#    if obsid is None:
#        calibrators = glob.glob(searchdir+'id*_-_*3[C|c]196*') + \
#                  glob.glob(searchdir+'id*_-_*3[C|c]295*') + \
#                  glob.glob(searchdir+'id*_-_*3[C|c]380*')
#    else:
#        calibrators = glob.glob(searchdir+'/id%i_-_*3[C|c]196*' % obsid) + \
#                  glob.glob(searchdir+'id%i_-_*3[C|c]295*' % obsid) + \
#                  glob.glob(searchdir+'id%i_-_*3[C|c]380*' % obsid)
#    if len(calibrators) == 0: return []
#    else: return calibrators


def update_status_db(field, status):
    with SurveysDB(survey='lba',readonly=False) as sdb:
        r = sdb.execute('UPDATE fields SET status="%s" WHERE id="%s"' % (status,field))

def check_done(pipename):
    """
    check if "Done" is written in the last line of the log file, otherwise quit with error.
    """
    logfile = sorted(glob.glob(pipename+'_*.logger'))[-1]
    with open(logfile, 'r') as f:
        last_line = f.readlines()[-1]
    if not "Done" in last_line:
        update_status_db(target, 'Error') 
        logger.error('Something went wrong in the last pipeline call.')
        sys.exit()

####################################################################################
# the the target from the db
logger.info('### Quering database...')
with SurveysDB(survey='lba',readonly=True) as sdb:
    sdb.execute('SELECT * FROM fields WHERE status="Downloaded" order by priority asc')
    r = sdb.cur.fetchall()
    if len(r) == 0:
        logger.warning('No field left in the db...')
        sys.exit()
    target = r[0]['id'] # here we set $target
    target_ra = r[0]['ra']
    target_dec = r[0]['decl']
    # save target name
    with open("target.txt", "w") as file:
        print('%s,%f,%f' % (target,target_ra,target_dec), file=file)

    sdb.execute('SELECT * FROM field_obs WHERE field_id="%s"' % target)
    r = sdb.cur.fetchall()
    obsids = [x['obs_id'] for x in r] # here we set $obsid

# add info, like cluster, node, user to the db...
logger.info(f"### Working on target: {target} (obsids: {obsids})")
username = getpass.getuser()
clustername = s.cluster
nodename = socket.gethostname()
with SurveysDB(survey='lba',readonly=False) as sdb:
    r = sdb.execute('UPDATE fields SET username="%s" WHERE id="%s"' % (username, target))
    r = sdb.execute('UPDATE fields SET clustername="%s" WHERE id="%s"' % (clustername, target))
    r = sdb.execute('UPDATE fields SET nodename="%s" WHERE id="%s"' % (nodename, target))

###################################################################################
# setup and copy

update_status_db(target, 'Copy') 

if not os.path.exists(working_dir):
    os.makedirs(working_dir)
if os.path.exists('lilf.config') and os.getcwd() != working_dir: 
    os.system('cp lilf.config '+working_dir)

os.chdir(working_dir)

# copy data
with w.if_todo('copy_data'):
            
    # get the targets
    for obsid in obsids:
        tgtdir = glob.glob(f'{tgtdirroot}/id{obsid}*{target}*')
        if len(tgtdir) == 0:
            logger.error(f'Missing target data for target: {target} (id: {obsid}).')
            sys.exit()
        logger.info(f'Found target data:', tgtdir[0])
        os.system('cp -r {tgtdir[0]} .')

    # now get the cals
    for obsid in obsids:
        caldir = glob.glob(f'{caldirroot}/id{obsid}*')
        if len(tgtdir) != 1:
            logger.error(f'Missing or too many calibrator data for target: {target} (id: {obsid}).')
            sys.exit()
        logger.info(f'Found calibrator data:', caldir[0])
        os.system('cp -r {caldir[0]} .')
### DONE

target_dirs = glob.glob('id*_-_*'+target+'*')
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
    for i, tc in enumerate(sorted(glob.glob('../id*_-_*'+target+'*/mss/TC*MS'))):
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
        logger.info(f'Self residual rms noise (cycle 0): %.1f mJy/b' % (qdict["self_c0_rms"] * 1e3))
        logger.info(f'Self residual rms noise (cycle 1): %.1f mJy/b' % (qdict["self_c1_rms"] * 1e3))
        logger.info('DDcal residual rms noise (cycle 0): %.1f mJy/b' % (qdict['ddcal_c0_rms'] * 1e3))
        logger.info('DDcal residual rms noise (cycle 1): %.1f mJy/b' % (qdict['ddcal_c1_rms'] * 1e3))
        logger.info('DDcal NVSS ratio (cycle 1): %.1f with %i matches' % (qdict['nvss_ratio'], qdict['nvss_match']))
        logger.info('DDcal total flags: %.1f%%' % (qdict['flag_frac']*100))

        with SurveysDB(survey='lba', readonly=False) as sdb:
            r = sdb.execute('UPDATE fields SET noise="%s", nvss_ratio="%s", nvss_match="%s", flag_frac="%s" WHERE id="%s"' \
                    % (qdict['ddcal_c1_rms'],qdict['nvss_ratio'], qdict['nvss_match'], qdict['flag_frac'],  target))
### DONE

################################################################################
# Save products
update_status_db(target, 'SaveProducts')

with w.if_todo('saveproducts_%s' % target):
    # copy images in herts
    logger.info('Copy ddcal products -> lofar.herts.ac.uk:/beegfs/lofar/lba/products/%s' % target)
    os.system('ssh herts "rm -rf /beegfs/lofar/lba/products/%s"' % target)
    os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s"' % target)
    os.system('scp -q self/images/wideP*.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q self/images/wideM-1-MFS-image.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q self/images/wide-largescale-MFS-image.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s/plots"' % target)
    os.system('scp -q -r self/plots/* herts:/beegfs/lofar/lba/products/%s/plots' % target)
    os.system('scp -q ddcal/c0*/images/wideDD-c*.MFS-image.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q ddcal/c0*/images/wideDD-c*.MFS-image-pb.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q ddcal/c0*/images/wideDD-c*.MFS-residual.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q ddcal/c01/solutions/interp.h5 herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q ddcal/c0*/skymodels/all*reg herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('scp -q ddcal/primarybeam.fits herts:/beegfs/lofar/lba/products/%s' % target)
    os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s/init"' % target)
    os.system('scp -q -r ddcal/init/*model.fits ddcal/init/wideM-1-sources.txt herts:/beegfs/lofar/lba/products/%s/init' % target)
    # logs
    os.system('ssh herts "mkdir /beegfs/lofar/lba/products/%s/logs"' % target)
    os.system('scp -q ../*logger ../*walker ../*%s*/*logger ../*%s*/*walker herts:/beegfs/lofar/lba/products/%s/logs' % (target, target, target))
    # copy ms in Bologna
    logger.info('Copy mss -> pleiadi:/iranet/lofarfs2/lofar2/fdg/surveytgts/%s' % target)
    os.system('tar zcf %s.tgz mss-avg' % target)
    os.system('ssh pleiadi "rm -rf /iranet/lofarfs2/lofar2/fdg/surveytgts/%s"' % target)
    os.system('ssh pleiadi "mkdir -p /iranet/lofarfs2/lofar2/fdg/surveytgts/%s"' % target)
    os.system('scp -q %s.tgz pleiadi:/iranet/lofarfs2/lofar2/fdg/surveytgts/%s' % (target,target))
### DONE


update_status_db(target, 'Done')
logger.info('### %s: Done. #####################################' % target)
w.alldone()
