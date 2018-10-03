#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation.
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import sys, os, glob, re
import numpy as np

if 'LBAsurvey' in os.getcwd():
    obs     = os.getcwd().split('/')[-2] # assumes .../c??-o??/3c196
    calname = os.getcwd().split('/')[-1] # assumes .../c??-o??/3c196
    data_dir = '../../download/%s/%s' % (obs, calname)
#    bl2flag = 'CS031LBA'
#    if 'c05' in os.getcwd(): bl2flag = 'RS409LBA'
#    if 'c09' in os.getcwd(): bl2flag = 'CS031LBA\;CS013LBA'

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
lib_log.set_logger('pipeline-cal.logger')
logger = lib_log.logger
s = lib_util.Scheduler(dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('cal','parset_dir')
data_dir = parset.get('cal','data_dir')
skymodel = parset.get('cal','skymodel')
imaging = parset.getboolean('cal','imaging')
bl2flag = parset.get('flag','stations')

#############################################################
MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s )
# copy data
logger.info('Copy data...')
for MS in MSs.getListObj():
    MS.move(MS.nameMS+'.MS', keepOrig=True)

# TEST HBA
#logger.warning('Copy data HBA...')
#MSs.run("DPPP DPPP-avg.parset msin=$pathMS msout=$nameMS.MS avg.freqstep=5", log="$nameMS_avg.log", commandType="DPPP")

MSs = lib_ms.AllMSs( glob.glob('*MS'), s )
calname = MSs.getListObj()[0].getNameField()
obsmode = MSs.getListObj()[0].getObsMode()
# find min freq
if min(MSs.getFreqs()) < 40.e6:
    iono3rd = True
    logger.debug('Include iono 3rd order.')
else: iono3rd = False

####################################################
## flag bad stations, flags will propagate
#logger.info("Flagging...")
#MSs.run("DPPP " + parset_dir + "/DPPP-flag.parset msin=$pathMS ant.baseline=\"" + bl2flag+"\"", log="$nameMS_flag.log", commandType="DPPP")
#
## predict to save time ms:MODEL_DATA
#logger.info('Add model to MODEL_DATA (%s)...' % calname)
#MSs.run("DPPP " + parset_dir + "/DPPP-predict.parset msin=$pathMS pre.sourcedb=" + skymodel + " pre.sources=" + calname, log="$nameMS_pre.log", commandType="DPPP")
#
####################################################
## 1: find PA
#
## Smooth data DATA -> SMOOTHED_DATA (BL-based smoothing)
#logger.info('BL-smooth...')
#MSs.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType ='python', maxThreads=20)
#
## Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#logger.info('Calibrating...')
#for MS in MSs.getListStr():
#    lib_util.check_rm(MS+'/pa.h5')
#MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/pa.h5 sol.mode=rotation+diagonal', log='$nameMS_solPA.log', commandType="DPPP")
#
#lib_util.run_losoto(s, 'pa', [ms+'/pa.h5' for ms in MSs.getListStr()], \
#        [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-pa.parset'])

# Pol align correction DATA -> CORRECTED_DATA
logger.info('Polalign correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DPPP")

####################################################
# 2: find FR

# Beam correction CORRECTED_DATA -> CORRECTED_DATA
logger.info('Beam correction...')
MSs.run("DPPP " + parset_dir + '/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType="DPPP")

# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
#logger.info('Converting to circular...')
#MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType ='python')

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/fr.h5')
MSs.run('DPPP ' + parset_dir + '/DPPP-soldd.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.mode=rotation+diagonal', log='$nameMS_solFR.log', commandType="DPPP")

lib_util.run_losoto(s, 'fr', [ms+'/fr.h5' for ms in MSs.getListStr()], \
        #[parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-fr.parset'])
        [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-rot.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-fr.parset'])

# Convert to linear CORRECTED_DATA -> CORRECTED_DATA
#logger.info('Converting to linear...')
#MSs.run('mslin2circ.py -r -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType ='python')

# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logger.info('Faraday rotation correction...')
MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DPPP")

#####################################################
# 3: find amplitude

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/amp.h5')
MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS sol.parmdb=$pathMS/amp.h5 sol.mode=fulljones', log='$nameMS_solG.log', commandType="DPPP")

lib_util.run_losoto(s, 'amp', [ms+'/amp.h5' for ms in MSs.getListStr()], \
        [parset_dir + '/losoto-flag.parset',parset_dir+'/losoto-plot-amp.parset',parset_dir+'/losoto-plot-ph.parset',parset_dir+'/losoto-bp.parset'])

sys.exit()

####################################################
# Re-do correcitons in right order

# Pol align correction DATA -> CORRECTED_DATA
logger.info('Polalign correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-pa.h5 cor.correction=polalign', log='$nameMS_corPA2.log', commandType="DPPP")
# Correct amp BP CORRECTED_DATA -> CORRECTED_DATA
logger.info('AmpBP correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-amp.h5 \
        cor.correction=amplitudeSmooth000 cor.updateweights=True', log='$nameMS_corAMP.log', commandType="DPPP")
# Beam correction CORRECTED_DATA -> CORRECTED_DATA
logger.info('Beam correction...')
MSs.run("DPPP " + parset_dir + '/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam2.log', commandType="DPPP")
# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logger.info('Faraday rotation correction...')
MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DPPP")
# Correct abs FR CORRECTED_DATA -> CORRECTED_DATA
MSs.run('createRMh5parm.py $pathMS $pathMS/amp.h5', log='$nameMS_RME.log', commandType="general")
MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.h5parm=$pathMS/amp.h5 cor.correction=rotationmeasure', log='$nameMS_corRME.log', commandType="DPPP")

#################################################
# 4: find iono

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth4.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/iono.h5')
MSs.run('DPPP '+parset_dir+'/DPPP-sol.parset msin=$pathMS sol.parmdb=$pathMS/iono.h5', log='$nameMS_solG2.log', commandType="DPPP")

# if field model available, subtract it
field_model = '/home/fdg/scripts/LiLF/models/calfields/'+calname+'-field.skydb'
if os.path.exists(field_model):
    logger.info('Removing field sources...')

    logger.info('Ft+corrupt model...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+field_model+' \
                pre.applycal.parmdb=$pathMS/iono.h5 pre.applycal.correction=phase000', log='$nameMS_field_pre.log', commandType="DPPP")

    # Remove the field sources CORRECTED_DATA -> CORRECTED_DATA - MODEL_DATA
    logger.info('Subtract model...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_field_taql1.log', commandType ='general')

    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_field_smooth.log', commandType ='python', maxThreads=20)

    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for MS in MSs.getListStr():
        lib_util.check_rm(ms+'/iono.h5')
    MSs.run('DPPP '+parset_dir+'/DPPP-sol.parset msin=$pathMS sol.parmdb=$pathMS/iono.h5', log='$nameMS_field_sol.log', commandType="DPPP")

if iono3rd:
    lib_util.run_losoto(s, 'iono', [ms+'/iono.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-flag.parset', parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-iono3rd.parset'])
else:
    lib_util.run_losoto(s, 'iono', [ms+'/iono.h5' for ms in MSs.getListStr()], \
            [parset_dir+'/losoto-flag.parset', parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-iono.parset'])

if 'survey' in os.getcwd():
    logger.info('Copy survey caltable...')
    cal = 'cal_'+os.getcwd().split('/')[-2]+'_'+calname
    logger.info('Copy: cal*h5 -> dsk:/disks/paradata/fdg/LBAsurvey/%s' % cal)
    os.system('ssh dsk "rm -rf /disks/paradata/fdg/LBAsurvey/%s"' % cal)
    os.system('ssh dsk "mkdir /disks/paradata/fdg/LBAsurvey/%s"' % cal)
    os.system('scp -q cal*h5 dsk:/disks/paradata/fdg/LBAsurvey/%s' % cal)

# a debug image
if imaging:
    logger.info("Imaging section:")

    if 'bootes2' in os.getcwd():
        fourth = int(len(glob.glob('./*MS'))/4.)
        MSs = lib_ms.AllMSs( sorted(glob.glob('./*MS'))[1*fourth:3*fourth], s ) # keep only mid band
    elif iono3rd:
        MSs = lib_ms.AllMSs( sorted(glob.glob('./*MS'))[int(len(glob.glob('./*MS'))/2.):], s ) # keep only upper band

    # Correct all CORRECTED_DATA (PA, beam, FR, BP corrected) -> CORRECTED_DATA
    logger.info('Amp/ph correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.steps=[ph,amp] cor.ph.parmdb=cal-iono.h5 cor.amp.parmdb=cal-iono.h5 \
        cor.ph.correction=phaseOrig000 cor.amp.correction=amplitude000 cor.amp.updateweights=False', log='$nameMS_corG.log', commandType="DPPP")

    logger.info('Subtract model...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql2.log', commandType ='general')

    logger.info('Flag...')
    MSs.run('DPPP '+parset_dir+'/DPPP-flag2.parset msin=$pathMS', log='$nameMS_flag2.log', commandType='DPPP')


    lib_util.check_rm('img')
    os.makedirs('img')

    if 'INNER' in obsmode: size = 8000
    elif 'SPARSE' in obsmode: size = 6000
    elif 'OUTER' in obsmode: size = 4000
    else: 
        logging.error('Observing mode not recognised.')
        sys.exit(1)

    logger.info('Cleaning low-res...')
    imagename = 'img/calLR'
    s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size/5) + ' ' + str(size/5) + ' -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 60arcsec -weight briggs 1.5 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
            -pol IUQV -join-channels -fit-spectral-pol 2 -channels-out 10 -apply-primary-beam -maxuv-l 1000 '+MSs.getStrWsclean(), \
            log='wscleanLR.log', commandType='wsclean', processors='max')
    s.run(check=True)

    logger.info('Cleaning normal...')
    imagename = 'img/cal'
    s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size) + ' ' + str(size) + ' -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.9 -minuv-l 100 \
            -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 20 '+MSs.getStrWsclean(), \
            log='wscleanA.log', commandType ='wsclean', processors='max')
    s.run(check = True)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits')
    im.makeMask(threshisl = 5)

    logger.info('Cleaning w/ mask...')
    imagename = 'img/calM'
    # -apply-primary-beam -use-differential-lofar-beam \
    s.add('wsclean -reorder -name ' + imagename + ' -size ' + str(size) + ' ' + str(size) + ' -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.8 -minuv-l 100 \
            -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 0.1 -save-source-list \
            -fits-mask '+im.maskname+' '+MSs.getStrWsclean(), \
            log='wscleanB.log', commandType='wsclean', processors = 'max')
    s.run(check = True)
    os.system('cat logs/wscleanB.log | grep "background noise"')

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits')
    im.makeMask(threshisl = 3)

    # apply mask
    import lsmtool
    logger.info('Predict (apply mask)...')
    lsm = lsmtool.load(imagename+'-sources.txt')
    lsm.select('%s == True' % (imagename+'-mask.fits'))
    cRA, cDEC = MSs.getListObj()[0].getPhaseCentre()
    lsm.select( lsm.getDistance(cRA, cDEC) > 0.1 ) # remove very centra part
    lsm.group('every')
    lsm.write(imagename+'-sources-cut.txt', format='makesourcedb', clobber = True)
    del lsm

logger.info("Done.")
