#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline to run on the calibrator observation.
# It isolates various systematic effects and
# prepare them for the transfer to the target field.

import sys, os, glob, re
import numpy as np

parset_dir = "/home/fdg/scripts/LiLF/parsets/LOFAR_cal"
imaging    = True

# Temporary!
if 'tooth' in os.getcwd(): # tooth 2013
    datadir = '../cals-bkp/'
    bl2flag = 'CS031LBA'
elif 'bootes' in os.getcwd(): # bootes 2013
    datadir = '../cals-bkp/'
    bl2flag = 'CS013LBA\;CS031LBA'
elif 'survey' in os.getcwd():
    obs     = os.getcwd().split('/')[-2] # assumes .../c??-o??/3c196
    calname = os.getcwd().split('/')[-1] # assumes .../c??-o??/3c196
    datadir = '../../download/%s/%s' % (obs, calname)
    bl2flag = 'CS031LBA\;RS310LBA\;RS210LBA\;RS409LBA'
    if 'c07-o00' in os.getcwd() or 'c07-o01' in os.getcwd() or 'c07-o02' in os.getcwd() or 'c07-o03' in os.getcwd() or 'c07-o04' in os.getcwd() or 'c07-o05' in os.getcwd() or 'c07-o06' in os.getcwd():
        bl2flag = 'CS031LBA\;RS310LBA\;RS210LBA\;RS409LBA\;RS407LBA'
    if 'c09' in os.getcwd(): bl2flag = 'CS031LBA\;CS013LBA'
else:
    datadir = '../cals-bkp/'
    bl2flag = ''

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
lib_log.set_logger('pipeline-cal.logger')
logger = lib_log.logger
lib_util.check_rm('logs')
s = lib_util.Scheduler(dry = False)
MSs = lib_ms.AllMSs( glob.glob(datadir+'/*MS'), s )

# copy data
logger.info('Copy data...')
for MS in MSs.getListObj():
    MS.move(MS.nameMS+'.MS', keepOrig=True)
MSs = lib_ms.AllMSs( glob.glob('*MS'), s )
calname = MSs.getListObj()[0].getNameField()

## flag bad stations, flags will propagate
#logger.info("Flagging...")
#MSs.run("DPPP " + parset_dir + "/DPPP-flag.parset msin=$pathMS flag1.baseline=" + bl2flag, log="$nameMS_flag.log", commandType="DPPP")
#
## predict to save time ms:MODEL_DATA
#logger.info('Add model to MODEL_DATA...')
#skymodel = "/home/fdg/scripts/LiLF/models/calib-simple.skydb"
#MSs.run("DPPP " + parset_dir + "/DPPP-predict.parset msin=$pathMS pre.sourcedb=" + skymodel + " pre.sources=" + calname, log="$nameMS_pre.log", commandType="DPPP")
##MSs.run("DPPP " + parset_dir + "/DPPP-predict.parset msin=$pathMS pre.sourcedb=/home/fdg/scripts/model/3C196-allfield.skydb", log="$nameMS_pre.log", commandType="DPPP")

##################################################
# 1: find the FR and remove it

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/fr-lin.h5')
MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS sol.caltype=diagonal sol.parmdb=$pathMS/fr-lin.h5', log='$nameMS_sol0a.log', commandType="DPPP")

lib_util.run_losoto(s, 'fr-lin', [ms+'/fr-lin.h5' for ms in MSs.getListStr()], [parset_dir+'/losoto-align.parset'])

logger.info('Polalign correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-fr-lin.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DPPP")

# Beam correction CORRECTED_DATA -> CORRECTED_DATA
logger.info('Beam correction...')
MSs.run("DPPP " + parset_dir + '/DPPP-beam.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA corrbeam.updateweights=True', log='$nameMS_beam1.log', commandType="DPPP")

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/fr-lin2.h5')
MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS sol.caltype=fulljones sol.parmdb=$pathMS/fr-lin2.h5', log='$nameMS_sol0b.log', commandType="DPPP")

lib_util.run_losoto(s, 'fr-lin2', [MS+'/fr-lin2.h5' for MS in MSs.getListStr()], [parset_dir + '/losoto-amp.parset', parset_dir+'/losoto-align.parset'])
sys.exit()

# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
logger.info('Converting to circular...')
MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType ='python')

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/fr.h5')
MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS sol.parmdb=$pathMS/fr.h5', log='$nameMS_sol1.log', commandType="DPPP")

lib_util.run_losoto(s, 'fr', [ms+'/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset', parset_dir + '/losoto-amp.parset'])

#####################################################
# 2: find amplitude + align

# Beam correction DATA -> CORRECTED_DATA
logger.info('Beam correction...')
MSs.run('DPPP ' + parset_dir + '/DPPP-beam.parset msin=$pathMS', log='$nameMS_beam2.log', commandType="DPPP")

# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logger.info('Faraday rotation correction...')
MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DPPP")

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/amp.h5')
MSs.run('DPPP ' + parset_dir + '/DPPP-sol.parset msin=$pathMS sol.parmdb=$pathMS/amp.h5', log='$nameMS_sol2.log', commandType="DPPP")

lib_util.run_losoto(s, 'amp', [ms+'/amp.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-align.parset'])

##################################################
# 3: find iono

# Beam correction (and update weight in case of imaging) DATA -> CORRECTED_DATA
logger.info('Beam correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS msin.datacolumn=DATA', log='$nameMS_beam3.log', commandType="DPPP")

# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logger.info('Faraday rotation correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR2.log', commandType="DPPP")

# Correct amp CORRECTED_DATA -> CORRECTED_DATA
logger.info('AmpBP correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-amp.h5 \
        cor.correction=amplitudeSmooth000 cor.updateweights=True', log='$nameMS_corAMP.log', commandType="DPPP")

# Correct PA CORRECTED_DATA -> CORRECTED_DATA
# could be done in the next losoto call but is a good debug to check if it is not present at calibration time
logger.info('Polalign correction...')
MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-amp.h5 cor.correction=polalign', log='$nameMS_corPA.log', commandType="DPPP")

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logger.info('BL-smooth...')
MSs.run('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3.log', commandType ='python', maxThreads=20)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logger.info('Calibrating...')
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/iono.h5')
MSs.run('DPPP '+parset_dir+'/DPPP-sol.parset msin=$pathMS sol.parmdb=$pathMS/iono.h5', log='$nameMS_sol3.log', commandType="DPPP")

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

lib_util.run_losoto(s, 'iono', [ms+'/iono.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-flag.parset', parset_dir + '/losoto-iono.parset'])

# TODO: bisogna tenere conto che possono essrci piu' calibratori in una run!
if 'survey' in os.getcwd():
    logger.info('Copy survey caltable...')
    cal = 'cal_'+os.getcwd().split('/')[-2]
    logger.info('Copy: cal*h5 -> dsk:/disks/paradata/fdg/LBAsurvey/%s' % cal)
    os.system('ssh dsk "rm -rf /disks/paradata/fdg/LBAsurvey/%s"' % cal)
    os.system('ssh dsk "mkdir /disks/paradata/fdg/LBAsurvey/%s"' % cal)
    os.system('scp -q cal*h5 dsk:/disks/paradata/fdg/LBAsurvey/%s' % cal)

# a debug image
if imaging:
    logger.info("Imaging seciotn:")

    if not 'survey' in os.getcwd():
        MSs = lib_ms.AllMSs( sorted(glob.glob('./*MS'))[int(len(glob.glob('./*MS'))/2.):], s ) # keep only upper band

    # Correct all CORRECTED_DATA (beam, CD, FR, BP corrected) -> CORRECTED_DATA
    logger.info('Amp/ph correction...')
    MSs.run("DPPP " + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-iono.h5 cor.steps=[ph,amp] \
        cor.ph.correction=phaseOrig000 cor.amp.correction=amplitude000 cor.amp.updateweights=False', log='$nameMS_corG.log', commandType="DPPP")

    logger.info('Subtract model...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql2.log', commandType ='general')

    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')
    imagename = 'img/wide'
    s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.9 -minuv-l 100 \
            -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 20 '+MSs.getStrWsclean(), \
            log='wscleanA.log', commandType ='wsclean', processors='max')
    s.run(check = True)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits')
    im.makeMask(threshisl = 5)

    logger.info('Cleaning w/ mask...')
    imagename = 'img/wideM'
    # -apply-primary-beam -use-differential-lofar-beam
    s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.8 -minuv-l 100 \
            -pol I -join-channels -fit-spectral-pol 2 -channels-out 10 -auto-threshold 0.1 -save-source-list \
            -fits-mask '+im.maskname+' '+MSs.getStrWsclean(), \
            log='wscleanB.log', commandType='wsclean', processors = 'max')
    s.run(check = True)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits')
    im.makeMask(threshisl = 3)

    # apply mask
    import lsmtool
    logger.info('Predict (apply mask)...')
    lsm = lsmtool.load(imagename+'-sources.txt')
    lsm.select('%s == True' % (imagename+'-mask.fits'))
    cRA, cDEC = MSs.getListObj[0].getPhaseCentre()
    lsm.select( lsm.getDistance(cRA, cDEC) > 0.1 ) # remove very centra part
    lsm.group('every')
    lsm.write(imagename+'-sources-cut.txt', format='makesourcedb', clobber = True)
    del lsm

logger.info("Done.")
