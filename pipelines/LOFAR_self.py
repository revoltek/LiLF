#!/usr/bin/env python
# -*- coding: utf-8 -*-
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal corrected source subtracted data in CORRECTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions
# last high/low resolution models are copied in the "self/models" dir
# last high/low resolution images + masks + empty images (CORRECTED_DATA) are copied in the "self/images" dir
# h5parm solutions and plots are copied in the "self/solutions" dir

import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *
import pyrap.tables as pt
from make_mask import make_mask
import lsmtool

parset_dir = '/home/fdg/scripts/autocal/parset_self/'
skymodel = '/home/fdg/scripts/model/calib-simple.skymodel'
niter = 3
user_mask = None
cc_predict = True

if 'tooth' in os.getcwd():
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skydb'
    apparent = True # no beam correction
    user_mask = '/home/fdg/scripts/autocal/regions/tooth.reg'
    multiepoch = False
elif 'bootes' in os.getcwd():
    sourcedb = '/home/fdg/scripts/model/Bootes_HBA.corr.skydb'
    apparent = False
    multiepoch = False
else:
    # Survey
    multiepoch = True
    obs = os.getcwd().split('/')[-1]
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/skymodels/%s.skydb' % obs
    apparent = False
    if not os.path.exists('mss'):
        os.makedirs('mss')
        for i, tc in enumerate(glob.glob('../../c*-o*/%s/mss/*' % obs)):
            tc_ren = 'TC%02i.MS' % i
            print 'cp -r %s mss/%s' % (tc,tc_ren)
            os.system('cp -r %s mss/%s' % (tc,tc_ren))

assert os.path.exists(sourcedb)

#############################################################################

def ft_model_wsclean(mss, imagename, c, user_mask = None, keep_in_beam=True, resamp = None, model_column='MODEL_DATA'):
    """
    mss : vector of mss
    imagename : root name for wsclean model images
    resamp : must be '10asec' or another pixels size to resample models
    keep_in_beam : if True remove everything outside primary beam, otherwise everything inside
    """
    logger.info('Predict with model image...')

    # remove CC not in mask
    logger.info('Predict (mask)...')
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 5, atrous_do=True)
    if user_mask is not None: 
        blank_image_reg(maskname, user_mask, inverse=False, blankval=1)
    blank_image_reg(maskname, 'self/beam.reg', inverse=keep_in_beam)
    for modelname in sorted(glob.glob(imagename+'*model.fits')):
        blank_image_fits(modelname, maskname, inverse=True)

    if resamp is not None:
        logger.info('Predict (resamp)...')
        for model in sorted(glob.glob(imagename+'*model.fits')):
            model_out = model.replace(imagename, imagename+'-resamp')
            s.add('/home/fdg/opt/src/nnradd/build/nnradd '+resamp+' '+model_out+' '+model, log='resamp-c'+str(c)+'.log', log_append=True, cmd_type='general')
        s.run(check=True)
        imagename = imagename+'-resamp'
 
    logger.info('Predict (ft)...')
    s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+' '.join(mss), \
            log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    if model_column != 'MODEL_DATA':
        logger.info('Predict (set %s = MODEL_DATA)...' % model_column)
        for ms in mss:
            s.add('taql "update '+ms+' set '+model_column+' = MODEL_DATA"', log=ms+'_taql0-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)


def ft_model_cc(mss, imagename, c, user_mask = None, keep_in_beam=True, model_column='MODEL_DATA'):
    """
    skymodel : cc-list made by wsclean
    keep_in_beam : if True remove everything outside primary beam, otherwise everything inside
    """
    logger.info('Predict with CC...')
    maskname = imagename+'-mask.fits'
    skymodel = imagename+'-sources.txt'
    skymodel_cut = imagename+'-sources-cut.txt'
    skydb = imagename+'-sources.skydb'

    # prepare mask
    if not os.path.exists(maskname):
        logger.info('Predict (make mask)...')
        make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 5, atrous_do=True)
    if user_mask is not None:
        blank_image_reg(maskname, user_mask, inverse=False, blankval=1) # set to 1 pixels into user_mask
    blank_image_reg(maskname, 'self/beam.reg', inverse=keep_in_beam, blankval=0) # if keep_in_beam set to 0 everything outside beam.reg

    # apply mask
    logger.info('Predict (apply mask)...')
    lsm = lsmtool.load(skymodel)
    lsm.select('%s == True' % maskname)
    fluxes = lsm.getColValues('I')
    #lsm.remove(np.abs(fluxes) < 5e-4) # TEST
    lsm.write(skymodel_cut, format='makesourcedb', clobber=True)
    del lsm

    # convert to skydb
    logger.info('Predict (makesourcedb)...')
    check_rm(skydb)
    s.add('makesourcedb outtype="blob" format="<" in="'+skymodel_cut+'" out="'+skydb+'"', log='makesourcedb-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    # predict
    logger.info('Predict (ft)...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' msout.datacolumn='+model_column+' pre.usebeammodel=false pre.sourcedb='+skydb, \
                log=ms+'_pre-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)


#############################################################################

logger = set_logger('pipeline-self.logger')
check_rm('logs')
s = Scheduler(dry=False)

##################################################
# Clear
logger.info('Cleaning...')

check_rm('img')
os.makedirs('img')
os.makedirs('logs/mss')

# here images, models, solutions for each group will be saved
check_rm('self')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

mss = sorted(glob.glob('mss/TC*[0-9].MS'))
concat_ms = 'mss/concat.MS'

# make beam
phasecentre = get_phase_centre(mss[0])
make_beam_reg(phasecentre[0], phasecentre[1], 12, 'self/beam.reg') # go to 12 deg, first null
#make_beam_reg(phasecentre[0], phasecentre[1], 8, 'self/beam.reg') # go to 7 deg, first null

###############################################################################################
# Create columns (non compressed)
logger.info('Creating MODEL_DATA_HIGHRES and SUBTRACTED_DATA...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
s.run(check=True)

##################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for ms in mss:
    check_rm(ms+'/'+sourcedb_basename)
    logger.debug('Copy: '+sourcedb+' -> '+ms)
    os.system('cp -r '+sourcedb+' '+ms)
logger.info('Add model to MODEL_DATA...')
for ms in mss:
    if apparent:
        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.usebeammodel=false pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP')
    else:
        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.usebeammodel=true pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP')
s.run(check=True)

###################################################################################
# Preapre fake FR parmdb
logger.info('Prepare fake FR parmdb...')
for ms in mss:
    if os.path.exists(ms+'/instrument-fr'): continue
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=replace(NAME,\':@MODEL_DATA\',\'\')"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

#####################################################################################################
# Self-cal cycle
for c in xrange(niter):

    logger.info('Start selfcal cycle: '+str(c))

    # Smooth DATA -> SMOOTHED_DATA
    # Re-done in case of new flags
    if c == 0:
        incol = 'DATA'
    else:
        incol = 'SUBTRACTED_DATA'

    logger.info('BL-based smoothing...')
    for ms in mss:
        s.add('BLsmooth.py -r -f 0.2 -i '+incol+' -o SMOOTHED_DATA '+ms, log=ms+'_smooth1-c'+str(c)+'.log', cmd_type='python')
    s.run(check=True, max_threads=6)

    logger.info('Concatenating TCs...')
    check_rm(concat_ms+'*')
    pt.msutil.msconcat(mss, concat_ms, concatTime=False)

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving TEC...')
    for ms in mss:
        check_rm(ms+'/instrument-tec')
        s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' sol.parmdb='+ms+'/instrument-tec', \
                log=ms+'_solTEC-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # LoSoTo plot
    if multiepoch:
        for i, ms in enumerate(mss):
            run_losoto(s, 'tec'+str(c)+'-ms'+str(i), [ms], [parset_dir+'/losoto-plot.parset'], ininstrument='instrument-tec', putback=False)
    else:
        run_losoto(s, 'tec'+str(c), mss, [parset_dir+'/losoto-plot.parset'], ininstrument='instrument-tec', putback=False)
    os.system('mv plots-tec'+str(c)+'* self/solutions/')
    os.system('mv cal-tec'+str(c)+'*.h5 self/solutions/')

    # correct TEC - group*_TC.MS:(SUBTRACTED_)DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correcting TEC...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn='+incol+' cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec', \
                log=ms+'_corTEC-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    #####################################################################################################
    # Cross-delay + Faraday rotation correction
    if c >= 1:

        # To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (circular)
        # TODO: check -w, is it ok?
        logger.info('Convert to circular...')
        for ms in mss:
            s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=4)
 
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=6)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Solving G...')
        for ms in mss:
            check_rm(ms+'/instrument-g')
            s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' sol.parmdb='+ms+'/instrument-g sol.solint=30 sol.nchan=8', \
                    log=ms+'_sol-g1-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        if multiepoch:
            for i, ms in enumerate(mss):
                run_losoto(s, 'fr'+str(c)+'-ms'+str(i), [ms], [parset_dir+'/losoto-fr.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
                outinstrument='instrument-fr', outglobaldb='globaldb-fr', outtab='rotationmeasure000', putback=True)
        else:
            run_losoto(s, 'fr'+str(c), mss, [parset_dir+'/losoto-fr.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
            outinstrument='instrument-fr', outglobaldb='globaldb-fr', outtab='rotationmeasure000', putback=True)
        os.system('mv plots-fr'+str(c)+'* self/solutions/')
        os.system('mv cal-fr'+str(c)+'*.h5 self/solutions/')
       
        # To linear - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (linear)
        logger.info('Convert to linear...')
        for ms in mss:
            s.add('/home/fdg/scripts/mslin2circ.py -r -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=4)
        
        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth3-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=6)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Solving G...')
        for ms in mss:
            check_rm(ms+'/instrument-g')
            s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' sol.parmdb='+ms+'/instrument-g sol.solint=30 sol.nchan=8', \
                    log=ms+'_sol-g2-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        if multiepoch:
            for i, ms in enumerate(mss):
                run_losoto(s, 'cd'+str(c)+'-ms'+str(i), [ms], [parset_dir+'/losoto-cd.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
                    outinstrument='instrument-cd', outglobaldb='globaldb', outtab='amplitude000,crossdelay', putback=True)
        else:
            run_losoto(s, 'cd'+str(c), mss, [parset_dir+'/losoto-cd.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
                outinstrument='instrument-cd', outglobaldb='globaldb', outtab='amplitude000,crossdelay', putback=True)
        os.system('mv plots-cd'+str(c)+'* self/solutions/')
        os.system('mv cal-cd'+(str(c))+'*.h5 self/solutions/')

        if multiepoch:
            for i, ms in enumerate(mss):
                run_losoto(s, 'amp'+str(c)+'-ms'+str(i), [ms], [parset_dir+'/losoto-amp.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
                    outinstrument='instrument-amp', outglobaldb='globaldb', outtab='amplitude000,phase000', putback=True)
        else:
            run_losoto(s, 'amp'+str(c), mss, [parset_dir+'/losoto-amp.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
                outinstrument='instrument-amp', outglobaldb='globaldb', outtab='amplitude000,phase000', putback=True)
        os.system('mv plots-amp'+str(c)+'* self/solutions/')
        os.system('mv cal-amp'+str(c)+'*.h5 self/solutions/')

        # Correct CD SB.MS:SUBTRACTED_DATA->CORRECTED_DATA
        logger.info('Cross-delay correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=SUBTRACTED_DATA cor.parmdb='+ms+'/instrument-cd cor.correction=Gain', log=ms+'_corCD-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
        # Correct beam amp SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logger.info('Beam amp correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cor.parmdb='+ms+'/instrument-amp cor.correction=Gain', log=ms+'_corAMP-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', \
                    log=ms+'_corFR-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # Finally re-calculate TEC
        logger.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLsmooth.py -r -f 0.2 -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth3-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=6)

        # solve TEC - group*_TC.MS:SMOOTHED_DATA
        logger.info('Solving TEC...')
        for ms in mss:
            check_rm(ms+'/instrument-tec')
            s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' sol.parmdb='+ms+'/instrument-tec', \
                    log=ms+'_solTEC-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # LoSoTo plot
        if multiepoch:
            for i, ms in enumerate(mss):
                run_losoto(s, 'tec'+str(c)+'b-ms'+str(i), [ms], [parset_dir+'/losoto-plot.parset'], ininstrument='instrument-tec', putback=False)
        else:
            run_losoto(s, 'tec'+str(c)+'b', mss, [parset_dir+'/losoto-plot.parset'], ininstrument='instrument-tec', putback=False)
        os.system('mv plots-tec'+str(c)+'b* self/solutions')
        os.system('mv cal-tec'+str(c)+'b*.h5 self/solutions')

        # correct TEC - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting TEC...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec', \
                    log=ms+'_corTECb-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

    # do beam-corrected+deeper image at last cycle
    if c == niter-1:
        # beam corrected: -use-differential-lofar-beam' - no baseline avg!
        logger.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -trim 3500 3500 -mem 90 -j '+str(s.max_processors)+' \
                -scale 8arcsec -weight briggs 0.0 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -apply-primary-beam -use-differential-lofar-beam -minuv-l 30 '+' '.join(mss), \
                log='wscleanBeam-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

        logger.info('Cleaning beam high-res (cycle: '+str(c)+')...')
        imagename = 'img/wideBeamHR'
        s.add('wsclean -reorder -name ' + imagename + ' -size 6000 6000 -trim 5500 5500 -mem 90 -j '+str(s.max_processors)+' \
                -scale 4arcsec -weight briggs -1.5 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -apply-primary-beam -use-differential-lofar-beam -minuv-l 30 '+' '.join(mss), \
                log='wscleanBeamHR-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

    # clean mask clean (cut at 5k lambda)
    # no MODEL_DATA update with -baseline-averaging
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -trim 3500 3500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.9 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 20 -minuv-l 30 '+' '.join(mss), \
            log='wsclean-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 3, atrous_do=True)
    if user_mask is not None: 
        blank_image_reg(maskname, user_mask, inverse=False, blankval=1)

    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    #TODO: -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,9 \
    if cc_predict:
        #s.add('wsclean -reorder -name ' + imagename + ' -size 3000 3000 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
        s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -trim 3500 3500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -niter 1000000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 0.1 -minuv-l 30 -save-source-list -fitsmask '+maskname+' '+' '.join(mss), \
            log='wscleanM-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    else:
        #s.add('wsclean -reorder -name ' + imagename + ' -size 3000 3000 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
        s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -trim 3500 3500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -niter 1000000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 0.1 -minuv-l 30 -fitsmask '+maskname+' '+' '.join(mss), \
            log='wscleanM-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-c'+str(c)+'.log | grep "background noise"')

    if c > 0 and c != niter:
        if cc_predict:
            ft_model_cc(mss, imagename, c, user_mask=user_mask, keep_in_beam=True, model_column='MODEL_DATA')
        else:
            ft_model_wsclean(mss, imagename, c, user_mask=user_mask, resamp='10asec', keep_in_beam=True)
    # do low-res first cycle and remove it from the data
    if c == 0:
        if cc_predict:
            ft_model_cc(mss, imagename, c, user_mask=user_mask, keep_in_beam=True, model_column='MODEL_DATA_HIGHRES')
        else:
            ft_model_wsclean(mss, imagename, c, user_mask=user_mask, resamp='10asec', keep_in_beam=True, model_column='MODEL_DATA_HIGHRES')

        # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_HIGHRES)...')
        s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_HIGHRES"', log='taql1-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        # reclean low-resolution
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        if cc_predict:
            #s.add('wsclean -reorder -name ' + imagename_lr + ' -size 4500 4500 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            s.add('wsclean -reorder -name ' + imagename_lr + ' -size 6000 6000 -trim 5500 5500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
                -scale 20arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 2000 -mgain 0.8 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 1 -minuv-l 100 -save-source-list '+' '.join(mss), \
                log='wsclean-lr.log', cmd_type='wsclean', processors='max')
        else:
            #s.add('wsclean -reorder -name ' + imagename_lr + ' -size 4500 4500 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            s.add('wsclean -reorder -name ' + imagename_lr + ' -size 6000 6000 -trim 5500 5500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
                -scale 20arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 2000 -mgain 0.8 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 1 -minuv-l 100 '+' '.join(mss), \
                log='wsclean-lr.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
       
        if cc_predict:
            ft_model_cc(mss, imagename_lr, 'lr', keep_in_beam=False, model_column='MODEL_DATA')
        else:
            ft_model_wsclean(mss, imagename_lr, 'lr', user_mask=None, resamp='10asec', keep_in_beam=False)

        # corrupt model with TEC solutions ms:MODEL_DATA -> ms:MODEL_DATA
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                cor1.parmdb='+ms+'/instrument-tec cor1.invert=false cor2.parmdb='+ms+'/instrument-tec cor2.invert=false', \
                log=ms+'_corrupt.log', cmd_type='NDPPP')
        s.run(check=True)
    
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
        logger.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='taql2-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)

        # Restore best model
        logger.info('Restoring high-res model (MODEL_DATA = MODEL_DATA_HIGHRES)...')
        s.add('taql "update '+concat_ms+' set MODEL_DATA = MODEL_DATA_HIGHRES"', log='taql3-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)


    ###############################################################################################################
    # Flag on residuals (CORRECTED_DATA)
    #logger.info('Flagging residuals...')
    #for ms in mss:
    #    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
    #s.run(check=True
    
# make beam
# TODO: remove when wsclean will produce a proper primary beam
os.system('~/opt/src/makeavgpb/build/wsbeam.py img/wideBeam')

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
if cc_predict: [ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in xrange(niter) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wideBeam-MFS-image.fits  img/wideBeam-MFS-image-pb.fits self/images')
os.system('mv img/wideBeamHR-MFS-image.fits  img/wideBeamHR-MFS-image-pb.fits self/images')
os.system('mv logs self')

logger.info("Done.")
