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
import pyrap.tables as pt
import lsmtool

parset_dir = '/home/fdg/scripts/autocal/parset_self'
niter = 3

# Temporary
if 'tooth' in os.getcwd():
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skydb'
    apparent = True # no beam correction
    user_mask = '/home/fdg/scripts/autocal/regions/tooth.reg'
    multiepoch = False
elif 'bootes' in os.getcwd():
    sourcedb = '/home/fdg/scripts/model/Bootes_HBA.corr.skydb'
    apparent = False
    user_mask = None
    multiepoch = False
else:
    # Survey
    obs = os.getcwd().split('/')[-1]
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/skymodels/%s.skydb' % obs
    apparent = False
    user_mask = None
    multiepoch = True
    if not os.path.exists('mss'):
        os.makedirs('mss')
        for i, tc in enumerate(glob.glob('../../c*-o*/%s/mss/*' % obs)):
            tc_ren = 'TC%02i.MS' % i
            print 'cp -r %s mss/%s' % (tc,tc_ren)
            os.system('cp -r %s mss/%s' % (tc,tc_ren))

assert os.path.exists(sourcedb)

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, make_mask
lib_log.set_logger('pipeline-self.logger')
logger = lib_log.logger
lib_util.check_rm('logs')
s = lib_util.Scheduler(dry = False)

#############################################################################

#def ft_model_wsclean(mss, imagename, c, user_mask = None, keep_in_beam=True, resamp = None, model_column='MODEL_DATA'):
#    """
#    mss : vector of mss
#    imagename : root name for wsclean model images
#    resamp : must be '10asec' or another pixels size to resample models
#    keep_in_beam : if True remove everything outside primary beam, otherwise everything inside
#    """
#    logger.info('Predict with model image...')
#
#    # remove CC not in mask
#    logger.info('Predict (mask)...')
#    maskname = imagename+'-mask.fits'
#    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 5, atrous_do=True)
#    if user_mask is not None: 
#        blank_image_reg(maskname, user_mask, inverse=False, blankval=1)
#    blank_image_reg(maskname, 'self/beam.reg', inverse=keep_in_beam)
#    for modelname in sorted(glob.glob(imagename+'*model.fits')):
#        blank_image_fits(modelname, maskname, inverse=True)
#
#    if resamp is not None:
#        logger.info('Predict (resamp)...')
#        for model in sorted(glob.glob(imagename+'*model.fits')):
#            model_out = model.replace(imagename, imagename+'-resamp')
#            s.add('/home/fdg/opt/src/nnradd/build/nnradd '+resamp+' '+model_out+' '+model, log='resamp-c'+str(c)+'.log', log_append=True, cmdType='general')
#        s.run(check=True)
#        imagename = imagename+'-resamp'
# 
#    logger.info('Predict (ft)...')
#    s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+' '.join(mss), \
#            log='wscleanPRE-c'+str(c)+'.log', cmdType='wsclean', processors='max')
#    s.run(check=True)
#
#    if model_column != 'MODEL_DATA':
#        logger.info('Predict (set %s = MODEL_DATA)...' % model_column)
#        for ms in mss:
#            s.add('taql "update '+ms+' set '+model_column+' = MODEL_DATA"', log=ms+'_taql0-c'+str(c)+'.log', cmdType='general')
#        s.run(check=True)


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
    lib_util.check_rm(skydb)
    s.add('makesourcedb outtype="blob" format="<" in="'+skymodel_cut+'" out="'+skydb+'"', log='makesourcedb-c'+str(c)+'.log', cmdType='general')
    s.run(check=True)

    # predict
    logger.info('Predict (ft)...')
    for ms in mss:
        s.add('DPPP '+parset_dir+'/DPPP-predict.parset msin='+ms+' msout.datacolumn='+model_column+' pre.usebeammodel=false pre.sourcedb='+skydb, \
                log=ms+'_pre-c'+str(c)+'.log', cmdType='DPPP')
    s.run(check=True)

#############################################################################
# Clear
logger.info('Cleaning...')

lib_util.check_rm('img')
os.makedirs('img')
os.makedirs('logs/mss')

# here images, models, solutions for each group will be saved
lib_util.check_rm('self')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

MSs = lib_ms.MSs( glob.glob('mss/TC*[0-9].MS') )
concat_ms = 'mss/concat.MS'

# make beam
phasecentre = MSs.getListObs()[0].getPhaseCentre()
MSs.getListObs()[0].makeBeamReg('self/beam.reg') # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null

###############################################################################################
# Create columns (non compressed)
logger.info('Creating MODEL_DATA_HIGHRES and SUBTRACTED_DATA...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log='$nameMS_addcol.log', cmdType='python')

##################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr:
    lib_util.check_rm(MS+'/'+sourcedb_basename)
    logger.debug('Copy: '+sourcedb+' -> '+MS)
    os.system('cp -r '+sourcedb+' '+MS)

logger.info('Add model to MODEL_DATA...')
if apparent:
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.usebeammodel=false pre.sourcedb=$pathMS/'+sourcedb_basename, log'$nameMS_pre.log', cmdType='DPPP')
else:
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.usebeammodel=true pre.sourcedb=$pathMS/'+sourcedb_basename, log='$nameMS_pre.log', cmdType='DPPP')

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
    MSs.run('BLsmooth.py -r -f 0.2 -i '+incol+' -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1-c'+str(c)+'.log', cmdType='python', maxThread=6)

    logger.info('Concatenating TCs...')
    lib_util.check_rm(concat_ms+'*')
    pt.msutil.msconcat(MSs.getListStr(), concat_ms, concatTime=False)

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving TEC...')
    for MS in MSs.get_list_str():
        lib_util.check_rm(MS+'/tec.h5')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS sol.parmdb=$pathMS/tec.h5', \
                log='$nameMS_solTEC-c'+str(c)+'.log', cmdType='DPPP')

    # LoSoTo plot
    if multiepoch:
        for i, MS in enumerate(MSs.getListStr()):
            lib_util.run_losoto(s, 'tec'+str(c)+'-ms'+str(i), [MS+'/tec.h5'], [parset_dir+'/losoto-plot.parset'])
    else:
        lib_util.run_losoto(s, 'tec'+str(c), [ms+'/tec.h5' for ms in MSs.get_list_str()], [parset_dir+'/losoto-plot.parset'], concat='time')
    os.system('mv plots-tec'+str(c)+'* self/solutions/')
    os.system('mv cal-tec'+str(c)+'*.h5 self/solutions/')

    # correct TEC - group*_TC.MS:(SUBTRACTED_)DATA -> group*_TC.MS:CORRECTED_DATA
    logger.info('Correcting TEC...')
    MSs.add('DPPP '+parset_dir+'/DPPP-corTEC.parset msin=$pathMS msin.datacolumn='+incol+' cor1.parmdb=$pathMS/tec.h5 cor2.parmdb=$pathMS/tec.h5', \
                log='$nameMS_corTEC-c'+str(c)+'.log', cmdType='DPPP')

    #####################################################################################################
    # Cross-delay + Faraday rotation correction
    if c >= 1:

        # To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (circular)
        # TODO: check -w, is it ok?
        logger.info('Convert to circular...')
        MSs.run('/home/fdg/scripts/mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin-c'+str(c)+'.log', cmdType='python', maxThread=4)
 
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth2-c'+str(c)+'.log', cmdType='python', maxThread=6)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Solving G...')
        for MS in MSs.get_list_str():
            lib_util.check_rm(MS+'/fr.h5')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.parmdb=$pathMS/fr.h5 sol.solint=30 sol.nchan=8', \
                    log='$nameMS_sol-g1-c'+str(c)+'.log', cmdType='DPPP')

        if multiepoch:
            for i, MS in enumerate(MSs.get_list_str()):
                lib_util.run_losoto(s, 'fr'+str(c)+'-ms'+str(i), [MS+'/fr.h5'], [parset_dir+'/losoto-fr.parset'])
        else:
            lib_util.run_losoto(s, 'fr'+str(c), [ms+'/fr.h5' for ms in MSs.get_list_str()], [parset_dir+'/losoto-fr.parset'])
        os.system('mv plots-fr'+str(c)+'* self/solutions/')
        os.system('mv cal-fr'+str(c)+'*.h5 self/solutions/')
       
        # To linear - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (linear)
        logger.info('Convert to linear...')
        MSs.run('/home/fdg/scripts/mslin2circ.py -r -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin-c'+str(c)+'.log', cmdType='python', maxThreads=4)
        
        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        if multiepoch: h5 = '$pathMS/fr.h5'
        else: h5 = 'cal-fr'+str(c)+'.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb='+h5+' cor.correction=rotationmeasure000', \
                    log='$nameMs_corFR-c'+str(c)+'.log', cmdType='DPPP')

        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -r -f 0.5 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3-c'+str(c)+'.log', cmdType='python', maxThreads=6)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logger.info('Solving G...')
        for MS in MSs.get_list_str():
            lib_util.check_rm(ms+'/amp.h5')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.parmdb=$pathMS/amp.h5 sol.solint=30 sol.nchan=8', \
                    log='$nameMS_sol-g2-c'+str(c)+'.log', cmdType='DPPP')

        if multiepoch:
            for i, MS in enumerate(MSs.get_list_str()):
                lib_util.run_losoto(s, 'amp'+str(c)+'-ms'+str(i), [MS+'/amp.h5'], [parset_dir+'/losoto-align.parset',parset_dir+'/losoto-amp.parset'])
        else:
            lib_util.run_losoto(s, 'amp'+str(c), [ms+'/amp.h5' for ms in MSs.get_list_str()], [parset_dir+'/losoto-align.parset'])
        os.system('mv plots-amp'+str(c)+'* self/solutions/')
        os.system('mv cal-amp'+(str(c))+'*.h5 self/solutions/')

        # Correct CD SB.MS:SUBTRACTED_DATA->CORRECTED_DATA
        logger.info('Pol-align correction...')
        if multiepoch: h5 = '$pathMS/amp.h5'
        else: h5 = 'cal-amp'+str(c)+'.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA cor.parmdb='+h5+' cor.correction=polalign', log=ms+'_corPA-c'+str(c)+'.log', cmdType='DPPP')
        # Correct beam amp SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logger.info('Beam amp correction...')
        if multiepoch: h5 = '$pathMS/amp.h5'
        else: h5 = 'cal-amp'+str(c)+'.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmd='+h5' cor.correction=amplitude000', log=ms+'_corAMP-c'+str(c)+'.log', cmdType='DPPP')
        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logger.info('Faraday rotation correction...')
        if multiepoch: h5 = '$pathMS/fr.h5'
        else: h5 = 'cal-fr'+str(c)+'.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb='+h5+' cor.correction=rotationmeasure000', \
                    log=ms+'_corFR-c'+str(c)+'.log', cmdType='DPPP')

        # Finally re-calculate TEC
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -r -f 0.2 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth3-c'+str(c)+'.log', cmdType='python', maxThreads=6)

        # solve TEC - group*_TC.MS:SMOOTHED_DATA
        logger.info('Solving TEC...')
        for MS in MSs.get_list_str():
            lib_util.check_rm(ms+'/tec.h5')
        MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS sol.parmdb=$pathMS/tec.h5', \
                    log=ms+'_solTEC-c'+str(c)+'.log', cmdType='DPPP')

        # LoSoTo plot
        if multiepoch:
            for i, MS in enumerate(MSs.get_list_str()):
                lib_util.run_losoto(s, 'tec'+str(c)+'b-ms'+str(i), [MS+'/tec.h5'], [parset_dir+'/losoto-plot.parset'])
        else:
            lib_util.run_losoto(s, 'tec'+str(c)+'b', [ms+'/tec.h5' for ms in MSs.get_list_str()], [parset_dir+'/losoto-plot.parset'])
        os.system('mv plots-tec'+str(c)+'b* self/solutions')
        os.system('mv cal-tec'+str(c)+'b*.h5 self/solutions')

        # correct TEC - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting TEC...')
        if multiepoch: h5 = '$pathMS/tec.h5'
        else: h5 = 'cal-tec'+str(c)+'b.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-corTEC.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor1.parmdb='+h5+' cor2.parmdb='+h5, \
                    log=ms+'_corTECb-c'+str(c)+'.log', cmdType='DPPP')

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

    # do beam-corrected+deeper image at last cycle
    if c == niter-1:
        # beam corrected: -use-differential-lofar-beam' - no baseline avg!
        logger.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -mem 90 -j '+str(s.max_processors)+' \
                -scale 8arcsec -weight briggs 0.0 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -apply-primary-beam -use-differential-lofar-beam -minuv-l 30 '+' '.join(mss), \
                log='wscleanBeam-c'+str(c)+'.log', cmdType='wsclean', processors='max')
        s.run(check=True)

        logger.info('Cleaning beam high-res (cycle: '+str(c)+')...')
        imagename = 'img/wideBeamHR'
        s.add('wsclean -reorder -name ' + imagename + ' -size 5500 5500 -mem 90 -j '+str(s.max_processors)+' \
                -scale 4arcsec -weight briggs -1.5 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -apply-primary-beam -use-differential-lofar-beam -minuv-l 30 '+' '.join(mss), \
                log='wscleanBeamHR-c'+str(c)+'.log', cmdType='wsclean', processors='max')
        s.run(check=True)

    # clean mask clean (cut at 5k lambda)
    # no MODEL_DATA update with -baseline-averaging
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.9 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 20 -minuv-l 30 '+' '.join(mss), \
            log='wsclean-c'+str(c)+'.log', cmdType='wsclean', processors='max')
    s.run(check=True)

    maskname = imagename+'-mask.fits'
    make_mask.make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 3, atrous_do=True)
    if user_mask is not None: 
        blank_image_reg(maskname, user_mask, inverse=False, blankval=1)

    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    #TODO: -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,9 \
    #s.add('wsclean -reorder -name ' + imagename + ' -size 3000 3000 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
    s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -niter 1000000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,3,9 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 0.1 -minuv-l 30 -save-source-list -fitsmask '+maskname+' '+' '.join(mss), \
            log='wscleanM-c'+str(c)+'.log', cmdType='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-c'+str(c)+'.log | grep "background noise"')

    if c > 0 and c != niter:
        ft_model_cc(mss, imagename, c, user_mask=user_mask, keep_in_beam=True, model_column='MODEL_DATA')
    # do low-res first cycle and remove it from the data
    if c == 0:
        ft_model_cc(mss, imagename, c, user_mask=user_mask, keep_in_beam=True, model_column='MODEL_DATA_HIGHRES')

        # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_HIGHRES)...')
        s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA_HIGHRES"', log='taql1-c'+str(c)+'.log', cmdType='general')
        s.run(check=True)
    
        # reclean low-resolution
        logger.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        #s.add('wsclean -reorder -name ' + imagename_lr + ' -size 4500 4500 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
        s.add('wsclean -reorder -name ' + imagename_lr + ' -size 6000 6000 -trim 5500 5500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
                -scale 20arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 2000 -mgain 0.8 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 1 -minuv-l 100 -save-source-list '+' '.join(mss), \
                log='wsclean-lr.log', cmdType='wsclean', processors='max')
        s.run(check=True)
       
        ft_model_cc(mss, imagename_lr, 'lr', keep_in_beam=False, model_column='MODEL_DATA')

        # corrupt model with TEC solutions ms:MODEL_DATA -> ms:MODEL_DATA
        if multiepoch: h5 = '$pathMS/tec.h5'
        else: h5 = 'cal-tec'+str(c)+'b.h5'
        MSs.run('DPPP '+parset_dir+'/DPPP-corTEC.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                cor1.parmdb='+h5+' cor1.invert=false cor2.parmdb='+h5+' cor2.invert=false', \
                log=ms+'_corrupt.log', cmdType='DPPP')
    
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
        logger.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='taql2-c'+str(c)+'.log', cmdType='general')
        s.run(check=True)

        # Restore best model
        logger.info('Restoring high-res model (MODEL_DATA = MODEL_DATA_HIGHRES)...')
        s.add('taql "update '+concat_ms+' set MODEL_DATA = MODEL_DATA_HIGHRES"', log='taql3-c'+str(c)+'.log', cmdType='general')
        s.run(check=True)

    ###############################################################################################################
    # Flag on residuals (CORRECTED_DATA)
    #logger.info('Flagging residuals...')
    #for ms in mss:
    #    s.add('DPPP '+parset_dir+'/DPPP-flag.parset msin='+ms, log=ms+'_flag-c'+str(c)+'.log', cmdType='DPPP')
    #s.run(check=True
    
# make beam
# TODO: remove when wsclean will produce a proper primary beam
os.system('~/opt/src/makeavgpb/build/wsbeam.py img/wideBeam')

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources.txt self/images') for c in xrange(niter) ]
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wideBeam-MFS-image.fits  img/wideBeam-MFS-image-pb.fits self/images')
os.system('mv img/wideBeamHR-MFS-image.fits  img/wideBeamHR-MFS-image-pb.fits self/images')
os.system('mv logs self')

logger.info("Done.")
