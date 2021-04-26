#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# they need to be in "./mss/"

import sys, os, glob, re
import numpy as np
import casacore.tables as pt
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-self.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_self','parset_dir')
skydb_demix = parset.get('LOFAR_self','demix_model')
sourcedb = parset.get('model','sourcedb')
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')

#############################################################################
# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')

    # here images, models, solutions for each group will be saved
    lib_util.check_rm('self')
    if not os.path.exists('self/images'): os.makedirs('self/images')
    if not os.path.exists('self/solutions'): os.makedirs('self/solutions')
    if not os.path.exists('self/plots'): os.makedirs('self/plots')

### DONE

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )
try:
    MSs.print_HAcov()
except:
    logger.error('Problem with HAcov, continue anyway.')

# make beam to the first mid null
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg', freq='mid', to_null=True)
beamReg = 'self/beam.reg'

# set image size
imgsizepix = int(2.1*MSs.getListObj()[0].getFWHM(freq='mid')*3600/10.)
if imgsizepix%2 != 0: imgsizepix += 1 # prevent odd img sizes

# set clean componet fit order (use 5 for large BW)
bandwidth = MSs.getBandwidth()
if bandwidth > 25e6: cc_fit_order = 5
else: cc_fit_order = 3

#################################################################
# Get online model
if sourcedb == '':
    if not os.path.exists('tgts.skydb'):
        fwhm = MSs.getListObj()[0].getFWHM(freq='min')
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<1')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')
        apparent = False
    sourcedb = 'tgts.skydb'

#################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS + '/' + sourcedb_basename)
    logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
    os.system('cp -r ' + sourcedb + ' ' + MS)

with w.if_todo('init_model'):

    # note: do not add MODEL_DATA or the beam is transported from DATA, while we want it without beam applied
    logger.info('Creating CORRECTED_DATA...')
    MSs.addcol('CORRECTED_DATA', 'DATA')

    logger.info('Add model to MODEL_DATA...')
    if apparent:
        MSs.run(
            'DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=false pre.sourcedb=$pathMS/' + sourcedb_basename,
            log='$nameMS_pre.log', commandType='DP3')
    else:
        MSs.run(
            'DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=true pre.sourcedb=$pathMS/' + sourcedb_basename,
            log='$nameMS_pre.log', commandType='DP3')
### DONE

#####################################################################################################
# Self-cal cycle
for c in range(2):

    logger.info('Start selfcal cycle: '+str(c))

    if c == 0:
        with w.if_todo('set_corrected_data'):
            logger.info('Set CORRECTED_DATA = DATA...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE
    else:
        with w.if_todo('cor_g_c%02i' % c):
            # initial correction for slow gain after first cycle
            # correct G - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
            logger.info('Correcting G...')
            MSs.run(
                'DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=self/solutions/cal-g-c0.h5 cor.correction=amplitudeSmooth', \
                log='$nameMS_corG-c' + str(c) + '.log', commandType='DP3')
    if c == 0:
        with w.if_todo('solve_fr_c%02i' % c):
            logger.info('Add column CIRC_PHASEDIFF_DATA...')
            MSs.addcol('CIRC_PHASEDIFF_DATA', 'CORRECTED_DATA', usedysco=False)
            # Probably we do not need smoothing since we have long time intervals and smoothnessconstraint?
            # logger.info('BL-smooth...')
            # MSs.run('BLsmooth.py -r -c 1 -n 8 -i CIRC_PHASEDIFF_DATA -o CIRC_PHASEDIFF_DATA $pathMS',
            #         log='$nameMS_smooth2-c' + str(c) + '.log',
            #         commandType='python', maxThreads=8)

            logger.info('Converting to circular...')
            MSs.run('mslin2circ.py -s -i $pathMS:CIRC_PHASEDIFF_DATA -o $pathMS:CIRC_PHASEDIFF_DATA',
                    log='$nameMS_lincirc.log',
                    commandType='python', maxThreads=2)

            # Get circular phase diff CIRC_PHASEDIFF_DATA -> CIRC_PHASEDIFF_DATA
            logger.info('Get circular phase difference...')
            MSs.run('taql "UPDATE $pathMS SET\
             CIRC_PHASEDIFF_DATA[,0]=0.5*EXP(1.0i*(PHASE(CIRC_PHASEDIFF_DATA[,0])-PHASE(CIRC_PHASEDIFF_DATA[,3]))), \
             CIRC_PHASEDIFF_DATA[,3]=CIRC_PHASEDIFF_DATA[,0], \
             CIRC_PHASEDIFF_DATA[,1]=0+0i, \
             CIRC_PHASEDIFF_DATA[,2]=0+0i"', log='$nameMS_taql_phdiff.log', commandType='general')

            logger.info('Creating FR_MODEL_DATA...')  # take from MODEL_DATA but overwrite
            MSs.addcol('FR_MODEL_DATA', 'MODEL_DATA', usedysco=False)
            MSs.run('taql "UPDATE $pathMS SET FR_MODEL_DATA[,0]=0.5+0i, FR_MODEL_DATA[,1]=0.0+0i, FR_MODEL_DATA[,2]=0.0+0i, \
             FR_MODEL_DATA[,3]=0.5+0i"', log='$nameMS_taql_frmodel.log', commandType='general')

            # Solve cal_SB.MS:CIRC_PHASEDIFF_DATA against FR_MODEL_DATA (only solve)
            logger.info('Solving circ phase difference ...')
            MSs.run('DP3 ' + parset_dir + '/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5',
                    log='$nameMS_solFR-c' + str(c) + '.log', commandType="DP3")
            lib_util.run_losoto(s, f'fr-c{c}', [ms + '/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
            os.system('mv cal-fr-c' + str(c) + '.h5 self/solutions/')
            os.system('mv plots-fr-c' + str(c) + ' self/plots/')
            # Delete cols again to not waste space
            MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, FR_MODEL_DATA"',
                    log='$nameMS_taql_delcol.log', commandType='general')
        ### DONE

    with w.if_todo('cor_fr_c%02i' % c):
        # Correct FR with results of cycle 0 - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting FR...')
        MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                cor.parmdb=self/solutions/cal-fr-c0.h5 cor.correction=rotationmeasure000',
                log='$nameMS_corFR-c' + str(c) + '.log', commandType='DP3')
    ### DONE

    with w.if_todo('solve_tec1_c%02i' % c):
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -c 8 -n 8 -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')
        MSs.run('BLsmooth.py -c 8 -n 8 -r -i MODEL_DATA -o MODEL_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')

        # solve TEC - ms:SMOOTHED_DATA
        logger.info('Solving TEC1...')
        MSs.run('DP3 '+parset_dir+'/DP3-solTEC.parset msin=$pathMS sol.h5parm=$pathMS/tec1.h5 \
                msin.baseline="[CR]*&&;!RS208LBA;!RS210LBA;!RS307LBA;!RS310LBA;!RS406LBA;!RS407LBA;!RS409LBA;!RS508LBA;!RS509LBA" \
                sol.antennaconstraint=[[CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA]] \
           	    sol.solint=15 sol.nchan=8', \
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec1-c'+str(c), [ms+'/tec1.h5' for ms in MSs.getListStr()], [parset_dir+'/losoto-resetremote.parset', parset_dir+'/losoto-plot-tec.parset'])
        os.system('mv cal-tec1-c'+str(c)+'.h5 self/solutions/')
        os.system('mv plots-tec1-c'+str(c)+' self/plots/')
    ### DONE

    with w.if_todo('cor_tec1_c%02i' % c):
        # correct TEC - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting TEC1...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA\
                cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=tec000 cor.missingantennabehavior=unit',
                log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DP3')
    ### DONE

    with w.if_todo('solve_tec2_c%02i' % c):
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -c 8 -n 8 -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')

        # solve TEC - ms:SMOOTHED_DATA
        logger.info('Solving TEC2...')
        MSs.run('DP3 '+parset_dir+'/DP3-solTEC.parset msin=$pathMS sol.h5parm=$pathMS/tec2.h5 \
                sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LBA]] \
                sol.solint=1 sol.nchan=4',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec2-c'+str(c), [ms+'/tec2.h5' for ms in MSs.getListStr()], [parset_dir+'/losoto-plot-tec.parset'])
        os.system('mv cal-tec2-c'+str(c)+'.h5 self/solutions/')
        os.system('mv plots-tec2-c'+str(c)+' self/plots/')
    ### DONE

    with w.if_todo('cor_tec2_c%02i' % c):
        # correct TEC - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting TEC2...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA\
                cor.parmdb=self/solutions/cal-tec2-c'+str(c)+'.h5 cor.correction=tec000',
                log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DP3')
    ### DONE

    # AMP DIE correction
    if c == 0:
        with w.if_todo('solve_g_c%02i' % c):
            # DIE Calibration - ms:CORRECTED_DATA
            logger.info('Solving slow G...')
            MSs.run('DP3 '+parset_dir+'/DP3-solG.parset msin=$pathMS sol.h5parm=$pathMS/g.h5',
                    log='$nameMS_solG-c'+str(c)+'.log', commandType='DP3')
            lib_util.run_losoto(s, 'g-c'+str(c), [MS+'/g.h5' for MS in MSs.getListStr()],
                    [parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-amp.parset'])
            os.system('mv plots-g-c'+str(c)+' self/plots/')
            os.system('mv cal-g-c'+str(c)+'.h5 self/solutions/')
        ### DONE

        with w.if_todo('cor_g_c%02i' % c):
            # correct G - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
            logger.info('Correcting G...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                    cor.parmdb=self/solutions/cal-g-c'+str(c)+'.h5 cor.correction=amplitudeSmooth',
                    log='$nameMS_corG-c'+str(c)+'.log', commandType='DP3')
        ### DONE

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA

    imagename = 'img/wide-0'
    maskname = imagename + '-mask.fits'
    imagenameM = 'img/wideM-'+str(c)
    with w.if_todo('imaging_c%02i' % c):
        logger.info('Cleaning (cycle: '+str(c)+')...')
        if c == 0:
            # make temp mask for cycle 0, in cycle 1 use the maske made from cycle 0 image
            lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), name=imagename,
                                 size=imgsizepix, scale='10arcsec',
                                 weight='briggs -0.3', niter=1000000, no_update_model_required='', minuv_l=30,
                                 parallel_gridding=2, baseline_averaging='', maxuv_l=4500, mgain=0.85,
                                 parallel_deconvolution=512, local_rms='', auto_threshold=4,
                                 join_channels='', fit_spectral_pol=cc_fit_order, channels_out=MSs.getChout(4.e6),
                                 deconvolution_channels=cc_fit_order)
            im = lib_img.Image(imagename + '-MFS-image.fits', userReg=userReg)
            im.makeMask(threshpix=5)

            kwargs = {'do_predict':True, 'reuse_dirty':imagename, 'reuse_psf':imagename}
        else: 
            kwargs = {}

        lib_util.run_wsclean(s, 'wscleanM-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM, save_source_list='',
                size=imgsizepix, scale='10arcsec',
                weight='briggs -0.3', niter=1000000, no_update_model_required='', minuv_l=30,
                parallel_gridding=2, baseline_averaging='', maxuv_l=4500, mgain=0.85,
                parallel_deconvolution=512, auto_threshold=3., fits_mask=maskname,
                join_channels='', fit_spectral_pol=cc_fit_order, channels_out=MSs.getChout(4.e6),
                multiscale='', multiscale_scale_bias=0.6,
                deconvolution_channels=cc_fit_order, **kwargs)

        os.system('cat logs/wscleanM-c'+str(c)+'.log | grep "background noise"')
    ### DONE

    if c == 0:

        with w.if_todo('lowres_setdata_c%02i' % c):
            # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
            logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE
    
        with w.if_todo('lowres_img_c%02i' % c):
            # Making beam mask
            logger.info('Preparing mask for low-res clean...')
            lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name='img/tmp', size=imgsizepix, scale='30arcsec')
            os.system('mv img/tmp-image.fits img/wide-lr-mask.fits')
            lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 0.)
            lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 1., inverse=True)
    
            # reclean low-resolution
            logger.info('Cleaning low-res...')
            imagename_lr = 'img/wide-lr'
            lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=False,
                    parallel_gridding=4, temp_dir='./', size=imgsizepix, scale='30arcsec',
                    weight='briggs -0.3', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
                    taper_gaussian='200arcsec', mgain=0.85, parallel_deconvolution=512, baseline_averaging='',
                    local_rms='', auto_mask=3, auto_threshold=1.5, fits_mask='img/wide-lr-mask.fits',
                    join_channels='', channels_out=MSs.getChout(2.e6))

            s.add('wsclean -predict -name '+imagename_lr+' -j '+str(s.max_processors)+' -channels-out '+str(MSs.getChout(2e6))+' '+MSs.getStrWsclean(), \
                  log='wscleanLR-PRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
            s.run(check=True)
        ### DONE

        with w.if_todo('lowres_sub_c%02i' % c):
            # Subtract low-res model - CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA
            logger.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE

        with w.if_todo('lowres_lsimg_c%02i' % c):
            logger.info('Cleaning large-scale...')
            imagename_ls = 'img/wide-largescale'
            #                     intervals_out=len(MSs.mssListObj)*4,
            #use_idg = '', aterm_kernel_size = 16, aterm_config = parset_dir + '/aconfig.txt',
            lib_util.run_wsclean(s, 'wscleanLS.log', MSs.getStrWsclean(), name=imagename_ls, do_predict=False,
                                 temp_dir='./', size=2000, scale='20arcsec',
                                 no_fit_beam='', circular_beam='', beam_size='180.0arcsec',
                                 multiscale='', multiscale_scales='0,4,8,16,32,64',
                                 weight='briggs -0.3', niter=10000, no_update_model_required='', minuv_l=20,
                                 maxuvw_m=5000, taper_gaussian='180arcsec', mgain=0.85,
                                 parallel_deconvolution=512, baseline_averaging='', local_rms='', auto_mask=1.5,
                                 auto_threshold=0.5, join_channels='', channels_out=MSs.getChout(4.e6))
        ### DONE

        with w.if_todo('lowres_flag_c%02i' % c):
            # Flag on residuals (CORRECTED_DATA)
            logger.info('Flagging residuals...')
            MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.lua',
                    log='$nameMS_flag-c'+str(c)+'.log', commandType='DP3')
        ### DONE

        with w.if_todo('lowres_corrupt_c%02i' % c):
            ##############################################
            # Prepare SUBTRACTED_DATA
    
            # corrupt model with TEC+FR+Beam2ord solutions - ms:MODEL_DATA -> ms:MODEL_DATA
            logger.info('Corrupt low-res model: TEC1...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                    cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
            logger.info('Corrupt low-res model: TEC2...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                    cor.parmdb=self/solutions/cal-tec2-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                    cor.parmdb=self/solutions/cal-fr-c'+str(c)+'.h5 cor.correction=rotationmeasure000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
            logger.info('Corrupt low-res model: G...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                    cor.parmdb=self/solutions/cal-g-c'+str(c)+'.h5 cor.correction=amplitudeSmooth cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
        ### DONE

        with w.if_todo('lowres_subtract_c%02i' % c):
            # Subtract low-res model - CORRECTED_DATA = DATA - MODEL_DATA
            logger.info('Subtracting low-res model (CORRECTED_DATA = DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE

        with w.if_todo('lowres_predict_c%02i' % c):
            # Recreate MODEL_DATA for next calibration cycle
            logger.info('Predict model...')
            s.add('wsclean -predict -name img/wideM-'+str(c)+' -j '+str(s.max_processors)+' -channels-out '+str(MSs.getChout(4e6))+' '+MSs.getStrWsclean(), \
                   log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
            s.run(check=True)
        ### DONE

# polarisation imaging
with w.if_todo('imaging-pol'):
    logger.info('Cleaning (Pol)...')
    imagenameP = 'img/wideP'
    lib_util.run_wsclean(s, 'wscleanP.log', MSs.getStrWsclean(), name=imagenameP, pol='QUV',
        size=imgsizepix, scale='10arcsec', weight='briggs -0.3', niter=0, no_update_model_required='',
        parallel_gridding=2, baseline_averaging='', minuv_l=30, maxuv_l=4500,
        join_channels='', channels_out=MSs.getChout(4.e6))

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image*.fits self/images') for c in range(2) ]
[ os.system('mv img/wideM-'+str(c)+'-MFS-residual*.fits self/images') for c in range(2) ]
[ os.system('mv img/wideM-'+str(c)+'-sources*.txt self/images') for c in range(2) ]
os.system('mv img/wideP-MFS-*-image.fits self/images')
os.system('mv img/wideM-1-*-model.fits self/images')
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wide-largescale-MFS-image.fits self/images')
os.system('mv logs self')

logger.info("Done.")
