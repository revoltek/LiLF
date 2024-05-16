#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# they need to be in "./mss/"

import sys, os, glob
import numpy as np
from regions import Regions
import lsmtool
from astropy.coordinates import Angle

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-self')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-self.walker')

parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_self'])))
parset_dir = parset.get('LOFAR_self','parset_dir')
subfield = parset.get('LOFAR_self','subfield') # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
maxIter = parset.getint('LOFAR_self','maxIter') # default = 2 (try also 3)
phaseSolMode = parset.get('LOFAR_self', 'ph_sol_mode') # tecandphase, tec, phase
#sourcedb = parset.get('model','sourcedb')
sourcedb = 'tgts.skymodel'
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')

#############################################################################

def clean_self(c, MSs, MSsClean, imgsizepix, cc_fit_order, subfield_kwargs, do_predict=True):
    """ Do normal (hres) imaging of the self-calibrated data. """
    imagename = 'img/wide-' + str(c)
    imagenameM = 'img/wideM-' + str(c)

    # first clean using local_rms and higher thresholds, -> update MODEL_DATA # concat_mss=True, baseline_averaging='',
    lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSsClean.getStrWsclean(), concat_mss=True, keep_concat=True,
                         name=imagenameM, size=imgsizepix, scale='4arcsec',
                         weight='briggs -0.3', niter=1000000, update_model_required='', minuv_l=30,
                         parallel_gridding=6,  mgain=0.8, parallel_deconvolution=1024,
                         auto_mask=6, auto_threshold=4.5, local_rms='', join_channels='', fit_spectral_pol=cc_fit_order,
                         channels_out=MSsClean.getChout(4.e6), multiscale='', multiscale_scale_bias=0.6,
                         deconvolution_channels=cc_fit_order, **subfield_kwargs)
    if userReg is not None:
        s.add('breizorro.py -t 6.5 -r %s -b 50 -o %s --merge %s' % (
        imagenameM + '-MFS-image.fits', imagename + '-mask.fits', userReg),
              log=f'makemask-{c}.log', commandType='python')
        s.run()
    else:
        s.add('breizorro.py -t 6.5 -r %s -b 50 -o %s' % (imagenameM + '-MFS-image.fits', imagename + '-mask.fits'),
              log=f'makemask-{c}.log', commandType='python')
        s.run()
    os.system(f'cp {imagenameM}-MFS-image.fits {imagename}-MFS-image.fits')
    os.system(f'cp {imagenameM}-MFS-residual.fits {imagename}-MFS-residual.fits')
    os.system(f'cp {imagenameM}-MFS-model.fits {imagename}-MFS-model.fits')
    # second continue clean using mask and lower threshold -> update MODEL_DATA
    # TODO check save_source_listbaseline_averaging='',
    lib_util.run_wsclean(s, 'wscleanM-c' + str(c) + '.log', MSsClean.getStrWsclean(), concat_mss=True, fits_mask=imagename+'-mask.fits',
                         name=imagenameM, reuse_psf=imagenameM, reuse_dirty=imagenameM, cont='', save_source_list='', size=imgsizepix,
                         scale='4arcsec', weight='briggs -0.3', niter=1000000, no_update_model_required='', minuv_l=30,
                         parallel_gridding=6,  mgain=0.8, parallel_deconvolution=1024,
                         auto_threshold=3., join_channels='', fit_spectral_pol=cc_fit_order, channels_out=MSsClean.getChout(4.e6),
                         multiscale='', multiscale_scale_bias=0.6, deconvolution_channels=cc_fit_order, **subfield_kwargs)

    os.system('cat ' + logger_obj.log_dir + '/wscleanM-c' + str(c) + '.log | grep "background noise"')

    # when wsclean allows station selection, then we can remove MSsClean and this predict can go in the previous call with do_predict=True
    if do_predict:
        logger.info('Predict model...')
        pred_str = f'wsclean -predict -padding 1.8 -name img/wideM-{c} -j {s.max_processors} -channels-out {MSs.getChout(4e6)}'
        if 'shift' in subfield_kwargs:
            pred_str += f' -shift {shift}'
        s.add(f'{pred_str} {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
        s.run(check=True)

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
    if not os.path.exists('self/skymodel'): os.makedirs('self/skymodel')

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
fov_maj = 2*MSs.getListObj()[0].getFWHM(freq='mid', elliptical=True)[0]

# set image size
imgsizepix_wide = int(2.1*MSs.getListObj()[0].getFWHM(freq='mid')*3600/4.)
imgsizepix_lr = int(5*MSs.getListObj()[0].getFWHM(freq='mid')*3600/30.)
imgsizepix_p = int(2.1*MSs.getListObj()[0].getFWHM(freq='mid')*3600/10.)

# set clean componet fit order (use 5 for large BW)
if MSs.getChout(4.e6) >= 7:  # Bandwidth of 28 MHz or more
    cc_fit_order = 5
else: cc_fit_order = 3

fullband = MSs.getBandwidth()
nchan = MSs.mssListObj[0].getNchan()
tint = MSs.mssListObj[0].getTimeInt()
if int(np.rint(fullband / nchan < 195.3e3/4)):
    base_nchan = int(np.rint((195.3e3/4)/(fullband/nchan))) # this is 1 for ducth observations, and larger (2,4) for IS observations
else: base_nchan = 1
if MSs.hasIS:
    base_solint = 1
elif tint < 4:
    base_solint = int(np.rint(4/tint)) # this is 2 for dutch SPARSE observations
else: base_solint = 1

#################################################################
# make filter mask
imagename_fov = 'img/fov_mask'
# lib_util.run_wsclean(s, 'wsclean_temp.log', MSs.getStrWsclean(), name='img/fov_mask', maxuv_l=1000,
#                      size=int(1.2*fov_maj*3600/30), scale='30arcsec')
# os.system('mv img/fov_mask-dirty.fits img/fov_mask.fits')
# lib_img.blank_image_reg(imagename_fov+'.fits', beamReg, blankval=1.)
# lib_img.blank_image_reg(imagename_fov+'.fits', beamReg, blankval=0., inverse=True)

# prefix = '/beegfs/p1uy068/virgo/models/LBA'
# sm = lsmtool.load(f'{prefix}/LVCS_20as_gaul_filtered_freqscaled.skymodel', beamMS=MSs.getListStr()[int(len(MSs.getListStr())/2)])
# sm.select('img/fov_mask.fits=1') # remove distant sources
# sm.select('I>0.7', applyBeam=True) # keep only reasonably bright sources
# sm.write('tgts-pre.skymodel', clobber=True, applyBeam=True, adjustSI=True)
# sm = lsmtool.load('tgts-pre.skymodel')
# sm.group('cluster', numClusters=10) #, applyBeam=True)
# sm.setPatchPositions(method='wmean') #, applyBeam=True)
# sm.plot('self/skymodel/patches.png', 'patch')
# sm.write(sourcedb, clobber=True)

# m87_dist = lib_util.distanceOnSphere(ra, dec, 187.705930, 12.391123)
# if m87_dist > 4:
#     print(f'm87_dist-{m87_dist}deg. Assume M87 demixed and not put it in model.')
# else:
#     m87 = lsmtool.load(f'{prefix}/VirA.skymodel', beamMS=args.MS)
#     field_hba.concatenate(m87)
#
# print(field_hba.info())
# field_hba.group('single', root='target')
#
# field_hba.write('tgts.skymodel', clobber=True, applyBeam=False)
#################################################################################################
# goes down to 8 seconds and multiple of 48 chans (this should be already the case as it's done in timesplit)

# Here the model is added only to CS+RS, IS used only for FR and model is not needed
# with w.if_todo('init_model'):
#     logger.info('Add model to MODEL_DATA...')
#     if apparent:
#         MSs.run('DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=false pre.sourcedb=$pathMS/' + sourcedb_basename,
#             log='$nameMS_pre.log', commandType='DP3')
#     else:
#         MSs.run('DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=true pre.usechannelfreq=True \
#                  pre.beammode=array_factor pre.onebeamperpatch=True pre.sourcedb=$pathMS/' + sourcedb_basename,
#                 log='$nameMS_pre.log', commandType='DP3')
# ### DONE

# TODO MSs_avg
with w.if_todo('solve_cor_fr'):
    logger.info('Add column CIRC_PHASEDIFF_DATA...')
    MSs.addcol('CIRC_PHASEDIFF_DATA', 'DATA', usedysco=False)  # No dysco here as off diag elements are 0 and dysco does not like 0s
    # Probably we do not need smoothing since we have long time intervals and smoothnessconstraint?

    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -s -i $pathMS:CIRC_PHASEDIFF_DATA -o $pathMS:CIRC_PHASEDIFF_DATA',
            log='$nameMS_lincirc.log', commandType='python', maxThreads=2)

    # Get circular phase diff CIRC_PHASEDIFF_DATA -> CIRC_PHASEDIFF_DATA
    logger.info('Get circular phase difference...')
    MSs.run('taql "UPDATE $pathMS SET\
         CIRC_PHASEDIFF_DATA[,0]=0.5*EXP(1.0i*(PHASE(CIRC_PHASEDIFF_DATA[,0])-PHASE(CIRC_PHASEDIFF_DATA[,3]))), \
         CIRC_PHASEDIFF_DATA[,3]=CIRC_PHASEDIFF_DATA[,0], \
         CIRC_PHASEDIFF_DATA[,1]=0+0i, \
         CIRC_PHASEDIFF_DATA[,2]=0+0i"', log='$nameMS_taql_phdiff.log', commandType='general')

    logger.info('Creating FR_MODEL_DATA...')  # take from MODEL_DATA but overwrite
    MSs.addcol('FR_MODEL_DATA', 'DATA', usedysco=False)
    MSs.run('taql "UPDATE $pathMS SET FR_MODEL_DATA[,0]=0.5+0i, FR_MODEL_DATA[,1]=0.0+0i, FR_MODEL_DATA[,2]=0.0+0i, \
         FR_MODEL_DATA[,3]=0.5+0i"', log='$nameMS_taql_frmodel.log', commandType='general')

    # Solve cal_SB.MS:CIRC_PHASEDIFF_DATA against FR_MODEL_DATA (only solve - solint=2m nchan=0 as it has the smoothnessconstrain)
    logger.info('Solving circ phase difference ...')
    MSs.run('DP3 ' + parset_dir + '/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.solint=' + str(30 * base_solint),
            log='$nameMS_solFR.log', commandType="DP3")
    lib_util.run_losoto(s, f'fr', [ms + '/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
    os.system('mv cal-fr.h5 self/solutions/')
    os.system('mv plots-fr self/plots/')
    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, FR_MODEL_DATA"',
            log='$nameMS_taql_delcol.log', commandType='general')

    # Correct FR with results of solve - group*_TC.MS:DATA -> group*_TC.MS:FR_CORRECTED_DATA
    logger.info('Correcting FR...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=FR_CORRECTED_DATA \
                cor.parmdb=self/solutions/cal-fr.h5 cor.correction=rotationmeasure000',
            log='$nameMS_corFR.log', commandType='DP3')
### DONE
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
# sourcedb_basename = sourcedb.split('/')[-1]
# for MS in MSs_avg.getListStr():
#     lib_util.check_rm(MS + '/' + sourcedb_basename)
#     logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
#     os.system('cp -r ' + sourcedb + ' ' + MS)
#####################################################################################################
# Self-cal cycle
for c in range(maxIter):
    logger.info('Start selfcal cycle: '+str(c))
    if c == 0:
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Creating CORRECTED_DATA from FR_CORRECTED_DATA...')
            MSs.addcol('CORRECTED_DATA', 'FR_CORRECTED_DATA')
    else:
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Set CORRECTED_DATA = SUBFIELD_DATA...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBFIELD_DATA"', log='$nameMS_taql-c' + str(c) + '.log',
                    commandType='general')
    ### DONE

    # with w.if_todo('smooth_model_c%02i' % c):
    #     # Smooth MODEL_DATA -> MODEL_DATA
    #     MSs.run_Blsmooth('MODEL_DATA', 'MODEL_DATA', logstr=f'smooth-c{c}')
    ### DONE

    with w.if_todo('sol_avg_%02i' % c):
        lib_util.check_rm('mss-sol')
        os.system('mkdir mss-sol')
        timeint = MSs.getListObj()[0].getTimeInt()
        avgtimeint = int(round(8 / timeint))  # to 16 seconds
        nchan_init = MSs.getListObj()[0].getNchan()
        # chan: avg (x8) sol (x6) - we need a multiple of 8x6=48, the largest that is <nchan
        # survey after avg (x8): 60, final number of sol 10
        # pointed after avg (x8): 120, final number of sol 20
        nchan = nchan_init - nchan_init % 48
        logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (
            timeint, timeint * avgtimeint, nchan_init, nchan/8))
        MSs.run(
            'DP3 ' + parset_dir + '/DP3-avg.parset msin=$pathMS msout=mss-sol/$nameMS.MS msin.datacolumn=CORRECTED_DATA msin.nchan=' + str(
                nchan) + ' \
                avg.timestep=' + str(avgtimeint) + ' avg.freqstep=8',
            log='$nameMS_initavg.log', commandType='DP3')

    MSs_sol = lib_ms.AllMSs(glob.glob('mss-sol/TC*[0-9].MS'), s, check_flags=True)

    with w.if_todo('solve_tec_c%02i' % c):
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        MSs_sol.run_Blsmooth('DATA', logstr=f'smooth-c{c}')

        # solve ionosphere phase - ms:SMOOTHED_DATA (1m 2SB)
        logger.info('Solving TEC1...')
        if phaseSolMode == 'phase': #phase
            solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=0.5e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
            losoto_parsets = [parset_dir+'/losoto-plot-scalar.parset']
        else: # TEC or TecAndPhase
            solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
            losoto_parsets = [parset_dir+'/losoto-plot-tec.parset']

        MSs_sol.run(f"DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec1.h5 sol.solint=1 {solver_params} sol.sourcedb=tgts.skymodel",
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec1-c'+str(c), [ms+'/tec1.h5' for ms in MSs_sol.getListStr()], losoto_parsets)
        os.system('mv cal-tec1-c'+str(c)+'.h5 self/solutions/')
        os.system('mv plots-tec1-c'+str(c)+' self/plots/')
    ### DONE

    facetregname = 'self/solutions/facets.reg'
    # TDOO temp
    MSs = MSs_sol
    with w.if_todo('c%02i-imaging' % c):
        logger.info('Preparing region file...')
        s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+f' --h5 self/solutions/cal-tec1-c{c}.h5 --imsize '+str(imgsizepix_wide)+' \
                --pixelscale 4 --writevoronoipoints --output '+facetregname,
              log='facet_generator.log', commandType='python')
        s.run()

        imagename = 'img/wide-' + str(c)
        imagenameM = 'img/wideM-' + str(c)
        # might want to add non-circ beam at low dec eventually

        # update_model=True to make continue after masking  dd_psf_grid='25 25', beam_size=15,CORRECTED_
        logger.info('Cleaning 1...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM, data_column='DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder', parallel_gridding=6, update_model_required='', minuv_l=30, mgain=0.8, parallel_deconvolution=1024,
                             auto_threshold=3.0, auto_mask=5.0, join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='', pol='i', nmiter=3,
                             facet_regions=facetregname,apply_facet_solutions=f'self/solutions/cal-tec1-c{c}.h5 phase000' )

        # masking
        if userReg is not None:
            s.add('breizorro.py -t 6.5 -r %s -b 50 -o %s --merge %s' % (
                imagenameM + '-MFS-image.fits', imagename + '-mask.fits', userReg),
                  log=f'makemask-{c}.log', commandType='python')
            s.run()
        else:
            s.add('breizorro.py -t 6.5 -r %s -b 50 -o %s' % (imagenameM + '-MFS-image.fits', imagename + '-mask.fits'),
                  log=f'makemask-{c}.log', commandType='python')
            s.run()
        os.system(f'cp {imagenameM}-MFS-image.fits {imagename}-MFS-image.fits')
        os.system(f'cp {imagenameM}-MFS-residual.fits {imagename}-MFS-residual.fits')
        os.system(f'cp {imagenameM}-MFS-model.fits {imagename}-MFS-model.fits')

        # TODO: Add force-reuse psf (and beam) once wsclean bug is fixed.
        # HE: What is optimal choice of subimage size and parallel gridding? Is cleaning to 3sigma enough? # CORRECTED_DATA
        logger.info('Cleaning 2...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM, fits_mask=imagename+'-mask.fits', data_column='DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder', cont=True, reuse_psf=imagenameM, reuse_dirty=imagenameM, parallel_gridding=6,
                             update_model_required='', minuv_l=30, mgain=0.8, parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=5.0,
                             join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='', pol='i', facet_regions=facetregname, apply_facet_solutions=f'self/solutions/cal-tec1-c{c}.h5 phase000' )
        sys.exit()

    imagenameM = 'img/wideM-'+str(c)
    if (c == 0) or (c == maxIter - 1): # make wide image
        imgsizepix = imgsizepix_wide
        subfield_kwargs = {}
    else: # make small image (field subtracted)
        imgsizepix = int(1.3 * field_size * 3600 / 4.)
        img_ra, img_dec = Angle(f'{field_center[0]}deg').hms, Angle(f'{field_center[1]}deg').dms
        shift = f'{int(img_ra.h)}h{int(img_ra.m)}m{img_ra.s:.4f}s {int(img_dec.d)}d{int(img_dec.m)}m{img_dec.s:.4f}s'
        subfield_kwargs = {'shift': shift}

    if c < maxIter - 1:
        with w.if_todo('imaging_c%02i' % c):
            clean_self(c, MSs, MSsClean, imgsizepix, cc_fit_order, subfield_kwargs)
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
            lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name='img/tmp', size=imgsizepix_lr, scale='30arcsec')
            os.system('mv img/tmp-image.fits img/wide-lr-mask.fits')
            lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 0.)
            lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 1., inverse=True)
    
            # reclean low-resolution
            logger.info('Cleaning low-res...')
            imagename_lr = 'img/wide-lr'
            lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=True,
                                 parallel_gridding=4, temp_dir='../', size=imgsizepix_lr, scale='30arcsec',
                                 weight='briggs -0.3', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
                                 taper_gaussian='200arcsec', mgain=0.85, parallel_deconvolution=512, baseline_averaging='',
                                 local_rms='', auto_mask=3, auto_threshold=1.5, fits_mask='img/wide-lr-mask.fits',
                                 join_channels='', channels_out=MSs.getChout(2.e6))
            # Test of we cam just do a do_predict
            # s.add('wsclean -predict -padding 1.8 -name '+imagename_lr+' -j '+str(s.max_processors)+' -channels-out '+str(MSs.getChout(2e6))+' '+MSs.getStrWsclean(), \
            #       log='wscleanLR-PRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
            # s.run(check=True)
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
                                 temp_dir='../', size=2000, scale='20arcsec',
                                 no_fit_beam='', circular_beam='', beam_size='200arcsec',
                                 multiscale='', multiscale_scales='0,4,8,16,32,64',
                                 weight='briggs -0.3', niter=10000, no_update_model_required='', minuv_l=20,
                                 maxuvw_m=5000, taper_gaussian='200arcsec', mgain=0.85,
                                 parallel_deconvolution=512, baseline_averaging='', local_rms='', auto_mask=1.5,
                                 auto_threshold=0.5, join_channels='', channels_out=MSs.getChout(4.e6))
        ### DONE

        with w.if_todo('lowres_corrupt_c%02i' % c):    
            # corrupt model with TEC+FR+Beam2ord solutions - ms:MODEL_DATA -> ms:MODEL_DATA
            if phaseSolMode in ['tec', 'tecandphase']:
                logger.info('Corrupt low-res model: TEC1...')
                MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                        cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False',
                        log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                logger.info('Corrupt low-res model: ph1...')
                MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                    cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=phase000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
            # logger.info('Corrupt low-res model: TEC+Ph 2...')
            # MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
            #         cor.parmdb=self/solutions/cal-tec2-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False',
            #         log='$nameMS_corrupt.log', commandType='DP3')
            # MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
            #         cor.parmdb=self/solutions/cal-tec2-c'+str(c)+'.h5 cor.correction=phase000 cor.invert=False',
            #         log='$nameMS_corrupt.log', commandType='DP3')

        with w.if_todo('lowres_subtract_c%02i' % c):
            # Permanently subtract low-res sidelobe model - FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA.
            logger.info('Subtracting low-res sidelobe model (FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA"',
                    log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE

        # Prepare region and models for subfield
        if subfield:
            subfield_path = subfield
            if len(Regions.read(subfield_path)) > 1:
                raise ValueError(f'Manual subfield region {subfield} contains more than one region')
        else:
            subfield_path = 'self/skymodel/subfield.reg'

        with w.if_todo('extreg_preapre_c%02i' % c):
            if not subfield: # automatically find subfield
                sm = lsmtool.load(f'img/wideM-{c}-sources.txt')
                sm.remove('img/wide-lr-mask.fits=1')  # remove sidelobe sources that were subtracted
                sm.remove('MajorAxis > 80')  # remove largest scales
                field_center1, field_size1 = lib_dd.make_subfield_region(subfield_path, MSs.getListObj()[0], sm,
                                                                         subfield_min_flux, debug_dir='img/')
            # prepare model of central/external regions
            logger.info('Blanking central region of model files and reverse...')
            for im in glob.glob('img/wideM-0*model.fits'):
                wideMint = im.replace('wideM','wideMint')
                os.system('cp %s %s' % (im, wideMint))
                lib_img.blank_image_reg(wideMint, subfield_path, blankval = 0., inverse=True)
                wideMext = im.replace('wideM','wideMext')
                os.system('cp %s %s' % (im, wideMext))
                lib_img.blank_image_reg(wideMext, subfield_path, blankval = 0.)
        # DONE
        subfield_reg = Regions.read(subfield_path)[0]
        field_center = subfield_reg.center.ra.deg, subfield_reg.center.dec.deg
        field_size = np.max([subfield_reg.width.to_value('deg'), subfield_reg.height.to_value('deg')])

        with w.if_todo('extreg_predict_corrupt_subtract_c%02i' % c):
            # Recreate MODEL_DATA of external region for subtraction
            logger.info('Predict model of external region...')
            s.add('wsclean -predict -padding 1.8 -name img/wideMext-'+str(c)+' -j '+str(s.max_processors)+' -channels-out '+str(MSs.getChout(4e6))+' '+MSs.getStrWsclean(), \
                  log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
            s.run(check=True)

            # corrupt model with TEC+FR+Beam2ord solutions - ms:MODEL_DATA -> ms:MODEL_DATA
            if phaseSolMode in ['tec', 'tecandphase']:
                logger.info('Corrupt low-res model: TEC1...')
                MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                        cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False',
                        log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                logger.info('Corrupt low-res model: ph1...')
                MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                        cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=phase000 cor.invert=False',
                        log='$nameMS_corrupt.log', commandType='DP3')
            # logger.info('Corrupt low-res model: TEC+Ph 2...')
            # MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
            #         cor.parmdb=self/solutions/cal-tec2-c'+str(c)+'.h5 cor.correction=tec000 cor.invert=False',
            #         log='$nameMS_corrupt.log', commandType='DP3')
            # MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
            #         cor.parmdb=self/solutions/cal-tec2-c'+str(c)+'.h5 cor.correction=phase000 cor.invert=False',
            #         log='$nameMS_corrupt.log', commandType='DP3')
            # logger.info('Corrupt low-res model: G...')
            # MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
            #         cor.parmdb=self/solutions/cal-g-c'+str(c)+'.h5 cor.correction=amplitudeSmooth cor.invert=False',
            #         log='$nameMS_corrupt.log', commandType='DP3')

            # subtract external region from FR_CORRECTED_DATA (sidelobe subtracted) to create SUBFIELD_DATA
            MSs.addcol('SUBFIELD_DATA','FR_CORRECTED_DATA')
            logger.info('Subtracting external region model (SUBFIELD_DATA = FR_CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = FR_CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE
        
        with w.if_todo('flag_c%02i' % c):
            # Flag on residuals (SUBFIELD_DATA)
            logger.info('Flagging residuals...')
            MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.lua',
                    log='$nameMS_flag-c'+str(c)+'.log', commandType='DP3')
        ### DONE

        with w.if_todo('centralreg_predict_c%02i' % c):
            # Recreate MODEL_DATA of internal region for next calibration cycle
            logger.info('Predict model of internal region...')
            s.add('wsclean -predict -padding 1.8 -name img/wideMint-'+str(c)+' -j '+str(s.max_processors)+' -channels-out '+str(MSs.getChout(4e6))+' '+MSs.getStrWsclean(), \
                   log='wscleanPRE-c'+str(c)+'.log', commandType='wsclean', processors='max')
            s.run(check=True)
        ### DONE

with w.if_todo('final_correct'):
    # correct model with TEC+Beam2ord solutions - ms:FR_CORRECTED_DATA -> ms:CORRECTED_DATA
    #logger.info('Correcting G...')
    #MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=FR_CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
    #        cor.parmdb=self/solutions/cal-g-c{c}.h5 cor.correction=fulljones cor.soltab=[amplitudeSmooth,phase000]',
    #        log='$nameMS_finalcor.log', commandType='DP3')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = FR_CORRECTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    if phaseSolMode in ['tec', 'tecandphase']:
        logger.info('Correcting TEC...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
            cor.parmdb=self/solutions/cal-tec1-c{c}.h5 cor.correction=tec000',
            log='$nameMS_finalcor.log', commandType='DP3')
    if phaseSolMode in ['phase', 'tecandphase']:
        logger.info('Correcting ph...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
            cor.parmdb=self/solutions/cal-tec1-c{c}.h5 cor.correction=phase000',
            log='$nameMS_finalcor.log', commandType='DP3')
    # logger.info('Correct low-res model: TEC+Ph 2...')
    # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
    #         cor.parmdb=self/solutions/cal-tec2-c{c}.h5 cor.correction=tec000',
    #         log='$nameMS_finalcor.log', commandType='DP3')
    # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
    #         cor.parmdb=self/solutions/cal-tec2-c{c}.h5 cor.correction=phase000',
    #         log='$nameMS_finalcor.log', commandType='DP3')
### DONE

with w.if_todo('imaging-final'):
    clean_self(c, MSs, MSsClean, imagenameM, imgsizepix, cc_fit_order, subfield_kwargs, do_predict=False)

# polarisation imaging
with w.if_todo('imaging-pol'):
    logger.info('Cleaning (Pol)...')
    imagenameP = 'img/wideP'
    lib_util.run_wsclean(s, 'wscleanP.log', MSs.getStrWsclean(), name=imagenameP, pol='QUV',
        size=imgsizepix_p, scale='10arcsec', weight='briggs -0.3', niter=0, no_update_model_required='',
        parallel_gridding=2, baseline_averaging='', minuv_l=30, maxuv_l=4500,
        join_channels='', channels_out=MSs.getChout(4.e6))

MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN SUBFIELD_DATA, FR_CORRECTED_DATA"',
        log='$nameMS_taql_delcol.log', commandType='general')

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image*.fits self/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-MFS-residual*.fits self/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources*.txt self/images') for c in range(maxIter) ]
os.system('mv img/wideP-MFS-*-image.fits self/images')
os.system('mv img/wide-lr-MFS-image.fits self/images')
os.system('mv img/wide-largescale-MFS-image.fits self/images')

# Copy model
os.system(f'mv img/wideM-{maxIter-1}-*-model.fits self/skymodel')

logger.info("Done.")
