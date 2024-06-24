#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# they need to be in "./mss/"

import sys, os, glob
import numpy as np
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, h5_merger, lib_h5
logger_obj = lib_log.Logger('pipeline-facet')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-facet.walker')

parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_facet'])))
parset_dir = parset.get('LOFAR_facet','parset_dir')
maxIter = parset.getint('LOFAR_facet','maxIter') # default = 2 (try also 3)
phaseSolMode = parset.get('LOFAR_facet', 'ph_sol_mode') # tecandphase, tec, phase
#sourcedb = parset.get('model','sourcedb')
sourcedb = 'tgts.skymodel'
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
if tint < 4:
    base_solint = int(np.rint(4/tint)) # this is 2 for dutch SPARSE observations
else: base_solint = 1

#################################################################
if not os.path.exists(sourcedb):
    # make filter mask
    imagename_fov = 'img/fov_mask'
    lib_util.run_wsclean(s, 'wsclean_temp.log', MSs.getStrWsclean(), name='img/fov_mask', maxuv_l=1000,
                         size=int(1.2*fov_maj*3600/30), scale='30arcsec')
    os.system('mv img/fov_mask-dirty.fits img/fov_mask.fits')
    lib_img.blank_image_reg(imagename_fov+'.fits', beamReg, blankval=1.)
    lib_img.blank_image_reg(imagename_fov+'.fits', beamReg, blankval=0., inverse=True)

    prefix = '/beegfs/p1uy068/virgo/models/LBA'
    sm = lsmtool.load(f'{prefix}/LVCS_20as_gaul_filtered_freqscaled.skymodel', beamMS=MSs.getListStr()[int(len(MSs.getListStr())/2)])
    sm.select('img/fov_mask.fits=1') # remove distant sources
    sm.select('I>0.4', applyBeam=True) # keep only reasonably bright sources
    sm.write('tgts-pre.skymodel', clobber=True, applyBeam=True, adjustSI=True)
    sm = lsmtool.load('tgts-pre.skymodel')
    sm.group('cluster', numClusters=20, root='MODEL_DATA')
    sm.setPatchPositions(method='wmean')
    sm.plot('self/skymodel/patches.png', 'patch')
    #sm.getColValues('I', aggregate='sum')
    sm.write(sourcedb, clobber=True)
else:
    logger.info(f'Load existing skymodel {sourcedb}')
    sm = lsmtool.load(sourcedb)


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
#####################################################################################################
# Self-cal cycle
for c in range(maxIter):
    logger.info('Start selfcal cycle: '+str(c))
    if c == 0:
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Creating CORRECTED_DATA from FR_CORRECTED_DATA...')
            MSs.addcol('CORRECTED_DATA', 'FR_CORRECTED_DATA')
    else:
        sourcedb = f'tgts-c{c}.skymodel'
        # with w.if_todo('set_corrected_data_c%02i' % c):
        #     logger.info('Set CORRECTED_DATA = SUBFIELD_DATA...')
        #     MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBFIELD_DATA"', log='$nameMS_taql-c' + str(c) + '.log',
        #             commandType='general')

        if not os.path.exists(sourcedb):
            # make filter mask
            sm = lsmtool.load(f'img/wideM-{c-1}-sources.txt')
            sm.select('img/fov_mask.fits=1')  # remove distant sources
            sm.select('I>0.05')  # remove the faintest sources
            sm.group('tessellate', targetFlux=15, root='MODEL_DATA')
            sm.setPatchPositions(method='wmean')
            sm.plot(f'self/skymodel/patches-c{c}.png', 'patch')
            # sm.getColValues('I', aggregate='sum')
            sm.write(sourcedb, clobber=True)
        else:
            logger.info(f'Load existing skymodel {sourcedb}')
            sm = lsmtool.load(sourcedb)
        ### DONE

    # with w.if_todo('smooth_model_c%02i' % c):
    #     # Smooth MODEL_DATA -> MODEL_DATA
    #     MSs.run_Blsmooth('MODEL_DATA', 'MODEL_DATA', logstr=f'smooth-c{c}')
    ### DONE

    with w.if_todo('sol_avg_%02i' % c):
        MSs = lib_ms.AllMSs(glob.glob('mss/TC*[0-9].MS'), s)
        lib_util.check_rm('mss-sol')
        os.system('mkdir mss-sol')
        timeint = MSs.getListObj()[0].getTimeInt()
        avgtimeint = int(round(8 / timeint))
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
    # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
    sourcedb_basename = sourcedb.split('/')[-1]
    for MS in MSs_sol.getListStr():
        lib_util.check_rm(MS + '/' + sourcedb_basename)
        logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
        os.system('cp -r ' + sourcedb + ' ' + MS)
    # Add model to MODEL_DATA
    skymodel = lsmtool.load(sourcedb)
    patches = skymodel.getPatchNames()
    with w.if_todo(f'init_model_c{c}'):
        for patch in patches:
            logger.info(f'Add model to {patch}...')
            MSs_sol.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.usebeammodel=false pre.sourcedb=$pathMS/{sourcedb_basename} pre.sources={patch} msout.datacolumn={patch}',
            log='$nameMS_pre.log', commandType='DP3')
            # Smooth CORRECTED_DATA -> SMOOTHED_DATA
            # MSs_sol.run_Blsmooth(patch, patch, logstr=f'smooth-c{c}')
        ### DONE

    with w.if_todo('solve_tecRS_c%02i' % c):
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        MSs_sol.run_Blsmooth('DATA', logstr=f'smooth-c{c}')
        # solve ionosphere phase - ms:SMOOTHED_DATA
        logger.info('Solving TEC (RS)...')
        if phaseSolMode == 'phase': #phase
            solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=1.0e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
            losoto_parsets = [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-resetph2-CSRS.parset', parset_dir+'/losoto-plot-scalar.parset']
        else: # TEC or TecAndPhase
            solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
            losoto_parsets = [parset_dir+'/losoto-plot-tec.parset']

        MSs_sol.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec2.h5 sol.solint=1 {solver_params} sol.modeldatacolumns="[{",".join(patches)}]"',
                    log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec2-c'+str(c), [ms+'/tec2.h5' for ms in MSs_sol.getListStr()], losoto_parsets)
        os.system('mv cal-tec2-c'+str(c)+'.h5 self/solutions/')
        os.system('mv plots-tec2-c'+str(c)+' self/plots/')

    with w.if_todo('corrupt_tecRS_c%02i' % c):
        ### CORRUPT the MODEL_DATA columns for all patches
        logger.info(f'Corrupt models: {phaseSolMode}2...')
        for patch in patches:
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs_sol.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={patch} msout.datacolumn={patch} cor.direction=[{patch}]  \
                        cor.parmdb=self/solutions/cal-tec2-c' + str(c) + '.h5 cor.correction=tec000 cor.invert=False',
                            log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs_sol.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={patch} msout.datacolumn={patch} cor.direction=[{patch}] \
                        cor.parmdb=self/solutions/cal-tec2-c' + str(c) + '.h5 cor.correction=phase000 cor.invert=False',
                            log='$nameMS_corrupt.log', commandType='DP3')

    with w.if_todo('solve_tecCS_c%02i' % c):

        # solve ionosphere phase - ms:SMOOTHED_DATA
        logger.info('Solving TEC (CS)...')
        if phaseSolMode == 'phase': #phase
            solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=3.0e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
            losoto_parsets = [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-resetph2-CS.parset', parset_dir+'/losoto-plot-scalar.parset']
        else: # TEC or TecAndPhase
            solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
            losoto_parsets = [parset_dir+'/losoto-plot-tec.parset']

        MSs_sol.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec1.h5 sol.solint=8 {solver_params} sol.modeldatacolumns="[{",".join(patches)}]"',
                    log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec1-c'+str(c), [ms+'/tec1.h5' for ms in MSs_sol.getListStr()], losoto_parsets)
        os.system('mv cal-tec1-c'+str(c)+'.h5 self/solutions/')
        os.system('mv plots-tec1-c'+str(c)+' self/plots/')
    ### DONE

    with w.if_todo('corrupt_tecCS_c%02i' % c):
        ### CORRUPT the MODEL_DATA columns for all patches
        logger.info(f'Corrupt models: {phaseSolMode}1...')
        for patch in patches:
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs_sol.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={patch} msout.datacolumn={patch} cor.direction=[{patch}]  \
                        cor.parmdb=self/solutions/cal-tec1-c' + str(c) + '.h5 cor.correction=tec000 cor.invert=False',
                            log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs_sol.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={patch} msout.datacolumn={patch} cor.direction=[{patch}] \
                        cor.parmdb=self/solutions/cal-tec1-c' + str(c) + '.h5 cor.correction=phase000 cor.invert=False',
                            log='$nameMS_corrupt.log', commandType='DP3')

    with w.if_todo('solve_tecCS0_c%02i' % c):

        # solve ionosphere phase - ms:SMOOTHED_DATA
        logger.info('Solving TEC (slowCS)...')
        if phaseSolMode == 'phase': #phase
            solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=6.0e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
            losoto_parsets = [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-plot-scalar.parset']
        else: # TEC or TecAndPhase
            solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
            losoto_parsets = [parset_dir+'/losoto-plot-tec.parset']

        MSs_sol.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec0.h5 sol.solint=64 {solver_params} sol.modeldatacolumns="[{",".join(patches)}]"',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec0-c'+str(c), [ms+'/tec0.h5' for ms in MSs_sol.getListStr()], losoto_parsets, plots_dir=f'self/plots/plots-tec0-c{c}', h5_dir=f'self/solutions/')
    ### DONE

   # need: sol1, sol2, sm,
    with w.if_todo(f'merge_h5_c{c}'):
        sol_dir = 'self/solutions'
        lib_util.check_rm(f'{sol_dir}/cal-tec-merged-c{c}.h5')
        # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec0-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec1-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec2-c{c}.h5', sourcedb)
        # reference, unflag and reset the added stations
        # upsample the CS h5parm and merge it with the RS h5parm
        # TODO mss-sol
        h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-tec-merged-c{c}.h5", min_distance=10/3600,
                           h5_tables=[f'{sol_dir}/cal-tec0-c{c}.h5',f'{sol_dir}/cal-tec1-c{c}.h5',f'{sol_dir}/cal-tec2-c{c}.h5'],
                           h5_time_freq=f'{sol_dir}/cal-tec2-c{c}.h5', no_pol=True, ms_files='mss/TC*.MS',)
        lib_util.run_losoto(s, f'tec-merged-c{c}', f'{sol_dir}/cal-tec-merged-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
                            plots_dir=f'self/plots/plots-tec-merged-c{c}', h5_dir='self/solutions')

    facetregname = 'self/solutions/facets.reg'
    # TODO temp CORRECTED_DATA
    MSs = MSs_sol

    with w.if_todo('c%02i-imaging' % c):
        logger.info('Preparing region file...')
        s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+f' --h5 self/solutions/cal-tec-merged-c{c}.h5 --imsize '+str(imgsizepix_wide)+' \
                --pixelscale 4 --writevoronoipoints --output '+facetregname,
              log='facet_generator.log', commandType='python')
        s.run()

        imagename = 'img/wide-' + str(c)
        imagenameM = 'img/wideM-' + str(c)
        # make quick image with -local-rms to get a mask
        #  dd_psf_grid='25 25', beam_size=15,CORRECTED_ dd_psf_grid='12 12', apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='',
        logger.info('Cleaning 1...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, data_column='DATA', size=imgsizepix_wide, scale='4arcsec', concat_mss=True, keep_concat=True,
                             weight='briggs -0.3', niter=1000000, gridder='wgridder', parallel_gridding=6, no_update_model_required='', minuv_l=30, mgain=0.85, parallel_deconvolution=1024,
                             auto_threshold=5.0, auto_mask=8.0, join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='', pol='i', nmiter=6,  facet_regions=facetregname, apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000' )
        # # masking
        if userReg:
            logger.info(f'Making mask with userReg {userReg}...')
            s.add('breizorro.py -t 6.5 -r %s -b 50 -o %s --merge %s' % (
                imagename + '-MFS-image.fits', imagename + '-mask.fits', userReg),
                  log=f'makemask-{c}.log', commandType='python')
            s.run()
        else:
            logger.info('Making mask...')
            s.add('breizorro.py -t 6.5 -r %s -b 50 -o %s' % (imagename + '-MFS-image.fits', imagename + '-mask.fits'),
                  log=f'makemask-{c}.log', commandType='python')
            s.run()

        # clean again, with mask nowe,
        logger.info('Cleaning 2...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM,  fits_mask=imagename+'-mask.fits', data_column='DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder',  parallel_gridding=6, save_source_list='', concat_mss=True,# reuse_concat=True,
                             no_update_model_required='', minuv_l=30, mgain=0.8, parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=5.0,
                             join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='', pol='i', facet_regions=facetregname, apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000' )

sys.exit()



if False: # add only later
        # full_image.nantozeroModel()
        # s.add('wsclean -predict -padding 1.8 -name ' + full_image.root + ' -j ' + str(
        #     s.max_processors) + ' -channels-out ' + str(ch_out) + ' \
        #         -apply-facet-beam -use-differential-lofar-beam -facet-beam-update 120 \
        #         -facet-regions ' + facetregname + ' -diagonal-solutions \
        #         -apply-facet-solutions ' + interp_h5parm_old + ' ' + correct_for + ' \
        #         -reorder -parallel-reordering 4 ' + MSs.getStrWsclean(),
        #       log='wscleanPRE-c' + str(cmaj) + '.log', commandType='wsclean', processors='max')
        # s.run(check=True)
        # TODO
        # is the model at this point corrupted for phases?
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

            logger.info('Preparing region file...')
            facetregnameLR = 'self/solutions/facetsLR.reg'

            s.add('ds9_facet_generator.py --ms ' + MSs.getListStr()[
                0] + f' --h5 self/solutions/cal-tec-merged-c{c}.h5 --imsize ' + str(imgsizepix_lr) + ' \
                    --pixelscale 4 --writevoronoipoints --output ' + facetregnameLR,
                  log='facet_generator.log', commandType='python')
            s.run()
            # reclean low-resolution
            logger.info('Cleaning low-res...')
            imagename_lr = 'img/wide-lr'
            lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=True,
                                 parallel_gridding=4, temp_dir='../', size=imgsizepix_lr, scale='30arcsec',
                                 weight='briggs -0.3', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
                                 taper_gaussian='200arcsec', mgain=0.85, parallel_deconvolution=512, baseline_averaging='',
                                 local_rms='', auto_mask=3, auto_threshold=1.5, fits_mask='img/wide-lr-mask.fits',
                                 join_channels='', channels_out=MSs.getChout(2.e6), facet_regions=facetregnameLR,
                                 apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000')


        with w.if_todo('lowres_subtract_c%02i' % c):
            # Permanently subtract low-res sidelobe model - FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA.
            logger.info('Subtracting low-res sidelobe model (FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA"',
                    log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
            # For flagging: subtract low-res sidelobe model - CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA.
            logger.info('Subtracting low-res sidelobe model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"',
                    log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
        ### DONE

        with w.if_todo('flag_c%02i' % c):
            # Flag on residuals (SUBFIELD_DATA)
            logger.info('Flagging residuals...')
            MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.lua',
                    log='$nameMS_flag-c'+str(c)+'.log', commandType='DP3')
        ### DONE

        # with w.if_todo('centralreg_predict_c%02i' % c):
        #     # Recreate MODEL_DATA of internal region for next calibration cycle
        #     logger.info('Predict back model...')
        #     s.add(f'wsclean -predict -padding 1.8 -name {imagenameM} -j {s.max_processors} -channels-out {MSs.getChout(4e6)} \
        #            -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged.h5 phase000 \
        #            {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
        #     s.run(check=True)
        # ### DONE

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

# Copy model
os.system(f'mv img/wideM-{maxIter-1}-*-model.fits self/skymodel')

logger.info("Done.")
