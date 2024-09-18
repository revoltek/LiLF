#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Perform self-calibration on a group of SBs concatenated in TCs.
# The algorithm works as follows:
# 1. Solve & apply direction-independent FR correction from LL-RR
# 2. Solve in a small number of directions (~5) against GSM or HBA input model. Solving happens in an iterative way,
#    using longer time intervals for the central stations.
# 3. Image the sky using wsclean facet-mode.
# 4. Find a subfield containing at least subfield_min_flux Jy or flux (default=40Jy). Subtract all but this subfield
#    using the dd-solutions
# 5. Solve for the sub-field and apply the solutions to the data
# 6. Repeat dd self-cal cycles with a growing number of directions
# they need to be in "./mss/"

# TODO figure out if it is better to properly predict with the beam
# TODO add LoTSS query for statring model, only if not available use GSM
# TODO different solints for different brightness facets?
# TODO implement sidelob subtraction + peeling of bright sidelobe sourced
# TODO do we need amplitude corrections(?)? If so, DI/DD/only for bright sources?
# TODO flag on empty data
# TODO add PBcorr + final imaging products

import sys, os, glob
import numpy as np
from regions import Regions
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, h5_merger, lib_h5, lib_dd_parallel
logger_obj = lib_log.Logger('pipeline-self')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-self.walker')

parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_self'])))
parset_dir = parset.get('LOFAR_parallelself','parset_dir')
subfield_min_flux = parset.getfloat('LOFAR_self','subfield_min_flux') # default 40 Jy
subfield = parset.get('LOFAR_self','subfield') # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
maxIter = parset.getint('LOFAR_self','maxIter') # default = 2 (try also 3)
phaseSolMode = parset.get('LOFAR_self', 'ph_sol_mode') # tecandphase, tec, phase
#sourcedb = parset.get('model','sourcedb')
sourcedb = 'tgts.skymodel'
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')

#############################################################################

def clean_empty(MSs, name, col='CORRECTED_DATA'):
    """ For testing/debugging only"""
    lib_util.run_wsclean(s, 'wsclean-empty.log', MSs.getStrWsclean(), name=name,
                         data_column=col, size=5000, scale='8arcsec', niter=0, nmiter=0,
                         weight='briggs 0.0', gridder='wgridder', parallel_gridding=6,
                         no_update_model_required='')

def corrupt_model_dirs(c, tc, patches):
    """ CORRUPT the MODEL_DATA columns for all patches
    Parameters
    ----------
    c: int, cycle
    tc: tec/phase step (e.g. 2 -
    patches: list of patch names
    """
    logger.info(f'Corrupt models: {phaseSolMode}{tc}...')
    for patch in patches:
        if phaseSolMode in ['tec', 'tecandphase']:
            MSs_sol.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={patch} msout.datacolumn={patch} cor.direction=[{patch}]  \
                    cor.parmdb=self/solutions/cal-tec{tc}-c{c}.h5 cor.correction=tec000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        if phaseSolMode in ['phase', 'tecandphase']:
            MSs_sol.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={patch} msout.datacolumn={patch} cor.direction=[{patch}] \
                    cor.parmdb=self/solutions/cal-tec{tc}-c{c}.h5 cor.correction=phase000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')

def solve_iono(c, tc, patches, smMHz, solint, resetant_parset=None):
    """
    Parallel solve for ionosphere systematics
    Parameters
    ----------
    c: major cycle
    tc: tec/phase solve step
    patches: list of patch names
    smMHz: smoothnessconstraint kernel size at 54 MHz in MHz
    solint: solution interval in timesteps (e.g. multiples of 8s)
    resetant_parset: string, optional. Extra losoto parset to reset solutions for some stations
    """

    if phaseSolMode == 'phase': #phase
        solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint={smMHz}e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
        if resetant_parset is not None:
            losoto_parsets = [parset_dir+'/losoto-refph.parset', resetant_parset, parset_dir+'/losoto-plot-scalar.parset']
        else:
            losoto_parsets = [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-plot-scalar.parset']
    else: # TEC or TecAndPhase
        solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
        losoto_parsets = [parset_dir+'/losoto-plot-tec.parset']
        if resetant_parset is not None:
            raise NotImplementedError('Resetant for TEC not yet implemented.')

    MSs_sol.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec{tc}.h5 sol.solint={solint} {solver_params} sol.modeldatacolumns="[{",".join(patches)}]"',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

    lib_util.run_losoto(s, f'tec{tc}-c{c}', [ms+f'/tec{tc}.h5' for ms in MSs_sol.getListStr()], losoto_parsets, plots_dir=f'self/plots/plots-tec{tc}-c{c}', h5_dir=f'self/solutions/')

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
# define list of facet fluxes per iteration -> this can go into the config
facet_fluxes = [15, 6, 3.5, 2.0, 1.5] # still needs to be tuned, maybe also depends on the field
#################################################################
# Create initial sourcedb
if not os.path.exists(sourcedb):
    fwhm = MSs.getListObj()[0].getFWHM(freq='min')
    radeg = phasecentre[0]
    decdeg = phasecentre[1]
    os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm))
    sm = lsmtool.load(f'tgts.skymodel', beamMS=MSs.getListStr()[int(len(MSs.getListStr())/2)])
    sm.select('I>0.05', applyBeam=True) # keep only reasonably bright sources
    sm.write('tgts-pre.skymodel', clobber=True, applyBeam=True, adjustSI=True)
    sm = lsmtool.load('tgts-pre.skymodel')
    sm.group('tessellate', targetFlux=facet_fluxes[0], root='MODEL_DATA') # once we have LoTSS query implemented, we can probaly start with a higher number of patches here?
    sm.setPatchPositions(method='wmean')
    sm = lib_dd_parallel.merge_faint_facets(sm, 0.5*facet_fluxes[0])  # merge all facets below threshold
    sm.plot('self/skymodel/patches.png', 'patch')
    sm.write(sourcedb, clobber=True)
else:
    logger.info(f'Load existing skymodel {sourcedb}')
    sm = lsmtool.load(sourcedb)

#################################################################################################

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

    logger.info('Creating MODEL_DATA_FR...')  # take from MODEL_DATA but overwrite
    MSs.addcol('MODEL_DATA_FR', 'DATA', usedysco=False)
    MSs.run('taql "UPDATE $pathMS SET MODEL_DATA_FR[,0]=0.5+0i, MODEL_DATA_FR[,1]=0.0+0i, MODEL_DATA_FR[,2]=0.0+0i, \
         MODEL_DATA_FR[,3]=0.5+0i"', log='$nameMS_taql_frmodel.log', commandType='general')

    # Solve cal_SB.MS:CIRC_PHASEDIFF_DATA against MODEL_DATA_FR (only solve - solint=2m nchan=0 as it has the smoothnessconstrain)
    logger.info('Solving circ phase difference ...')
    MSs.run('DP3 ' + parset_dir + '/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.solint=' + str(30 * base_solint),
            log='$nameMS_solFR.log', commandType="DP3")
    lib_util.run_losoto(s, f'fr', [ms + '/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
    os.system('mv cal-fr.h5 self/solutions/')
    os.system('mv plots-fr self/plots/')
    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, MODEL_DATA_FR"',
            log='$nameMS_taql_delcol.log', commandType='general')

    # Correct FR with results of solve - group*_TC.MS:DATA -> group*_TC.MS:FR_CORRECTED_DATA
    logger.info('Correcting FR...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=FR_CORRECTED_DATA \
                cor.parmdb=self/solutions/cal-fr.h5 cor.correction=rotationmeasure000',
            log='$nameMS_corFR.log', commandType='DP3')
### DONE

with w.if_todo('template_h5'): # create a template h5parm file
    if phaseSolMode == 'phase':  # phase
        solver_params = f'sol.mode=scalarphase sol.nchan=1'
    else:  # TEC or TecAndPhase
        solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'

    MSs.run(
        f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=FR_CORRECTED_DATA sol.h5parm=$pathMS/template.h5 sol.solint=1 sol.maxiter=1 {solver_params} sol.modeldatacolumns="[DATA]"',
        log='$nameMS_soltemplate.log', commandType='DP3')
    lib_util.run_losoto(s, 'template', [ms + '/template.h5' for ms in MSs.getListStr()], [])
    os.system('mv cal-template.h5 self/solutions/')
### DONE

#####################################################################################################
# Self-cal cycle
for c in range(maxIter):
    logger.info('Start selfcal cycle: '+str(c))
    if c == 0:
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Creating CORRECTED_DATA from FR_CORRECTED_DATA...')
            MSs.addcol('CORRECTED_DATA', 'FR_CORRECTED_DATA')
    elif c == 1 :
        # Only in first iteration: apply the subfield solutions to the data
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Set CORRECTED_DATA = FR_CORRECTED_DATA...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = FR_CORRECTED_DATA"', log='$nameMS_taql-c' + str(c) + '.log',
                    commandType='general')
            logger.info('Correct subfield iono...')
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        cor.parmdb=self/solutions/cal-sf-tec2-c' + str(c-1) + '.h5 cor.correction=tec000 ',
                    log='$nameMS_corrupt.log', commandType='DP3')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        cor.parmdb=self/solutions/cal-sf-tec1-c' + str(c-1) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_corrupt.log', commandType='DP3')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        cor.parmdb=self/solutions/cal-sf-tec0-c' + str(c-1) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA  \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(c-1) + '.h5 cor.correction=phase000 \
                        cor.invert=True',log='$nameMS_corrupt.log', commandType='DP3')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA  \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(c-1) + '.h5 cor.correction=phase000 \
                        cor.invert=True',log='$nameMS_corrupt.log', commandType='DP3')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA  \
                        cor.parmdb=self/solutions/cal-tec0-sf-c' + str(c-1) + '.h5 cor.correction=phase000 \
                        cor.invert=True',log='$nameMS_corrupt.log', commandType='DP3')
    if c > 0:
        sourcedb = f'tgts-c{c}.skymodel'
        if not os.path.exists(sourcedb):
            sm = lsmtool.load(f'img/wideM-{c-1}-sources.txt')
            sm.select('I>0.01')  # remove the faintest sources
            sm.group('tessellate', targetFlux=facet_fluxes[c], root='MODEL_DATA')
            sm.setPatchPositions(method='wmean')
            sm = lib_dd_parallel.merge_faint_facets(sm, 0.5*facet_fluxes[c]) # merge all facets below targetFlux
            sm.plot(f'self/skymodel/patches-c{c}.png', 'patch')
            sm.write(sourcedb, clobber=True)
        else:
            logger.info(f'Load existing skymodel {sourcedb}')
            sm = lsmtool.load(sourcedb)

    # with w.if_todo('smooth_model_c%02i' % c):
    #     # Smooth MODEL_DATA -> MODEL_DATA
    #     MSs.run_Blsmooth('MODEL_DATA', 'MODEL_DATA', logstr=f'smooth-c{c}')
    ### DONE

    # Create an averaged data set, only used for the DP3 solve. This speeds up the solve and also it strongly reduces the
    # disk footprint (each facet get's its own MODEL_DATA!)
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
        # How many channels for solve
        nchan = nchan_init - nchan_init % 48
        freqstep = 2 if c == 0 else 4 # can average more after DI-correction is applied!
        logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (
            timeint, timeint * avgtimeint, nchan_init, nchan/freqstep))
        MSs.run(
            'DP3 ' + parset_dir + '/DP3-avg.parset msin=$pathMS msout=mss-sol/$nameMS.MS msin.datacolumn=CORRECTED_DATA msin.nchan=' + str(
                nchan) + ' \
                avg.timestep=' + str(avgtimeint) + f' avg.freqstep={freqstep}',
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
        smMHz = 0.8 if c == 0 else 1.6
        # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for all BUT most distant RS!
        logger.info('Solving TEC (fastRS)...')
        solve_iono(c, 2, patches, smMHz, 1, resetant_parset=parset_dir+'/losoto-resetph2-CSRS.parset')
    ### DONE

    with w.if_todo('corrupt_tecRS_c%02i' % c):
        ### CORRUPT the MODEL_DATA columns for all patches
        corrupt_model_dirs(c, 2, patches)
    ### DONE

    with w.if_todo('solve_tecCS_c%02i' % c):
        smMHz = 3.0 if c == 0 else 5.0
        # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for central CS
        logger.info('Solving TEC (midRS)...')
        solve_iono(c, 1, patches, smMHz, 4, resetant_parset=parset_dir+'/losoto-resetph2-CS.parset')
    ### DONE

    with w.if_todo('corrupt_tecCS_c%02i' % c):
        ### CORRUPT the MODEL_DATA columns for all patches
        corrupt_model_dirs(c, 1, patches)
    ### DONE

    with w.if_todo('solve_tecCS0_c%02i' % c):
        # solve ionosphere phase - ms:SMOOTHED_DATA
        smMHz = 6.0 if c == 0 else 8.0
        logger.info('Solving TEC (slowCS)...')
        solve_iono(c, 0, patches, smMHz, 16)
    ### DONE

    # merge solutions into one h5parms for large scale image
    with w.if_todo(f'merge_h5_c{c}'):
        sol_dir = 'self/solutions'
        lib_util.check_rm(f'{sol_dir}/cal-tec-merged-c{c}.h5')
        # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec0-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec1-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec2-c{c}.h5', sourcedb)
        # reference, unflag and reset the added stations
        h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-tec-merged-c{c}.h5", min_distance=10/3600,
                           h5_tables=[f'{sol_dir}/cal-tec0-c{c}.h5',f'{sol_dir}/cal-tec1-c{c}.h5',f'{sol_dir}/cal-tec2-c{c}.h5'],
                           h5_time_freq=f'{sol_dir}/cal-template.h5', no_pol=True, ms_files='mss/TC*.MS',no_antenna_crash=True)
        lib_util.run_losoto(s, f'tec-merged-c{c}', f'{sol_dir}/cal-tec-merged-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
                            plots_dir=f'self/plots/plots-tec-merged-c{c}', h5_dir='self/solutions')
    facetregname = f'self/solutions/facets-c{c}.reg'

    with w.if_todo('c%02i-imaging' % c):
        logger.info('Preparing region file...')
        s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+f' --h5 self/solutions/cal-tec-merged-c{c}.h5 --imsize '+str(imgsizepix_wide)+' \
                --pixelscale 4 --writevoronoipoints --output '+facetregname,
              log='facet_generator.log', commandType='python')
        s.run()

        imagename = 'img/wide-' + str(c)
        imagenameM = 'img/wideM-' + str(c)
        # make quick image with -local-rms to get a mask
        # TODO make this faster - experiment with increased parallel-gridding as well as shared facet reads option
        logger.info('Cleaning 1...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, data_column='CORRECTED_DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', local_rms='', niter=1000000, gridder='wgridder',  parallel_gridding=6, no_update_model_required='', minuv_l=30, mgain=0.85, parallel_deconvolution=1024,
                             auto_threshold=5.0, auto_mask=8.0, beam_size=15, join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='', pol='i', nmiter=6,   facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000' )
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
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM,  fits_mask=imagename+'-mask.fits', reuse_psf=imagename, reuse_dirty=imagename, data_column='CORRECTED_DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder',  parallel_gridding=6, save_source_list='',
                             no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=5.0,
                             join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='',  multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80',   pol='i', facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000',
                             apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='')

        logger.info('Predict corrupted MODEL_DATA...')
        s.add(f'wsclean -predict -padding 1.8 -name {imagenameM} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
               -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged-c{c}.h5 phase000 \
                -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam {MSs.getStrWsclean()}',
               log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
        s.run(check=True)

    if c == 0:
        # TODO: we could do a second round of calibration here before calibrating the sub-field to get a more accurate subtraction
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
                # sm.remove('img/wide-lr-mask.fits=1')  # remove sidelobe sources that were subtracted
                sm.remove('MajorAxis > 80')  # remove largest scales
                field_center1, field_size1 = lib_dd.make_subfield_region(subfield_path, MSs.getListObj()[0], sm,
                                                                         subfield_min_flux, debug_dir='img/')
            # prepare model of central/external regions
            logger.info('Blanking central region of model files and reverse...')
            for im in glob.glob('img/wideM-0*model*.fits'):
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

        with w.if_todo('extreg_subtract_c%02i' % c):
            # Recreate MODEL_DATA of external region for subtraction
            logger.info('Predict corrupted model of external region...')
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMext-{c} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                   -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam -facet-regions {facetregname} \
                   -apply-facet-solutions self/solutions/cal-tec-merged-c{c}.h5 phase000 {MSs.getStrWsclean()}',
                log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
            s.run(check=True)

            # subtract external region from FR_CORRECTED_DATA to create SUBFIELD_DATA
            MSs.addcol('SUBFIELD_DATA','FR_CORRECTED_DATA')
            logger.info('Subtracting external region model (SUBFIELD_DATA = FR_CORRECTED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = FR_CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
            clean_empty(MSs,'only_subfield', 'SUBFIELD_DATA')
        ### DONE

        with w.if_todo('intreg_predict%02i' % c):
            # Recreate MODEL_DATA of internal region for solve
            logger.info('Predict model of internal region...')
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMint-{c} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                   -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam {MSs.getStrWsclean()}',
                  log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
            s.run(check=True)

        with w.if_todo('subfield_solve_tecRS_c%02i' % c):
            # Smooth CORRECTED_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving TEC (RS)...')
            if phaseSolMode == 'phase':  # phase
                solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=0.7e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
                losoto_parsets = [parset_dir + '/losoto-refph.parset', parset_dir + '/losoto-resetph2-CSRS.parset',
                                  parset_dir + '/losoto-plot-scalar.parset']
            else:  # TEC or TecAndPhase
                solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
                losoto_parsets = [parset_dir + '/losoto-plot-tec.parset']

            MSs.run(
                f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec2.h5 sol.solint=1 {solver_params} sol.modeldatacolumns="[MODEL_DATA]"',
                log='$nameMS_solTEC-c' + str(c) + '.log', commandType='DP3')

            lib_util.run_losoto(s, 'tec2-sf-c' + str(c), [ms + '/tec2.h5' for ms in MSs.getListStr()], losoto_parsets)
            os.system('mv cal-tec2-sf-c' + str(c) + '.h5 self/solutions/')
            os.system('mv plots-tec2-sf-c' + str(c) + ' self/plots/')

        with w.if_todo('subfield_corrupt_tecRS_c%02i' % c):
            ### CORRUPT the MODEL_DATA columns for all patches
            logger.info(f'Corrupt MODEL_DATA...')
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(
                    f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(
                        c) + '.h5 cor.correction=tec000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(
                    f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(
                        c) + '.h5 cor.correction=phase000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')

        with w.if_todo('subfield_solve_tecCS_c%02i' % c):
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving TEC (CS)...')
            if phaseSolMode == 'phase':  # phase
                solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=3.0e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
                losoto_parsets = [parset_dir + '/losoto-refph.parset', parset_dir + '/losoto-resetph2-CS.parset',
                                  parset_dir + '/losoto-plot-scalar.parset']
            else:  # TEC or TecAndPhase
                solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
                losoto_parsets = [parset_dir + '/losoto-plot-tec.parset']

            MSs.run(
                f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec1.h5 sol.solint=8 {solver_params} sol.modeldatacolumns="[MODEL_DATA]"',
                log='$nameMS_solTEC-c' + str(c) + '.log', commandType='DP3')

            lib_util.run_losoto(s, 'tec1-sf-c' + str(c), [ms + '/tec1.h5' for ms in MSs.getListStr()], losoto_parsets)
            os.system('mv cal-tec1-sf-c' + str(c) + '.h5 self/solutions/')
            os.system('mv plots-tec1-sf-c' + str(c) + ' self/plots/')
        ### DONE

        with w.if_todo('subfield_corrupt_tecCS_c%02i' % c):
            ### CORRUPT the MODEL_DATA columns for all patches
            logger.info(f'Corrupt MODEL_DATA: {phaseSolMode}1...')
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(
                    f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(
                        c) + '.h5 cor.correction=tec000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(
                    f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(
                        c) + '.h5 cor.correction=phase000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')

        with w.if_todo('subfield_solve_tecCS0_c%02i' % c):
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving TEC (slowCS)...')
            if phaseSolMode == 'phase':  # phase
                solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=6.0e6 sol.smoothnessreffrequency=54e6 sol.nchan=1'
                losoto_parsets = [parset_dir + '/losoto-refph.parset', parset_dir + '/losoto-plot-scalar.parset']
            else:  # TEC or TecAndPhase
                solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
                losoto_parsets = [parset_dir + '/losoto-plot-tec.parset']

            MSs.run(
                f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec0.h5 sol.solint=16 {solver_params} sol.modeldatacolumns="[MODEL_DATA]"',
                log='$nameMS_solTEC-c' + str(c) + '.log', commandType='DP3')

            lib_util.run_losoto(s, 'tec0-sf-c' + str(c), [ms + '/tec0.h5' for ms in MSs.getListStr()], losoto_parsets,
                                plots_dir=f'self/plots/plots-tec0-sf-c{c}', h5_dir=f'self/solutions/')
        ### DONE
        #
        # with w.if_todo('lowres_setdata_c%02i' % c):
        #     # Subtract model from all TCs - ms:CORRECTED_DATA - MODEL_DATA -> ms:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        #     logger.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        #     MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        # ### DONE
        #
        # with w.if_todo('lowres_img_c%02i' % c):
        #     # Making beam mask
        #     logger.info('Preparing mask for low-res clean...')
        #     lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name='img/tmp', size=imgsizepix_lr, scale='30arcsec')
        #     os.system('mv img/tmp-image.fits img/wide-lr-mask.fits')
        #     lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 0.)
        #     lib_img.blank_image_reg('img/wide-lr-mask.fits', beamReg, blankval = 1., inverse=True)
        #
        #     # logger.info('Preparing region file...')
        #     # facetregnameLR = 'self/solutions/facetsLR.reg'
        #     #
        #     # s.add('ds9_facet_generator.py --ms ' + MSs.getListStr()[
        #     #     0] + f' --h5 self/solutions/cal-tec-merged-c{c}.h5 --imsize ' + str(imgsizepix_lr) + ' \
        #     #         --pixelscale 4 --writevoronoipoints --output ' + facetregnameLR,
        #     #       log='facet_generator.log', commandType='python')
        #     # s.run()
        #     # reclean low-resolution
        #     logger.info('Cleaning low-res...')
        #     imagename_lr = 'img/wide-lr'
        #     lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=True,
        #                          parallel_gridding=4, temp_dir='../', size=imgsizepix_lr, scale='30arcsec',
        #                          weight='briggs -0.3', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
        #                          taper_gaussian='200arcsec', mgain=0.85, parallel_deconvolution=512, baseline_averaging='',
        #                          local_rms='', auto_mask=3, auto_threshold=1.5, fits_mask='img/wide-lr-mask.fits',
        #                          join_channels='', channels_out=MSs.getChout(2.e6))
        #
        # with w.if_todo('lowres_subtract_c%02i' % c):
        #     # Permanently subtract low-res sidelobe model - FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA.
        #     logger.info('Subtracting low-res sidelobe model (FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA)...')
        #     MSs.run('taql "update $pathMS set FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA"',
        #             log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        #     # For flagging: subtract low-res sidelobe model - CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA.
        #     logger.info('Subtracting low-res sidelobe model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        #     MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"',
        #             log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
        # ### DONE
        #
        # with w.if_todo('flag_c%02i' % c):
        #     # Flag on residuals (SUBFIELD_DATA)
        #     logger.info('Flagging residuals...')
        #     MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.lua',
        #             log='$nameMS_flag-c'+str(c)+'.log', commandType='DP3')
        # ### DONE
        # sys.exit()
        # TODO from here on, make DI-solve for DI-phases and apply.
        # Then, again Facet-cal + imaging for final subtract model.

        # with w.if_todo('centralreg_predict_c%02i' % c):
        #     # Recreate MODEL_DATA of internal region for next calibration cycle
        #     logger.info('Predict back model...')
        #     s.add(f'wsclean -predict -padding 1.8 -name {imagenameM} -j {s.max_processors} -channels-out {MSs.getChout(4e6)} \
        #            -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged.h5 phase000 \
        #            {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
        #     s.run(check=True)
        # ### DONE

# with w.if_todo('final_correct'):
#     # correct model with TEC+Beam2ord solutions - ms:FR_CORRECTED_DATA -> ms:CORRECTED_DATA
#     #logger.info('Correcting G...')
#     #MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=FR_CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
#     #        cor.parmdb=self/solutions/cal-g-c{c}.h5 cor.correction=fulljones cor.soltab=[amplitudeSmooth,phase000]',
#     #        log='$nameMS_finalcor.log', commandType='DP3')
#     MSs.run('taql "update $pathMS set CORRECTED_DATA = FR_CORRECTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
#     if phaseSolMode in ['tec', 'tecandphase']:
#         logger.info('Correcting TEC...')
#         MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
#             cor.parmdb=self/solutions/cal-tec1-c{c}.h5 cor.correction=tec000',
#             log='$nameMS_finalcor.log', commandType='DP3')
#     if phaseSolMode in ['phase', 'tecandphase']:
#         logger.info('Correcting ph...')
#         MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
#             cor.parmdb=self/solutions/cal-tec1-c{c}.h5 cor.correction=phase000',
#             log='$nameMS_finalcor.log', commandType='DP3')
    # logger.info('Correct low-res model: TEC+Ph 2...')
    # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
    #         cor.parmdb=self/solutions/cal-tec2-c{c}.h5 cor.correction=tec000',
    #         log='$nameMS_finalcor.log', commandType='DP3')
    # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA  \
    #         cor.parmdb=self/solutions/cal-tec2-c{c}.h5 cor.correction=phase000',
    #         log='$nameMS_finalcor.log', commandType='DP3')
### DONE

# # polarisation imaging
# with w.if_todo('imaging-pol'):
#     logger.info('Cleaning (Pol)...')
#     imagenameP = 'img/wideP'
#     lib_util.run_wsclean(s, 'wscleanP.log', MSs.getStrWsclean(), name=imagenameP, pol='QUV',
#         size=imgsizepix_p, scale='10arcsec', weight='briggs -0.3', niter=0, no_update_model_required='',
#         parallel_gridding=2, baseline_averaging='', minuv_l=30, maxuv_l=4500,
#         join_channels='', channels_out=MSs.getChout(4.e6))
#
# MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN SUBFIELD_DATA, FR_CORRECTED_DATA"',
#         log='$nameMS_taql_delcol.log', commandType='general')

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image*.fits self/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-MFS-residual*.fits self/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources*.txt self/images') for c in range(maxIter) ]
# os.system('mv img/wideP-MFS-*-image.fits self/images')
# os.system('mv img/wide-lr-MFS-image.fits self/images')

# Copy model
os.system(f'mv img/wideM-{maxIter-1}-*-model.fits self/skymodel')
os.system(f'mv img/wideM-{maxIter-1}-*-model-pb.fits self/skymodel')

logger.info("Done.")
