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

# TODO add BDA
# TODO the subfield algorithm should not cut any sources ... how to best implement that? Something with mask islands?
# TODO do we need amplitude corrections(?)? If so, DI/DD/only for bright sources?
# TODO final imaging products

# Waiting for bug fixes in other software
# TODO add facet-beam in imaging and predict steps once wsclean bug is fixed!
# TODO add LoTSS query for statring model once bug is fixed! (Don't use for now, it crashes the VO server)

import sys, os, glob, random
import numpy as np
from regions import Regions
import astropy.units as u
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
subfield_min_flux = parset.getfloat('LOFAR_self','subfield_min_flux') # default 20 Jy
subfield = parset.get('LOFAR_self','subfield') # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
maxIter = parset.getint('LOFAR_self','maxIter') # default = 2 (try also 3)
phaseSolMode = parset.get('LOFAR_self', 'ph_sol_mode') # tecandphase, tec, phase
intrinsic = parset.getboolean('LOFAR_self', 'intrinsic') # True means using intrinsic sky model and DP3 applybeam predict, False means staying with apparent skymodel
start_sourcedb = parset.get('model','sourcedb')
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')

#############################################################################

def clean_empty(MSs, name, col='CORRECTED_DATA', size=5000):
    """ For testing/debugging only"""
    lib_util.run_wsclean(s, 'wsclean-empty.log', MSs.getStrWsclean(), name=f'img/{name}',
                         data_column=col, size=size, scale='8arcsec', niter=0, nmiter=0,
                         weight='briggs 0.0', gridder='wgridder', parallel_gridding=1,
                         no_update_model_required='')

def corrupt_model_dirs(MSs, c, tc, model_columns):
    """ CORRUPT the MODEL_DATA columns for model_columns
    Parameters
    ----------
    MSs: MSs object
    c: int, cycle
    tc: tec/phase step (e.g. 2 for fast RS, 1 for intermediate, 0 for slow CS
    model_columns: list of patch names
    """
    logger.info(f'Corrupt models: {phaseSolMode}{tc}...')
    for model_column in model_columns:
        if phaseSolMode in ['tec', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}]  \
                    cor.parmdb=self/solutions/cal-tec{tc}-c{c}.h5 cor.correction=tec000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        if phaseSolMode in ['phase', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                    cor.parmdb=self/solutions/cal-tec{tc}-c{c}.h5 cor.correction=phase000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')

def solve_iono(MSs, c, tc, model_columns, smMHz, solint, soltype='phase', resetant_parset=None, model_column_fluxes=None, variable_solint_threshold=None):
    """
    Parallel solve for ionosphere systematics
    Parameters
    ----------
    MSs: MSs object
    c: major cycle
    tc: tec/phase solve step
    model_columns: columns to solve against
    smMHz: smoothnessconstraint kernel size at 54 MHz in MHz
    solint: solution interval in timesteps (e.g. multiples of 8s)
    resetant_parset: string, optional. Extra losoto parset to reset solutions for some stations
    model_column_fluxes: list of float, optional. Default=None. List of flux densities per direction/model column, used for variable solint.
    variable_solint_threshold: float, optional. Default = None. Use twice as long solint for directions that have less than this threshold in Jy.
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

    solutions_per_direction = np.ones(len(model_columns), dtype=int)
    if variable_solint_threshold: # if actived, use twice the solint for fainter directions
        solint *= 2 # use twice as long solint
        # get two solutions per solint (i.e. one per time step) for bright directions
        solutions_per_direction += model_column_fluxes > variable_solint_threshold

    MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec{tc}.h5 sol.solint={solint} {solver_params} \
                  sol.modeldatacolumns="[{",".join(model_columns)}]" sol.solutions_per_direction={np.array2string(solutions_per_direction, separator=",")}',
                  log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

    lib_util.run_losoto(s, f'tec{tc}-c{c}', [ms+f'/tec{tc}.h5' for ms in MSs.getListStr()], losoto_parsets, plots_dir=f'self/plots/plots-tec{tc}-c{c}', h5_dir=f'self/solutions/')


def make_current_best_mask(imagename, threshold=6.5, userReg=None):
    current_best_mask = f'{imagename}-mask.fits'
    if userReg:
        logger.info(f'Making mask with userReg {userReg}...')
        s.add(f'breizorro.py -t {threshold} -r {imagename}-MFS-image.fits -b 50 -o {current_best_mask} --merge {userReg}',
              log=f'makemask-{c}.log', commandType='python')
    else:
        logger.info('Making mask...')
        s.add(f'breizorro.py -t {threshold} -r {imagename}-MFS-image.fits -b 50 -o {current_best_mask}',
              log=f'makemask-{c}.log', commandType='python')
    s.run(check=True)
    return current_best_mask

def make_source_regions(sm, c):
    lib_util.check_rm(f'self/skymodel/regions_c{c}')
    lib_util.check_rm(f'self/skymodel/sources_c{c}.reg')
    os.makedirs(f'self/skymodel/regions_c{c}')
    for p in sm.getPatchNames():
        sm_p = sm.copy()
        sm_p.select(f'patch=={p}')
        sm_p.write(f'self/skymodel/regions_c{c}/{p}.reg', format='ds9', clobber=True)
        regs = Regions.read(f'self/skymodel/regions_c{c}/{p}.reg')
        col = '#{:06x}'.format(random.randint(0, 256 ** 3))
        for reg in regs:
            reg.visual['facecolor'] = col
            reg.visual['edgecolor'] = col
        regs.write(f'self/skymodel/regions_c{c}/{p}.reg',overwrite=True)
        os.system(f'cat self/skymodel/regions_c{c}/*.reg >> self/skymodel/sources_c{c}.reg')


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
beamMask = 'self/beam.fits'

# set image size
imgsizepix_wide = int(2.1*MSs.getListObj()[0].getFWHM(freq='mid')*3600/4.)
imgsizepix_lr = int(5*MSs.getListObj()[0].getFWHM(freq='mid')*3600/30.)
current_best_mask = None

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

mask_threshold = [6.5,5.5,4.5,4.0,4.0,4.0] # sigma values for beizorro mask in cycle c
# define list of facet fluxes per iteration -> this can go into the config
facet_fluxes = [5, 2, 1.4, 0.9, 0.6] # still needs to be tuned, maybe also depends on the field
# Try increasing number earlier
#facet_fluxes = [5, 1.5, 0.8, 0.7, 0.5] # still needs to be tuned, maybe also depends on the field
# define smoothness kernels in MHz (centered at 54 MHz)

# TODO try these kernels
# smMHz2 = [0.8,2.5,2.5,2.5,2.5,2.5]
# smMHz1 = [2.0,5.0,5.0,5.0,5.0,5.0]
# smMHz0 = [6.0,20.0,20.0,20.0,20.0,20.0]

smMHz2 = [1.0,5.0,5.0,5.0,5.0,5.0]
smMHz1 = [2.0,8.0,8.0,8.0,8.0,8.0]
smMHz0 = [6.0,20.0,20.0,20.0,20.0,20.0]
#################################################################

# Make beam mask
if not os.path.exists(beamMask):
    logger.info('Making mask of primary beam...')
    lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name=beamMask.replace('.fits',''), size=imgsizepix_lr, scale='30arcsec')
    os.system(f'mv {beamMask.replace(".fits","-image.fits")} {beamMask}')
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 1.)
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 0., inverse=True)

if not os.path.exists('self/solutions/cal-template.h5'):
    # create a template h5parm file
    if phaseSolMode == 'phase':  # phase
        solver_params = f'sol.mode=scalarphase sol.nchan=1'
    else:  # TEC or TecAndPhase
        solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
    MSs.run(
        f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/template.h5 sol.solint=1 sol.maxiter=1 {solver_params} sol.modeldatacolumns="[DATA]"',
        log='$nameMS_soltemplate.log', commandType='DP3')
    lib_util.run_losoto(s, 'template', [ms + '/template.h5' for ms in MSs.getListStr()], [])
    os.system('mv cal-template.h5 self/solutions/')
    os.system('rm -r plots-template')
#################################################################################################

with w.if_todo('solve_fr'):
    # Probably we do not need BLsmoothing since we have long time intervals and smoothnessconstraint?
    logger.info('Converting to circular (DATA -> CIRC_PHASEDIFF_DATA)...')
    # lin2circ conversion TCxx.MS:DATA -> CIRC_PHASEDIFF_DATA # use no dysco here!
    MSs.run(f'DP3 {parset_dir}/DP3-lin2circ.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CIRC_PHASEDIFF_DATA', log='$nameMS_lin2circ.log', commandType="DP3")

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

    # Solve TCXX.MS:CIRC_PHASEDIFF_DATA against MODEL_DATA_FR (only solve - solint=2m nchan=0 as it has the smoothnessconstrain)
    logger.info('Solving circ phase difference ...')
    MSs.run('DP3 ' + parset_dir + '/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.solint=' + str(30 * base_solint),
            log='$nameMS_solFR.log', commandType="DP3")
    lib_util.run_losoto(s, f'fr', [ms + '/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
    os.system('mv cal-fr.h5 self/solutions/')
    os.system('mv plots-fr self/plots/')
    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, MODEL_DATA_FR"',
            log='$nameMS_taql_delcol.log', commandType='general')

with w.if_todo('cor_fr'):
    # Correct FR with results of solve - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA_FR
    logger.info('Correcting FR (DATA -> CORRECTED_DATA_FR...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA_FR \
                cor.parmdb=self/solutions/cal-fr.h5 cor.correction=rotationmeasure000',
            log='$nameMS_corFR.log', commandType='DP3')
### DONE

#####################################################################################################
# Self-cal cycle
for c in range(maxIter):
    logger.info('Start selfcal cycle: '+str(c))
    if c == 0:
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Creating CORRECTED_DATA = CORRECTED_DATA_FR...')
            MSs.addcol('CORRECTED_DATA', 'CORRECTED_DATA_FR')
    elif c in [1,2] :
        # Only after the first two iterations: apply the subfield solutions to the data.
        with w.if_todo('set_corrected_data_c%02i' % c):
            logger.info('Set CORRECTED_DATA = CORRECTED_DATA_FR...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA_FR"', log='$nameMS_taql-c' + str(c) + '.log',
                    commandType='general')
            logger.info('Correct subfield ionosphere (CORRECTED_DATA -> CORRECTED_DATA)...')
            # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                        cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c-1) + '.h5 cor.correction=tec000 ',
                    log='$nameMS_sf-correct.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA  \
                        cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c-1) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')

            # if c == 2:
            #     logger.info('Correct subfield amplitude (CORRECTED_DATA -> CORRECTED_DATA)...')
            #     # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
            #     MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_CORRECTED_DATA \
            #             cor.parmdb=self/solutions/cal-amp-sf-c' + str(c-1) + '.h5 cor.correction=amplitude000',
            #             log='$nameMS_sf-correct.log', commandType='DP3')

    # get sourcedb
    sourcedb = f'tgts-c{c}.skymodel'
    beamMS = MSs.getListStr()[int(len(MSs.getListStr()) / 2)] # use central MS, should not make a big difference
    if not os.path.exists(sourcedb):
        logger.info(f'Create skymodel {sourcedb}')
        if c == 0:
            # if provided, use manual model
            if start_sourcedb:
                logger.info(f'Using input skymodel {start_sourcedb}')
                sm = lsmtool.load(start_sourcedb, beamMS=beamMS)
            else:
                # Create initial sourcedb from LoTSS or GSM
                fwhm = MSs.getListObj()[0].getFWHM(freq='min')
                try:
                    raise KeyError
                    # sm = lsmtool.load('LoTSS', VOPosition=phasecentre, VORadius=1.3 * fwhm / 2,
                    #                   beamMS=beamMS)
                except (KeyError, FileNotFoundError) as e:
                    logger.warning('Could not retrieve LoTSS data - resort to GSM.')
                    os.system(f'wget -O {sourcedb} "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord={phasecentre[0]},{phasecentre[1]}&radius={1.3*fwhm/2}&unit=deg"')
                    sm = lsmtool.load(sourcedb, beamMS=beamMS)
            # turn to apparent sky
            if not intrinsic:
                sm.write(sourcedb, clobber=True, applyBeam=True, adjustSI=True)
                sm = lsmtool.load(sourcedb)
        else:
            # get wsclean skymodel of last iteration
            wsc_src = f'img/wideM-{c-1}-pb-sources.txt' if intrinsic else f'img/wideM-{c-1}-sources.txt'
            sm = lsmtool.load(wsc_src, beamMS=beamMS if intrinsic else None)
        bright_sources_flux = facet_fluxes[c]
        # if using e.g. LoTSS, adjust for the frequency
        logger.info(f'Extrapolating input skymodel fluxes from {np.mean(MSs.getFreqs())/1e6:.0f}MHz to {sm.getDefaultValues()["ReferenceFrequency"]/1e6:.0f}MHz assuming si=-0.7')
        si_factor = (np.mean(MSs.getFreqs())/sm.getDefaultValues()['ReferenceFrequency'])**0.7 # S144 = si_factor * S54
        print(si_factor)
        sm.select(f'I>{0.05*si_factor}', applyBeam=intrinsic)  # keep only reasonably bright sources
        sm.select(f'{beamMask}==True')  # remove outside of FoV (should be subtracted (c>0) or not present (c==0)!)
        sm.group('threshold', FWHM=5/60, root='Src') # group nearby components to single source patch
        sm.setPatchPositions(method='wmean', applyBeam=intrinsic)
        sm = lib_dd_parallel.merge_nearby_bright_facets(sm, 1/60, 0.5, applyBeam=intrinsic)
        # TODO we need some logic here to avoid picking up very extended sources. Also case no bright sources in a field.
        patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=intrinsic)
        if sum(patch_fluxes/si_factor > bright_sources_flux) < 4:
            bright_sources_flux = np.sort(patch_fluxes)[-4]/si_factor
            logger.warning(f'Not enough bright sources flux! Using sources above {bright_sources_flux:.2f} Jy')
        bright_names = sm.getPatchNames()[patch_fluxes > bright_sources_flux*si_factor]
        bright_pos = sm.getPatchPositions(bright_names)
        sm.group('voronoi', targetFlux=bright_sources_flux*si_factor, applyBeam=intrinsic, root='', byPatch=True)
        sm.setPatchPositions(bright_pos)
        sm.plot(f'self/skymodel/patches-c{c}.png', 'patch')
        make_source_regions(sm, c)
        sm.write(sourcedb, clobber=True)
    else:
        logger.info(f'Load existing skymodel {sourcedb}')
        sm = lsmtool.load(sourcedb, beamMS=beamMS if intrinsic else None)

    # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
    sourcedb_basename = sourcedb.split('/')[-1]
    for MS in MSs.getListStr():
        lib_util.check_rm(MS + '/' + sourcedb_basename)
        logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
        os.system('cp -r ' + sourcedb + ' ' + MS)

    # Add model to MODEL_DATA
    patches = sm.getPatchNames()
    patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=intrinsic)

    with w.if_todo(f'init_model_c{c}'):
        for patch in patches:
            logger.info(f'Add model to {patch}...')
            pred_parset = 'DP3-predict-beam.parset' if intrinsic else 'DP3-predict.parset'
            MSs.run(f'DP3 {parset_dir}/{pred_parset} msin=$pathMS pre.sourcedb=$pathMS/{sourcedb_basename} pre.sources={patch} msout.datacolumn={patch}',
            log='$nameMS_pre.log', commandType='DP3')
            # Smooth CORRECTED_DATA -> SMOOTHED_DATA
            # MSs_sol.run_Blsmooth(patch, patch, logstr=f'smooth-c{c}')
            # pos = sm.getPatchPositions()[patch]
            # size = int((1.1*sm.getPatchSizes()[np.argwhere(sm.getPatchNames()==patch)]) // 4)
            # logger.info(f'Test image MODEL_DATA...')
            # lib_util.run_wsclean(s, 'wscleanMODEL-c' + str(c) + '.log', MSs.getStrWsclean(), name=f'self/skymodel/{patch}',
            #                      data_column=patch, size=size, scale='8arcsec', shift=f'{pos[0].to(u.hourangle).to_string()} {pos[1].to_string()}',
            #                      weight='briggs -0.3', niter=10000, gridder='wgridder', parallel_gridding=6, no_update_model_required='', minuv_l=30, mgain=0.9,
            #                      parallel_deconvolution=512, beam_size=15, join_channels='', fit_spectral_pol=3,
            #                      channels_out=MSs.getChout(4.e6), deconvolution_channels=3, pol='i', nmiter=5 )
        ### DONE

    with w.if_todo('solve_tecRS_c%02i' % c):
        # Smooth MSs:DATA -> SMOOTHED_DATA
        MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
        # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for all BUT most distant RS!
        logger.info('Solving TEC (fastRS)...')
        solve_iono(MSs, c, 2, patches, smMHz2[c], 2*base_solint, resetant_parset=parset_dir+'/losoto-resetph2-allCS.parset', model_column_fluxes=patch_fluxes, variable_solint_threshold=8.) # '/losoto-resetph2-CSRS.parset'
    ### DONE

    ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('corrupt_tecRS_c%02i' % c):
        corrupt_model_dirs(MSs, c, 2, patches)
    ### DONE

    with w.if_todo('solve_tecCS_c%02i' % c):
        # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for central CS
        logger.info('Solving TEC (midRS)...')
        solve_iono(MSs, c, 1, patches, smMHz1[c], 8*base_solint)#, resetant_parset=parset_dir+'/losoto-resetph2-CS.parset')
    ### DONE

    ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('corrupt_tecCS_c%02i' % c):
        corrupt_model_dirs(MSs, c, 1, patches)
    ### DONE

    # with w.if_todo('solve_tecCS0_c%02i' % c):
    #     # solve ionosphere phase - ms:SMOOTHED_DATA
    #     logger.info('Solving TEC (slowCS)...')
    #     solve_iono(MSs, c, 0, patches, smMHz0[c], 16*base_solint)
    # ### DONE


    # merge solutions into one h5parms for large scale image
    with w.if_todo(f'merge_h5_c{c}'):
        sol_dir = 'self/solutions'
        lib_util.check_rm(f'{sol_dir}/cal-tec-merged-c{c}.h5')
        # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
        # lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec0-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec1-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec2-c{c}.h5', sourcedb)
        # reference, unflag and reset the added stations # {sol_dir}/cal-tec0-c{c}.h5',
        h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-tec-merged-c{c}.h5", min_distance=1/3600, propagate_flags=True,
                          h5_tables=[f'{sol_dir}/cal-tec1-c{c}.h5',f'{sol_dir}/cal-tec2-c{c}.h5'],
                          h5_time_freq=f'{sol_dir}/cal-tec2-c{c}.h5', no_pol=True, ms_files='mss/TC*.MS',no_antenna_crash=True)
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
        reuse_kwargs = {}
        # make quick image with -local-rms to get a mask
        # TODO make this faster - experiment with increased parallel-gridding as well as shared facet reads option
        if c==0:
            logger.info('Making wide-field image for clean mask...')
            lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, data_column='CORRECTED_DATA', size=imgsizepix_wide, scale='4arcsec',
                                 weight='briggs -0.3', local_rms='', niter=1000000, gridder='wgridder',  parallel_gridding=32, no_update_model_required='', minuv_l=30, mgain=0.85, parallel_deconvolution=1024,
                                 auto_threshold=5.0, auto_mask=8.0, beam_size=15, join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                                 multiscale='', pol='i', nmiter=6,   facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000' )
            # make initial mask
            current_best_mask = make_current_best_mask(imagename, mask_threshold[c], userReg)
            # safe a bit of time by reusing psf and dirty in first iteration
            reuse_kwargs = {'reuse_psf':imagename, 'reuse_dirty':imagename}
        else:
            current_best_mask = f'img/wideM-{c}-mask.fits'


        if intrinsic or c>1:
            # if we are done with the things that require blanked pedict, we can also use beam keywords
            beam_kwargs = {'apply_facet_beam':'', 'facet_beam_update':120, 'use_differential_lofar_beam':''}
        else:
            beam_kwargs = {}
        # clean again, with mask now
        logger.info('Making wide field image ...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM,  fits_mask=current_best_mask, data_column='CORRECTED_DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder',  parallel_gridding=32, save_source_list='',
                             update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=4.0,
                             join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='',  multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80', pol='i',
                             facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000 ',
                             **reuse_kwargs, **beam_kwargs)
        # make a new mask from the image
        current_best_mask = make_current_best_mask(imagenameM, mask_threshold[c], userReg)

        # reset NaNs if present
        im = lib_img.Image(f'{imagenameM}-MFS-image.fits')
        im.nantozeroModel()

        # TODO wsclean facet pb predict seems bugged right now!
        # logger.info('Set MODEL_DATA=0...')
        # MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
        # logger.info('Predict corrupted MODEL_DATA...')
        # s.add(f'wsclean -predict -padding 1.8 -name {imagenameM} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
        #        -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged-c{c}.h5 phase000 \
        #         -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam {MSs.getStrWsclean()}',
        #        log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
        # s.run(check=True)
    #####################################################################################################
    ########################### TESTING FOR DD-AMP-CAL ####################################
    # if c == 2:
    #     with w.if_todo('solve_amp_c%02i' % c):
    #         logger.info('Solving amp (slowCS)...')
    #         MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset memory_logging=True numthreads=32 msin=$pathMS sol.datause=dual sol.nchan=24 \
    #                       sol.mode=diagonal sol.h5parm=$pathMS/amp4.h5 sol.solint={300 * base_solint} sol.modeldatacolumns="[{",".join(patches)}]"',
    #                 log='$nameMS_solamp-c' + str(c) + '.log', commandType='DP3', maxThreads=1)
    #
    #
    #         lib_util.run_losoto(s, f'amp4', [ms + f'/amp4.h5' for ms in MSs.getListStr()],
    #                             [f'{parset_dir}/losoto-plot-amp.parset', f'{parset_dir}/losoto-plot-ph.parset'],
    #                             plots_dir=f'self/plots/plots-amp4', h5_dir=f'self/solutions/')
    #
    #     with w.if_todo(f'merge_h5_2_c{c}'):
    #         sol_dir = 'self/solutions'
    #         lib_util.check_rm(f'{sol_dir}/cal-merged-c{c}.h5')
    #         # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
    #         lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-amp4.h5', sourcedb)
    #         # reference, unflag and reset the added stations
    #         h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-merged-c{c}.h5", min_distance=1/3600,
    #                            h5_tables=[f"{sol_dir}/cal-tec-merged-c{c}.h5",f"{sol_dir}/cal-amp4.h5"], propagate_flags=True,
    #                            h5_time_freq=f"{sol_dir}/cal-tec-merged-c{c}.h5",ms_files='mss/TC*.MS',no_antenna_crash=True)
    #         lib_util.run_losoto(s, f'merged-c{c}', f'{sol_dir}/cal-merged-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
    #                             plots_dir=f'self/plots/plots-merged-c{c}', h5_dir='self/solutions')
    #
    #     sys.exit()
    #     logger.info('Deleting model columns...')
    #     for patch in patches:
    #         MSs.run(f'taql "ALTER TABLE $pathMS DELETE COLUMN {patch}"',
    #                 log='$nameMS_taql_delcol.log', commandType='general')
    #
    #     imagenameM = 'img/wideM-2amp4'
    #     facetregname = f'self/solutions/facets-c{c}.reg'
    #     lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), name=imagenameM,
    #                          fits_mask='img/wideM-2-mask.fits', data_column='CORRECTED_DATA', size=imgsizepix_wide,
    #                          scale='4arcsec',
    #                          weight='briggs -0.3', niter=1000000, gridder='wgridder', parallel_gridding=32,
    #                          no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12,
    #                          parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=4.0,
    #                          join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6),
    #                          deconvolution_channels=3,
    #                          multiscale='', multiscale_scale_bias=0.7, multiscale_scales='0,10,20', pol='i',
    #                          apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='',
    #                          facet_regions=facetregname, apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000 self/solutions/cal-amp3.h5 amplitude000',
    #                          )
    #     sys.exit()
    ########################### TESTING FOR DI-AMP-CAL ####################################
    # if c == 2:
    #     MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
    #     # solve ionosphere phase - ms:SMOOTHED_DATA
    #     logger.info('Solving amp (slowCS)...')
    #     MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=dual sol.nchan=12 sol.modeldatacolumns=[MODEL_DATA] sol.mode=diagonal sol.h5parm=$pathMS/amp2.h5 sol.solint={300 * base_solint} ',
    #         log='$nameMS_solamp-c' + str(c) + '.log', commandType='DP3')
    #
    #     lib_util.run_losoto(s, f'amp2', [ms + f'/amp2.h5' for ms in MSs.getListStr()],
    #                         [f'{parset_dir}/losoto-plot-amp.parset', f'{parset_dir}/losoto-plot-ph.parset'],
    #                         plots_dir=f'self/plots/plots-amp2', h5_dir=f'self/solutions/')
    #
    #     logger.info('Correct subfield amplitude (CORRECTED_DATA -> CORRECTED_DATA)...')
    #     # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
    #     MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA2 \
    #             cor.parmdb=self/solutions/cal-amp2.h5 cor.correction=amplitude000',
    #             log='$nameMS_sf-correct.log', commandType='DP3')
    #     logger.info('Making wide field image ...')
    #     imagenameM = 'img/wideM-2amp2'
    #     lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), name=imagenameM,
    #                          fits_mask='img/wideM-2-mask.fits', data_column='CORRECTED_DATA2', size=imgsizepix_wide,
    #                          scale='4arcsec',
    #                          weight='briggs -0.3', niter=1000000, gridder='wgridder', parallel_gridding=32,
    #                          no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12,
    #                          parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=4.0,
    #                          join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6),
    #                          deconvolution_channels=3,
    #                          multiscale='', multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80', pol='i',
    #                          apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='',
    #                          facet_regions=facetregname, scalar_visibilities='',
    #                          apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000',
    #                          )
    #####################################################################################################
    # Find calibration solutions for subfield
    if c < 2:
        # Prepare region and models for subfield
        if subfield:
            subfield_path = subfield
            if len(Regions.read(subfield_path)) > 1:
                raise ValueError(f'Manual subfield region {subfield} contains more than one region')
        else:
            subfield_path = 'self/skymodel/subfield.reg'

        with w.if_todo('extreg_prepare_c%02i' % c):
            if not subfield and not os.path.exists(subfield_path): # automatically find subfield
                sm = lsmtool.load(f'img/wideM-{c}-sources.txt')
                # sm.remove('img/wide-lr-mask.fits=1')  # remove sidelobe sources that were subtracted
                sm.remove('MajorAxis > 80')  # remove largest scales
                field_center1, field_size1 = lib_dd.make_subfield_region(subfield_path, MSs.getListObj()[0], sm,
                                                                         subfield_min_flux, debug_dir='img/')
            # prepare model of central/external regions
            logger.info('Blanking central region of model files and reverse...')
            for im in glob.glob(f'img/wideM-{c}*model*.fits'):
                wideMint = im.replace('wideM','wideMint')
                os.system('cp %s %s' % (im, wideMint))
                lib_img.blank_image_reg(wideMint, subfield_path, blankval = 0., inverse=True)
                wideMext = im.replace('wideM','wideMext')
                os.system('cp %s %s' % (im, wideMext))
                lib_img.blank_image_reg(wideMext, subfield_path, blankval = 0.)
        # DONE
        subfield_reg = Regions.read(subfield_path)[0]
        field_center = subfield_reg.center.ra, subfield_reg.center.dec
        field_size = np.max([subfield_reg.width.to_value('deg'), subfield_reg.height.to_value('deg')])

        with w.if_todo('extreg_subtract_c%02i' % c):
            # Recreate MODEL_DATA of external region for subtraction
            logger.info('Set MODEL_DATA=0...')
            MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
            logger.info('Predict corrupted model of external region...')
            # TODO ignore facet beam for now due to bug
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMext-{c} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                    -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged-c{c}.h5 phase000 \
                    {MSs.getStrWsclean()}',
                log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
            s.run(check=True)

            # cycle >0: need to add DI-corruption on top (previous iteration sub-field)
            if c > 0:
                logger.info('Add previous iteration sub-field corruption on top of DD-corruption...')
                if phaseSolMode in ['tec', 'tecandphase']:
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                            cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c-1) + '.h5 cor.correction=tec000 cor.invert=False',
                            log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
                if phaseSolMode in ['phase', 'tecandphase']:
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                            cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c-1) + '.h5 cor.correction=phase000 cor.invert=False',
                            log='$nameMS_sidelobe_corrupt.log', commandType='DP3')

            # subtract external region from CORRECTED_DATA_FR to create SUBFIELD_DATA
            MSs.addcol('SUBFIELD_DATA','CORRECTED_DATA_FR')
            logger.info('Subtracting external region model (SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
            clean_empty(MSs,'only_subfield', 'SUBFIELD_DATA')
        ### DONE

        with w.if_todo('intreg_predict%02i' % c):
            # Recreate MODEL_DATA of internal region for solve
            logger.info('Predict model of internal region...')
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMint-{c} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                    {MSs.getStrWsclean()}',
                  log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
            s.run(check=True)
        ### DONE

        with w.if_todo('subfield_solve_tecRS_c%02i' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving TEC (fast RS)...')
            solve_iono(MSs, c, '2-sf', ['MODEL_DATA'], 0.7, base_solint, resetant_parset=parset_dir + '/losoto-resetph2-CSRS.parset')
        ### DONE

        with w.if_todo('subfield_corr_tecRS_c%02i' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        with w.if_todo('subfield_solve_tecCS_c%02i' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving TEC (mid RS)...')
            solve_iono(MSs, c, '1-sf', ['MODEL_DATA'], 3.0, 4*base_solint, resetant_parset=parset_dir + '/losoto-resetph2-CS.parset')
        ### DONE

        with w.if_todo('subfield_corr_tecCS_c%02i' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        with w.if_todo('subfield_solve_tecCS0_c%02i' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving TEC (slowCS)...')
            solve_iono(MSs, c, '0-sf', ['MODEL_DATA'], 6.0, 16*base_solint)
        ### DONE

        with w.if_todo('subfield_corr_tecCS0_c%02i' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec0-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec0-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        # merge solutions into one h5parms for plotting and apply
        with w.if_todo(f'merge_h5_subfield_c{c}'):
            sol_dir = 'self/solutions'
            lib_util.check_rm(f'{sol_dir}/cal-tec-sf-merged-c{c}.h5')
            # reference, unflag and reset the added stations
            h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-tec-sf-merged-c{c}.h5", min_distance=10 / 3600,
                               h5_tables=[f'{sol_dir}/cal-tec0-sf-c{c}.h5', f'{sol_dir}/cal-tec1-sf-c{c}.h5',
                                          f'{sol_dir}/cal-tec2-sf-c{c}.h5'],
                               h5_time_freq=f'{sol_dir}/cal-template.h5', no_pol=True, ms_files='mss/TC*.MS',
                               no_antenna_crash=True)
            lib_util.run_losoto(s, f'tec-sf-merged-c{c}', f'{sol_dir}/cal-tec-sf-merged-c{c}.h5',
                                [f'{parset_dir}/losoto-plot-scalar.parset'],
                                plots_dir=f'self/plots/plots-tec-sf-merged-c{c}', h5_dir='self/solutions')
        ### DONE

        # Do a quick debug image...
        with w.if_todo(f'image-subfield-c{c}'):
            logger.info('Test image subfield...')
            lib_util.run_wsclean(s, 'wscleanSF-c'+str(c)+'.log', MSs.getStrWsclean(), name=f'img/subfield-c{c}', data_column='SUBFIELD_DATA', size=3000, scale='4arcsec',
                                 weight='briggs -0.3', niter=100000, gridder='wgridder',  parallel_gridding=6, shift=f'{field_center[0].to(u.hourangle).to_string()} {field_center[1].to_string()}',
                                 no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                                 join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3, baseline_averaging='',
                                 multiscale='',  multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80',  pol='i')
        ### DONE
        if c>0:
            with w.if_todo('subfield_solve_amp%02i' % c):
                # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
                # MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
                # solve ionosphere phase - ms:SMOOTHED_DATA
                logger.info('Solving amp (slowCS)...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=dual sol.nchan=1 sol.modeldatacolumns=[MODEL_DATA] sol.mode=diagonal sol.h5parm=$pathMS/amp.h5 sol.solint={300*base_solint} ',
                log='$nameMS_solamp-c' + str(c) + '.log', commandType='DP3')

                lib_util.run_losoto(s, f'amp-sf-c{c}', [ms + f'/amp.h5' for ms in MSs.getListStr()], [f'{parset_dir}/losoto-plot-amp.parset',f'{parset_dir}/losoto-plot-ph.parset'],
                                    plots_dir=f'self/plots/plots-amp-sf-c{c}', h5_dir=f'self/solutions/')


            with w.if_todo('subfield_corr_amp%02i' % c):
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msout.datacolumn=SUBFIELD_CORRECTED_DATA \
                        cor.parmdb=self/solutions/cal-amp-sf-c' + str(c) + '.h5 cor.correction=amplitude000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
                # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_CORRECTED_DATA msout.datacolumn=SUBFIELD_CORRECTED_DATA \
                #         cor.parmdb=self/solutions/cal-amp-sf-c' + str(c) + '.h5 cor.correction=phase000',
                #         log='$nameMS_sf-correct.log', commandType='DP3')
                logger.info('Test image subfield...')
                lib_util.run_wsclean(s, 'wscleanSF-c'+str(c)+'.log', MSs.getStrWsclean(), name=f'img/subfield-amp-c{c}', data_column='SUBFIELD_CORRECTED_DATA', size=1000, scale='4arcsec',
                                     weight='briggs -0.3', niter=100000, gridder='wgridder',  parallel_gridding=6, shift=f'{field_center[0].to(u.hourangle).to_string()} {field_center[1].to_string()}',
                                     no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                                     join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3, baseline_averaging='',
                                     multiscale='',  multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80',  pol='i')
                sys.exit()
            ### DONE
        #####################################################################################################
        # Subtract side-lobe sources
        if c == 0:
            # Predict the main-lobe with DD-sols and subtract it
            with w.if_todo('subtract_mainlobe'):
                logger.info('Blanking central region of model files and reverse...')
                for im in glob.glob('img/wideM-0*model*.fits'):
                    wideMintpb = im.replace('wideM', 'wideMintpb')
                    os.system('cp %s %s' % (im, wideMintpb))
                    lib_img.blank_image_reg(wideMintpb, beamReg , blankval=0., inverse=True)
                    wideMextpb = im.replace('wideM', 'wideMextpb')
                    os.system('cp %s %s' % (im, wideMextpb))
                    lib_img.blank_image_reg(wideMextpb, beamReg, blankval=0.)
                logger.info('Set MODEL_DATA=0...')
                MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c' + str(c) + '.log')
                # Recreate MODEL_DATA of external region for subtraction -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam
                logger.info('Predict corrupted model of external region...')
                s.add(f'wsclean -predict -padding 1.8 -name img/wideMintpb-0 -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                       -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged-c{c}.h5 phase000 {MSs.getStrWsclean()}',
                    log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
                s.run(check=True)

                # TODO we use SUBFIELD_DATA here since it exists and is not needed, might be clearer to instead use SUBTRACTED_DATA or so?
                # subtract internal region from CORRECTED_DATA_FR to create SUBFIELD_DATA
                logger.info('Subtract main-lobe (SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA"',
                        log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
                clean_empty(MSs, 'only_sidelobe', 'SUBFIELD_DATA', size=10000)
            ### DONE

            # Do a rough correction of the sidelobe data using the subfield solutions
            with w.if_todo(f'correct-sidelobe-c{c}'): # just for testing/debug
                logger.info('Correct sidelobe data with subfield iono solutions...')
                if phaseSolMode in ['tec', 'tecandphase']:
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                            cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=tec000 ',
                            log='$nameMS_sf-correct.log', commandType='DP3')
                if phaseSolMode in ['phase', 'tecandphase']:
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                            cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000',
                            log='$nameMS_sf-correct.log', commandType='DP3')
            # DONE

            imagename_lr = 'img/wide-lr'
            # Image the sidelobe data
            with w.if_todo('image_sidelobe'):
                logger.info('Cleaning low-res...')
                lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=True, data_column='SUBFIELD_DATA',
                                     parallel_gridding=4, temp_dir='../', size=imgsizepix_lr, scale='30arcsec',
                                     weight='briggs -0.3', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
                                     taper_gaussian='200arcsec', mgain=0.85, channels_out=MSs.getChout(2.e6), parallel_deconvolution=512, baseline_averaging='',
                                     local_rms='', auto_mask=3, auto_threshold=1.5, join_channels='')

            # Subtract full low-resolution field (including possible large-scale emission within primary beam)
            # to get empty data set for flagging
            with w.if_todo('subtract_lr'):
                logger.info('Subtract low-resolution to get empty data set (SUBFIELD_DATA = SUBFIELD_DATA - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set SUBFIELD_DATA = SUBFIELD_DATA - MODEL_DATA"',
                    log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
                clean_empty(MSs, 'empty', 'SUBFIELD_DATA')

            # Flag on residuals (SUBFIELD_DATA)
            with w.if_todo('flag_residuals'):
                logger.info('Flagging residuals (SUBFIELD_DATA)...')
                MSs.run(
                    'DP3 ' + parset_dir + '/DP3-flag.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA aoflagger.strategy=' + parset_dir + '/LBAdefaultwideband.lua',
                    log='$nameMS_flag-c' + str(c) + '.log', commandType='DP3')

            # Now subtract the sidelobe
            with w.if_todo('subtract_sidelobe'):
                # blank within main-libe to not subtract anything from there
                for im in glob.glob(f'{imagename_lr}*model*fits'):
                    wideLRext = im.replace(imagename_lr, f'{imagename_lr}-blank')
                    os.system('cp %s %s' % (im, wideLRext))
                    lib_img.blank_image_reg(wideLRext, beamReg , blankval=0.)

                logger.info('Predict model of sidelobe region...')
                s.add(
                    f'wsclean -predict -padding 1.8 -name {imagename_lr}-blank -j {s.max_processors} -channels-out {MSs.getChout(2.e6)} {MSs.getStrWsclean()}',
                    log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
                s.run(check=True)

                logger.info('Corrupt sidelobe model with subfield solutions...')
                if phaseSolMode in ['tec', 'tecandphase']:
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                            cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=tec000 cor.invert=False',
                            log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
                if phaseSolMode in ['phase', 'tecandphase']:
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                            cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000 cor.invert=False',
                            log='$nameMS_sidelobe_corrupt.log', commandType='DP3')

                logger.info('Subtract corrupted sidelobe model (CORRECTED_DATA_FR = CORRECTED_DATA_FR - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set CORRECTED_DATA_FR = CORRECTED_DATA_FR - MODEL_DATA"', log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
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
# MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN SUBFIELD_DATA, CORRECTED_DATA_FR"',
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
