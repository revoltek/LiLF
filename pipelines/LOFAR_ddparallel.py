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

# TODO test the effect of scalar-visibilities on stoke I and V image quality
# TODO the subfield algorithm should not cut any sources ... how to best implement that? Something with mask islands?
# TODO final imaging products

# Waiting for bug fixes in other software
# TODO add facet-beam in imaging and predict steps once wsclean bug is fixed!
# TODO add LoTSS query for statring model once bug is fixed! (Don't use for now, it crashes the VO server)
# TODO add BDA

import sys, os, glob, random
import numpy as np
from regions import Regions
import astropy.units as u
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, lib_h5, lib_dd_parallel
logger_obj = lib_log.Logger('pipeline-ddparallel')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-ddparallel.walker')

parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_ddparallel'])))
parset_dir = parset.get('LOFAR_ddparallel','parset_dir')
subfield_min_flux = parset.getfloat('LOFAR_ddparallel','subfield_min_flux') # default 20 Jy
subfield = parset.get('LOFAR_ddparallel','subfield') # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
maxIter = parset.getint('LOFAR_ddparallel','maxIter') # default = 2 (try also 3)
phaseSolMode = parset.get('LOFAR_ddparallel', 'ph_sol_mode') # tecandphase, tec, phase
sf_phaseSolMode = 'phase' #'tec'
start_sourcedb = parset.get('model','sourcedb')
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')
subtract_predict_mode = 'wsclean'

#############################################################################

def clean_empty(MSs, name, col='CORRECTED_DATA', size=5000):
    """ For testing/debugging only"""
    lib_util.run_wsclean(s, 'wsclean-empty.log', MSs.getStrWsclean(), name=f'img/{name}',
                         data_column=col, size=size, scale=f'{int(pixscale*2)}arcsec', niter=0, nmiter=0,
                         weight='briggs 0.0', gridder='wgridder', parallel_gridding=1,
                         no_update_model_required='')

def corrupt_model_dirs(MSs, c, tc, model_columns, solmode='phase'):
    """ CORRUPT the MODEL_DATA columns for model_columns
    Parameters
    ----------
    MSs: MSs object
    c: int, cycle
    tc: tec/phase step (e.g. 2 for fast RS, 1 for intermediate, 0 for slow CS
    model_columns: list of patch names
    solmode: which solution mode to corrupt for (phase, tec, tecandphase)
    """
    logger.info(f'Corrupt models: {phaseSolMode}{tc}...')
    for model_column in model_columns:
        if solmode in ['tec', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}]  \
                    cor.parmdb=self/solutions/cal-tec{tc}-c{c}.h5 cor.correction=tec000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        elif solmode in ['phase', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                    cor.parmdb=self/solutions/cal-tec{tc}-c{c}.h5 cor.correction=phase000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        elif solmode in ['scalaramplitude']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                    cor.parmdb=self/solutions/cal-solamp-c{c}.h5 cor.correction=amplitude000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')

def solve_iono(MSs, c, tc, model_columns, smMHz, solint, solmode, resetant=None, constrainant=None, model_column_fluxes=None, variable_solint_threshold=None):
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
    solmode: srt, tec, phase, tecandphase
    resetant: string, optional. Reset either 'inner' or intermediate stations. Default None
    constrainant: string, None or RS -> constrain in solving
    model_column_fluxes: list of float, optional. Default=None. List of flux densities per direction/model column, used for variable solint.
    variable_solint_threshold: float, optional. Default = None. Use twice as long solint for directions that have less than this threshold in Jy.
    """
    # reset ants after solving if specified
    resetant_parset = None
    if solmode == 'phase':
        if resetant == 'intermediate':
            resetant_parset = parset_dir+'/losoto-resetph2-CSRS.parset'
        elif resetant == 'inner':
            resetant_parset = parset_dir+'/losoto-resetph2-CS.parset'
        elif resetant == 'CS':
            resetant_parset = parset_dir + '/losoto-resetph2-allCS.parset'
    elif solmode == 'tec':
        if resetant == 'intermediate':
            resetant_parset = parset_dir+'/losoto-resettec2-CSRS.parset'
        elif resetant == 'inner':
            resetant_parset = parset_dir+'/losoto-resettec2-CS.parset'
        elif resetant == 'CS':
            resetant_parset = parset_dir+'/losoto-resettec2-allCS.parset'

    if solmode == 'phase': #phase
        if resetant_parset is not None:
            losoto_parsets = [parset_dir+'/losoto-refph.parset', resetant_parset, parset_dir+'/losoto-plot-scalar.parset']
        else:
            losoto_parsets = [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-plot-scalar.parset']
    else: # TEC or TecAndPhase
        if resetant_parset is not None:
            losoto_parsets = [parset_dir+'/losoto-reftec.parset', resetant_parset, parset_dir+'/losoto-plot-tec.parset']
        else:
            losoto_parsets = [parset_dir+'/losoto-reftec.parset', parset_dir+'/losoto-plot-tec.parset']

    antennaconstraint = ''
    if constrainant is not None:
        if constrainant == 'RS':
            antennaconstraint = 'sol.antennaconstraint=[[RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]]'
        else:
            raise ValueError

    solutions_per_direction = np.ones(len(model_columns), dtype=int)
    if variable_solint_threshold and solmode == 'tec':
        raise ValueError
    elif variable_solint_threshold: # if actived, use twice the solint for fainter directions
        solint *= 2 # use twice as long solint
        # get two solutions per solint (i.e. one per time step) for bright directions
        solutions_per_direction += model_column_fluxes > variable_solint_threshold
    if solmode == 'phase':
        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec{tc}.h5 sol.solint={solint} \
                  sol.mode=scalarphase sol.smoothnessconstraint={smMHz}e6 sol.smoothnessreffrequency=54e6 sol.nchan=1 {antennaconstraint} \
                  sol.modeldatacolumns="[{",".join(model_columns)}]" sol.solutions_per_direction={np.array2string(solutions_per_direction, separator=",")}',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')
    else:
        MSs.run(f'DP3 {parset_dir}/DP3-solTEC.parset msin=$pathMS sol.h5parm=$pathMS/tec{tc}.h5 sol.solint={solint} {antennaconstraint} \
                  sol.modeldatacolumns="[{",".join(model_columns)}]" sol.mode={solmode}',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

    lib_util.run_losoto(s, f'tec{tc}-c{c}', [ms+f'/tec{tc}.h5' for ms in MSs.getListStr()], losoto_parsets, 
                        plots_dir=f'self/plots/plots-tec{tc}-c{c}', h5_dir=f'self/solutions/')
    
    
def solve_amplitude(MSs, c, model_columns, smMHz, solint, resetant=None, constrainant=None, model_column_fluxes=None, variable_solint_threshold=8.):
    solutions_per_direction = np.ones(len(model_columns), dtype=int)
    
    MSs.run(
        f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS \
            msin.datacolumn=SMOOTHED_DATA sol.mode=scalaramplitude sol.nchan=1\
            sol.h5parm=$pathMS/solamp_c{c}.h5 sol.modeldatacolumns="[{",".join(model_columns)}]" \
            sol.solutions_per_direction={np.array2string(solutions_per_direction.astype(int), separator=",")} \
            sol.smoothnessreffrequency=54e6 sol.solint={solint} sol.smoothnessconstraint={smMHz}e6',
        log=f'$nameMS_solamp_c{c}.log', 
        commandType="DP3"
    )
    
    losoto_parsets = [parset_dir+'/losoto-amp.parset', parset_dir+'/losoto-plot-amp.parset']
    
    lib_util.run_losoto(s, f'solamp-c{c}', [ms+f'/solamp_c{c}.h5' for ms in MSs.getListStr()], losoto_parsets, 
                        plots_dir=f'self/plots/plots-solamp-c{c}', h5_dir=f'self/solutions/')


def make_current_best_mask(imagename, threshold=6.5, userReg=None):
    current_best_mask = f'{imagename}-mask.fits'
    if userReg:
        logger.info(f'Making mask with userReg {userReg}...')
        s.add(f'breizorro.py -t {threshold} -r {imagename}-MFS-image.fits -b 70 -o {current_best_mask} --merge {userReg}',
              log=f'makemask-{c}.log', commandType='python')
    else:
        logger.info('Making mask...')
        s.add(f'breizorro.py -t {threshold} -r {imagename}-MFS-image.fits -b 70 -o {current_best_mask}',
              log=f'makemask-{c}.log', commandType='python')
    s.run(check=True)
    return current_best_mask

def add_3c_models(sm: lsmtool.skymodel.SkyModel, phasecentre=[0,0], fwhm=0, max_sep=30., threshold=1.):
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    
    Vizier.ROW_LIMIT = 200
    table = Vizier.get_catalogs("J/MNRAS/204/151")[0]
    all_3c = table['Name']
    
    phasecentre = SkyCoord(phasecentre[0], phasecentre[1], unit=(u.deg, u.deg))
    logger.info('Adding 3C models...')
    for i, source in enumerate(all_3c):
        pos = SkyCoord(table['_RA.icrs'][i], table['_DE.icrs'][i], unit=(u.hourangle, u.deg))
        if phasecentre.separation(pos).deg < fwhm/2 or phasecentre.separation(pos).deg > max_sep:
            continue 
        
        sourcedb = f'/home/local/work/j.boxelaar/LiLF//models/3CRR/{source.replace(" ","")}.txt'
        if not os.path.exists(sourcedb):
            logger.warning(f'No model found for {source} (seperation {phasecentre.separation(pos).deg:.2f} deg)')
            continue
        
        logger.info(f'Appending model from {sourcedb.split("/")[-1]} (seperation {phasecentre.separation(pos).deg:.2f} deg)...')
        sm_3c = lsmtool.load(sourcedb, beamMS=sm.beamMS)
        sm_3c.setColValues("Patch", ["source_"+source.replace(" ","")]*len(sm_3c.getColValues("I")))
        sm_3c.select(f'I>{threshold:.02f}', aggregate='sum', applyBeam=True)
        sm.concatenate(sm_3c)
        sm.setPatchPositions(method='wmean', applyBeam=True)
    return sm

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
    lib_util.check_rm('tgts*skymodel')

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

# make beam to the first mid null - outside of that do a rough subtraction and/or 3C peeling. Use sources inside for calibration
phasecentre = MSs.getListObj()[0].getPhaseCentre()

fwhm = MSs.getListObj()[0].getFWHM(freq='mid') * 1.8 # FWHM to null
MSs.getListObj()[0].makeBeamReg('self/beam.reg', freq='mid', to_pbval=0)

beamReg = 'self/beam.reg'
beamMask = 'self/beam.fits'

# set image size - this should be a bit more than the beam region used for calibration
pixscale = MSs.getListObj()[0].getPixelScale()
imgsizepix_wide = int(1.85*MSs.getListObj()[0].getFWHM(freq='mid')*3600/pixscale) # roughly to null
if imgsizepix_wide > 10000:
    imgsizepix_wide = 10000
imgsizepix_lr = int(5*MSs.getListObj()[0].getFWHM(freq='mid')*3600/(pixscale*8))
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
    base_solint = int(np.rint(4/tint)) # this is already 4 for dutch observations
else: base_solint = 1

mask_threshold = [5.0,4.5,4.0,4.0,4.0,4.0] # sigma values for beizorro mask in cycle c
# define list of facet fluxes per iteration -> this can go into the config
facet_fluxes = np.array([4, 1.8, 1.2, 0.8, 0.6])*(54e6/np.mean(MSs.getFreqs()))**0.7 # this is not the total flux, but the flux of bright sources used to construct the facets. still needs to be tuned, maybe also depends on the field
min_facets = [3, 6, 18, 24, 24, 24]

smMHz2 = [1.0,5.0,5.0,5.0,5.0,5.0]
smMHz1 = [5.0,8.0,8.0,8.0,8.0,8.0]
# smMHz0 = [6.0,10.0,10.0,10.0,10.0,10.0]
#################################################################

# Make beam mask
if not os.path.exists(beamMask):
    logger.info('Making mask of primary beam...')
    lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name=beamMask.replace('.fits',''), size=imgsizepix_lr, scale='30arcsec')
    os.system(f'mv {beamMask.replace(".fits","-image.fits")} {beamMask}')
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 1.)
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 0., inverse=True)

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
        with w.if_todo('c%02i_set_corrected_data' % c):
            logger.info('Creating CORRECTED_DATA = CORRECTED_DATA_FR...')
            MSs.addcol('CORRECTED_DATA', 'CORRECTED_DATA_FR')

    # get sourcedb
    sourcedb = f'tgts-c{c}.skymodel'
    beamMS = MSs.getListStr()[int(len(MSs.getListStr()) / 2)] # use central MS, should not make a big difference
    intrinsic = not c in [1,2] # these cycles we only have the apparent skymodel - predict w/o beam, find patches in apparent model
    if not os.path.exists(sourcedb):
        logger.info(f'Creating skymodel {sourcedb}...')
        if c == 0:
            # if provided, use manual model
            if start_sourcedb:
                logger.info(f'Using input skymodel {start_sourcedb}')
                sm = lsmtool.load(start_sourcedb, beamMS=beamMS)
            else:
                # Create initial sourcedb from GSM
                fwhm = MSs.getListObj()[0].getFWHM(freq='mid')*1.8 # to null
                logger.info('Get skymodel from GSM.')
                os.system(f'wget -O {sourcedb} "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord={phasecentre[0]},{phasecentre[1]}&radius={fwhm/2}&unit=deg"')
                sm = lsmtool.load(sourcedb, beamMS=beamMS)
        else:
            # get wsclean skymodel of previous iteration
            wsc_src = f'img/wideM-{c-1}-sources-pb.txt' if intrinsic else f'img/wideM-{c-1}-sources.txt'
            sm = lsmtool.load(wsc_src, beamMS=beamMS if intrinsic else None)
        # if using e.g. LoTSS, adjust for the frequency
        logger.debug(f'Extrapolating input skymodel fluxes from {sm.getDefaultValues()["ReferenceFrequency"]/1e6:.0f}MHz to {np.mean(MSs.getFreqs())/1e6:.0f}MHz assuming si=-0.7')
        si_factor = (np.mean(MSs.getFreqs())/sm.getDefaultValues()['ReferenceFrequency'])**0.7 # S144 = si_factor * S54
        # sm.select(f'I>{0.01*si_factor}', applyBeam=intrinsic)  # keep only reasonably bright sources
        sm.select(f'{beamMask}==True')  # remove outside of FoV (should be subtracted (c>0) or not present (c==0)!)
        sm.group('threshold', FWHM=5/60, root='Src') # group nearby components to single source patch
        sm.setPatchPositions(method='wmean', applyBeam=intrinsic)
        sm = lib_dd_parallel.merge_nearby_bright_facets(sm, 1/60, 0.5, applyBeam=intrinsic)
        # TODO we need some logic here to avoid picking up very extended sources. Also case no bright sources in a field.
        patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=intrinsic)
        if sum(patch_fluxes/si_factor > facet_fluxes[c]) < min_facets[c]:
            bright_sources_flux = np.sort(patch_fluxes)[-min_facets[c]] / si_factor
            logger.warning(f'Not enough bright sources above minimum flux {bright_sources_flux:.2f} Jy! Using sources above {bright_sources_flux:.2f} Jy')
        else:
            bright_sources_flux = facet_fluxes[c]
        bright_names = sm.getPatchNames()[patch_fluxes > bright_sources_flux*si_factor]
        bright_pos = sm.getPatchPositions(bright_names)
        sm.group('voronoi', targetFlux=bright_sources_flux*si_factor, applyBeam=intrinsic, root='', byPatch=True)
        sm.setPatchPositions(bright_pos)
        lib_dd_parallel.rename_skymodel_patches(sm, applyBeam=intrinsic)

        if c == 0:
            # Add models of bright 3c sources to the sky model. model will be subtracted from data before imaging.
            sm = add_3c_models(sm, phasecentre=phasecentre, fwhm=fwhm)
        
        sm.plot(f'self/skymodel/patches-c{c}.png', 'patch')
        make_source_regions(sm, c)
        sm.write(sourcedb, clobber=True)
        logger.info(f'Using {len(sm.getPatchNames())} patches.')
    else:
        logger.info(f'Load existing skymodel {sourcedb}')
        sm = lsmtool.load(sourcedb, beamMS=beamMS if intrinsic else None)

    # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
    sourcedb_basename = sourcedb.split('/')[-1]
    for MS in MSs.getListStr():
        lib_util.check_rm(MS + '/' + sourcedb_basename)
        logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
        os.system('cp -r ' + sourcedb + ' ' + MS)

    patches = sm.getPatchNames()
    patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=intrinsic)
    for patch in patches[np.argsort(patch_fluxes)[::-1]]:
        logger.info(f'{patch}: {patch_fluxes[patches==patch][0]:.1f}Jy')

    with w.if_todo('c%02i_init_model' % c):
        for patch in patches:
            # Add model to MODEL_DATA
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

    with w.if_todo('c%02i_solve_tecRS' % c):
        # Smooth MSs:CORRECTED_DATA -> SMOOTHED_DATA
        MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
        # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for all BUT most distant RS!
        logger.info('Solving TEC (RS)...')
        solve_iono(MSs, c, 2, patches, smMHz2[c], 2*base_solint, 'phase', resetant='CS', model_column_fluxes=patch_fluxes, variable_solint_threshold=8.)
    ### DONE

    ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('c%02i_corrupt_tecRS' % c):
        corrupt_model_dirs(MSs, c, 2, patches)
    ### DONE

    with w.if_todo('c%02i_solve_tecCS' % c):
        # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for central CS
        logger.info('Solving TEC (CS)...')
        solve_iono(MSs, c, 1, patches, smMHz1[c], 16*base_solint, 'phase', constrainant=None) # 'RS'
    ### DONE

    # ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('c%02i_corrupt_tecCS' % c):
        corrupt_model_dirs(MSs, c, 1, patches, 'phase')
    # ### DONE
    
    # with w.if_todo('c%02i_solve_amp' % c):
    #     # solve ionosphere phase - ms:SMOOTHED_DATA - > reset for central CS
    #     logger.info('Solving Amplitude...')
    #     solve_amplitude(MSs, c, patches, smMHz1[c], 16*base_solint, model_column_fluxes=patch_fluxes, variable_solint_threshold=8.)
    # ### DONE
    #
    # # ### CORRUPT the MODEL_DATA columns for all patches
    # with w.if_todo('c%02i_corrupt_tecCS' % c):
    #     corrupt_model_dirs(MSs, c, 1, patches, 'scalaramplitude')
    # ### DONE

    # # Only once in cycle 1: do di amp to capture element beam 2nd order effect
    if c == 1:
        with w.if_todo('solve_amp_di'):
            # TODO -> down the road this could be a fulljones-calibration to correct element-beam related leakage.
            MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
            logger.info('Setting MODEL_DATA to sum of corrupted patch models...')
            MSs.run(f'taql "UPDATE $pathMS SET MODEL_DATA={"+".join(patches)}"', log='$nameMS_taql_phdiff.log', commandType='general')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info('Solving amp-di...')
            MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=dual sol.nchan=12 sol.modeldatacolumns=[MODEL_DATA] \
                     sol.mode=diagonal sol.h5parm=$pathMS/amp-di.h5 sol.solint={150*base_solint} \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]',
                     log='$nameMS_solamp-c' + str(c) + '.log', commandType='DP3')

            lib_util.run_losoto(s, f'amp-di', [ms + f'/amp-di.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-amp.parset', f'{parset_dir}/losoto-plot-ph.parset', f'{parset_dir}/losoto-amp-di.parset'],
                                plots_dir=f'self/plots/plots-amp-di', h5_dir=f'self/solutions/')
            
        with w.if_todo('correct_amp_di'):
            # TODO add updateweights in production
            # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
            logger.info('Correct amp-di (CORRECTED_DATA -> CORRECTED_DATA)...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb=self/solutions/cal-amp-di.h5 cor.correction=amplitudeSmooth cor.updateweights=False',
                    log='$nameMS_sf-correct.log', commandType='DP3')

    # merge solutions into one h5parms for large scale image
    with w.if_todo('c%02i_merge_h5' % c):
        sol_dir = 'self/solutions'
        lib_util.check_rm(f'{sol_dir}/cal-tec-merged-c{c}.h5')
        # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
        # lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec0-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec1-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec2-c{c}.h5', sourcedb)
        # reference, unflag and reset the added stations f'{sol_dir}/cal-tec0-c{c}.h5',
        logger.info('Merge solutions...')
        s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-RS-c{c}.h5 --h5_tables {sol_dir}/cal-tec2-c{c}.h5 --h5_time_freq {sol_dir}/cal-tec2-c{c}.h5 \
              --no_antenna_crash --no_pol' , log='h5_merger.log', commandType='python')
        s.run(check=True)
        lib_util.run_losoto(s, f'tec-RS-c{c}', f'{sol_dir}/cal-tec-RS-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
                            plots_dir=f'self/plots/plots-tec-RS-c{c}', h5_dir='self/solutions')
        s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-CS-c{c}.h5 --h5_tables {sol_dir}/cal-tec1-c{c}.h5 --h5_time_freq {sol_dir}/cal-tec2-c{c}.h5 \
              --no_antenna_crash --no_pol' , log='h5_merger.log', commandType='python')
        s.run(check=True)
        lib_util.run_losoto(s, f'tec-CS-c{c}', f'{sol_dir}/cal-tec-CS-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
                            plots_dir=f'self/plots/plots-tec-CS-c{c}', h5_dir='self/solutions')
        s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-merged-c{c}.h5 --h5_tables {sol_dir}/cal-tec-RS-c{c}.h5 {sol_dir}/cal-tec-CS-c{c}.h5 \
               --h5_time_freq {sol_dir}/cal-tec2-c{c}.h5 --no_antenna_crash --no_pol' , log='h5_merger.log', commandType='python')
        s.run(check=True)
        lib_util.run_losoto(s, f'tec-merged-c{c}', f'{sol_dir}/cal-tec-merged-c{c}.h5',
                            [f'{parset_dir}/losoto-plot-scalar.parset'], plots_dir=f'self/plots/plots-tec-merged-c{c}',
                            h5_dir='self/solutions')
    facetregname = f'self/solutions/facets-c{c}.reg'
    
    
    with w.if_todo('c%02i_subtract_3Csources' % c):
        for patch in patches:
            if "patch" in patch:
                continue
            
            clean_empty(MSs, "empty-pre-subtract-"+patch, size=imgsizepix_wide, col="CORRECTED_DATA")
            MSs.run(
                f"taql 'UPDATE $pathMS SET CORRECTED_DATA_FR = CORRECTED_DATA_FR - {patch}'",
                log = f'$nameMS_subtract_{patch}.log',
                commandType = 'general'
            )

            MSs.run(
                f"taql 'UPDATE $pathMS SET CORRECTED_DATA = CORRECTED_DATA - {patch}'",
                log = f'$nameMS_subtract_{patch}.log', 
                commandType = 'general'
            )
            clean_empty(MSs, "empty-post-subtract-"+patch, size=imgsizepix_wide, col="CORRECTED_DATA")
            MSs.deletecol(patch)
            
            

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
                                 weight='briggs -0.3', niter=1000000, gridder='wgridder',  parallel_gridding=32, no_update_model_required='', minuv_l=30, mgain=0.9, parallel_deconvolution=1024,
                                 auto_threshold=5.0, auto_mask=8.0, beam_size=15, join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                                 multiscale='', pol='i', nmiter=6,   facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000' )
            # make initial mask
            current_best_mask = make_current_best_mask(imagename, mask_threshold[c], userReg)
            # safe a bit of time by reusing psf and dirty in first iteration
            reuse_kwargs = {'reuse_psf':imagename, 'reuse_dirty':imagename}
        else:
            current_best_mask = f'img/wideM-{c-1}-mask.fits'

        if c>1: # add earlier if bug is fixed
            # if we are done with the things that require blanked pedict, we can also use beam keywords
            beam_kwargs = {'apply_facet_beam':'', 'facet_beam_update':120, 'use_differential_lofar_beam':''}
        else:
            beam_kwargs = {}
        # clean again, with mask now
        logger.info('Making wide field image ...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM,  fits_mask=current_best_mask, data_column='CORRECTED_DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder',  parallel_gridding=32, save_source_list='',
                             update_model_required='', minuv_l=30, beam_size=15, mgain=0.9, nmiter=20, parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=4.0,
                             join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='',  multiscale_scale_bias=0.65, multiscale_max_scales=5, pol='i',
                             facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000 ',
                             **reuse_kwargs, **beam_kwargs)
        if c == 0: # imagename if c=0 else imagenameM
            reuse_kwargs = {'reuse_psf':imagename, 'reuse_dirty':imagename}
        # make a new mask from the image
        current_best_mask = make_current_best_mask(imagenameM, mask_threshold[c]-0.5, userReg)
        # rename source lists
        os.system(f'mv {imagenameM}-sources.txt {imagenameM}-clean1-sources.txt')
        os.system(f'mv {imagenameM}-sources-pb.txt {imagenameM}-clean1-sources-pb.txt')
        # Continue-clean faint sources using local-rms
        logger.info('Making wide field image ...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM, fits_mask=current_best_mask, data_column='CORRECTED_DATA', size=imgsizepix_wide, scale='4arcsec',
                             weight='briggs -0.3', niter=1000000, gridder='wgridder',  parallel_gridding=32, save_source_list='', local_rms='', cont='',
                             update_model_required='', minuv_l=30, beam_size=15, mgain=0.9, nmiter=10, parallel_deconvolution=1024, auto_threshold=2.0, auto_mask=4.0,
                             join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6), deconvolution_channels=3,
                             multiscale='',  multiscale_scale_bias=0.65, multiscale_max_scales=5, pol='i',
                             facet_regions=facetregname, scalar_visibilities='', apply_facet_solutions=f'self/solutions/cal-tec-merged-c{c}.h5 phase000 ', **reuse_kwargs, **beam_kwargs)
        # make a new mask from the image
        current_best_mask = make_current_best_mask(imagenameM, mask_threshold[c]-0.5, userReg)
        # merge source lists
        os.system(f'mv {imagenameM}-sources.txt {imagenameM}-clean2-sources.txt')
        os.system(f'mv {imagenameM}-sources-pb.txt {imagenameM}-clean2-sources-pb.txt')
        os.system(f'cat {imagenameM}-clean1-sources-pb.txt {imagenameM}-clean2-sources-pb.txt >> {imagenameM}-sources-pb.txt')
        os.system(f'cat {imagenameM}-clean1-sources.txt {imagenameM}-clean2-sources.txt >> {imagenameM}-sources.txt')

        # reset NaNs if present
        im = lib_img.Image(f'{imagenameM}-MFS-image.fits')
        im.nantozeroModel()

    ########################### TESTING FOR DD-AMP-CAL ####################################
    # if c > 3:
    #     with w.if_todo('solve_amp_c%02i' % c):
    #         MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
    #         logger.info('Solving amp (slowCS)...')
    #         MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset memory_logging=True numthreads=32 msin=$pathMS sol.nchan=24 \
    #                       sol.mode=scalaramplitude sol.h5parm=$pathMS/amp-dd-c{c}.h5 sol.solint={300 * base_solint} sol.modeldatacolumns="[{",".join(patches)}]"',
    #                 log='$nameMS_solamp-c' + str(c) + '.log', commandType='DP3', maxThreads=1)
    #
    #
    #         lib_util.run_losoto(s, f'amp-dd-c{c}', [ms + f'/amp-dd-c{c}.h5' for ms in MSs.getListStr()],
    #                             [f'{parset_dir}/losoto-plot-scalaramp.parset'],
    #                             plots_dir=f'self/plots/plots-amp-dd-c{c}', h5_dir=f'self/solutions/')
    #
    #     with w.if_todo(f'merge_h5_2_c{c}'):
    #         sol_dir = 'self/solutions'
    #         # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
    #         # lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-amp-dd-c{c}.h5', sourcedb)
    #         # # # reference, unflag and reset the added stations
    #         # h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-amp-merged-c{c}.h5", min_distance=1/3600,
    #         #                    h5_tables=[f"{sol_dir}/cal-amp-dd-c{c}.h5"], propagate_flags=True,
    #         #                    h5_time_freq=f"{sol_dir}/cal-tec-merged-c{c}.h5",ms_files='mss/TC*.MS',no_antenna_crash=True)
    #         # lib_util.run_losoto(s, f'amp-merged-c{c}', f'{sol_dir}/cal-amp-merged-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
    #         #                     plots_dir=f'self/plots/plots-amp-merged-c{c}', h5_dir='self/solutions')
    #         # # remove error000 - this is to avoid a bug in h5parm merger
    #         # s.add(f'losoto -d sol000/error000 self/solutions/cal-tec-merged-c{c}.h5', log='losoto-removeerror.log', commandType="python")
    #         # s.run(maxThreads=1, check=True)  # final check on losoto-final.log
    #         h5_merger.merge_h5(h5_out=f"{sol_dir}/cal-merged-c{c}.h5",
    #                            h5_tables=[f"{sol_dir}/cal-amp-merged-c{c}.h5",f"{sol_dir}/cal-tec-merged-c{c}.h5"], propagate_flags=True,
    #                            h5_time_freq=f"{sol_dir}/cal-tec-merged-c{c}.h5",ms_files='mss/TC*.MS',no_antenna_crash=True)
    #         lib_util.run_losoto(s, f'merged-c{c}', f'{sol_dir}/cal-merged-c{c}.h5', [f'{parset_dir}/losoto-plot-scalar.parset'],
    #                             plots_dir=f'self/plots/plots-merged-c{c}', h5_dir='self/solutions')
    #
    #     # logger.info('Deleting model columns...')
    #     # for patch in patches:
    #     #     MSs.run(f'taql "ALTER TABLE $pathMS DELETE COLUMN {patch}"',
    #     #             log='$nameMS_taql_delcol.log', commandType='general')
    #
    #     imagenameM = 'img/wideM-3ampdd'
    #     facetregname = f'self/solutions/facets-c{c}.reg'
    #     lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), name=imagenameM,
    #                          fits_mask='img/wideM-3-mask.fits', data_column='CORRECTED_DATA', size=imgsizepix_wide,
    #                          scale='4arcsec',
    #                          weight='briggs -0.3', niter=1000000, gridder='wgridder', parallel_gridding=32,
    #                          no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12,
    #                          parallel_deconvolution=1024, auto_threshold=3.0, auto_mask=4.0,
    #                          join_channels='', fit_spectral_pol=3, channels_out=MSs.getChout(4.e6),
    #                          deconvolution_channels=3,
    #                          multiscale='', multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80', pol='i',
    #                          apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='', scalar_visibilities='',
    #                          facet_regions=facetregname, apply_facet_solutions=f'self/solutions/cal-merged-c{c}.h5 phase000,amplitude000 ',
    #                          )
    #     sys.exit()
    #####################################################################################################
    # Find calibration solutions for subfield
    if c < 2:
        # Prepare region and models for subfield
        if subfield:
            subfield_path = subfield
            if len(Regions.read(subfield_path)) > 1:
                raise ValueError(f'Manual subfield region {subfield} contains more than one region.')
        else:
            subfield_path = 'self/skymodel/subfield.reg'

        with w.if_todo('c%02i_extreg_prepare' % c):
            if not subfield and not os.path.exists(subfield_path): # automatically find subfield
                sm = lsmtool.load(f'img/wideM-{c}-sources.txt')
                # sm.remove('img/wide-lr-mask.fits=1')  # remove sidelobe sources that were subtracted
                sm.remove('MajorAxis > 80')  # remove largest scales
                field_center1, field_size1 = lib_dd.make_subfield_region(subfield_path, MSs.getListObj()[0], sm,
                                                                         subfield_min_flux, debug_dir='img/')
            if subtract_predict_mode == 'wsclean':
                # prepare model of central/external regions
                logger.info('Blanking central region of model files and reverse...')
                for im in glob.glob(f'img/wideM-{c}*model*.fits'):
                    wideMint = im.replace('wideM','wideMint')
                    os.system('cp %s %s' % (im, wideMint))
                    lib_img.blank_image_reg(wideMint, subfield_path, blankval = 0., inverse=True)
                    wideMext = im.replace('wideM','wideMext')
                    os.system('cp %s %s' % (im, wideMext))
                    lib_img.blank_image_reg(wideMext, subfield_path, blankval = 0.)
            elif subtract_predict_mode == 'DP3':
                raise ValueError('Not implemented.')
            else: raise ValueError
        # DONE
        subfield_reg = Regions.read(subfield_path)[0]
        field_center = subfield_reg.center.ra, subfield_reg.center.dec
        field_size = np.max([subfield_reg.width.to_value('deg'), subfield_reg.height.to_value('deg')])

        with w.if_todo('c%02i_xtreg_subtract' % c):
            # Recreate MODEL_DATA of external region for subtraction
            MSs.run('taql "update $pathMS set MODEL_DATA=0"', log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
            logger.info('Predict corrupted model of external region (wsclean)...')
            # TODO ignore facet beam for now due to bug
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMext-{c} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                    -facet-regions {facetregname} -apply-facet-solutions self/solutions/cal-tec-merged-c{c}.h5 phase000 \
                    {MSs.getStrWsclean()}',
                log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
            s.run(check=True)

            # # Add model to MODEL_DATA
            # logger.info('Predict corrupted model of external region (DP3)...')
            # MSs.run(f'DP3 {parset_dir}/DP3-h5parmpredict.parset msin=$pathMS pre.sourcedb=$pathMS/wideM-{c}-sources-pb.txt pre.applycal.parmdb=cal-tec-merged-c{c}.h5',
            #         log='$nameMS_h5pre.log', commandType='DP3')

            # cycle > 0: need to add DI-corruption on top (previous iteration sub-field)
            if c > 0:
                logger.info('Add previous iteration sub-field corruption on top of DD-corruption...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c-1) + '.h5 cor.correction=phase000 cor.invert=False',
                        log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
                logger.info('Add DI amplitude corruption on top of DD-corruption...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=self/solutions/cal-amp-di.h5 cor.correction=amplitudeSmooth cor.invert=False',
                        log='$nameMS_sidelobe_corrupt.log', commandType='DP3')

            # subtract external region from CORRECTED_DATA_FR to create SUBFIELD_DATA
            MSs.addcol('SUBFIELD_DATA','CORRECTED_DATA_FR')
            logger.info('Subtracting external region model (SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
            if c > 0:
                logger.info('Correct subfield DI amplitude...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-amp-di.h5 cor.correction=amplitudeSmooth cor.invert=True',
                        log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
            clean_empty(MSs,f'only_subfield-{c}', 'SUBFIELD_DATA') # DEBUG
        ### DONE

        with w.if_todo('c%02i_intreg_predict' % c):
            # Recreate MODEL_DATA of internal region for solve
            logger.info('Predict model of internal region...')
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMint-{c} -j {s.max_processors} -channels-out {MSs.getChout(4.e6)} \
                    {MSs.getStrWsclean()}',
                  log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean', processors='max')
            s.run(check=True)
        ### DONE

        with w.if_todo('c%02i_subfield_solve_tecRS' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info(f'Solving {sf_phaseSolMode} (fast RS)...')
            solve_iono(MSs, c, '2-sf', ['MODEL_DATA'], 1.0, base_solint, sf_phaseSolMode, resetant='intermediate')
        ### DONE

        with w.if_todo('c%02i_subfield_corr_tecRS' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if sf_phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if sf_phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec2-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        with w.if_todo('c%02i_subfield_solve_tecCS' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info(f'Solving {sf_phaseSolMode} (mid RS)...')
            solve_iono(MSs, c, '1-sf', ['MODEL_DATA'], 3.0, 4*base_solint, sf_phaseSolMode, resetant='inner')
        ### DONE

        with w.if_todo('c%02i_subfield_corr_tecCS' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if sf_phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if sf_phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec1-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        with w.if_todo('c%02i_subfield_solve_tecCS0' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info(f'Solving {sf_phaseSolMode} (slowCS)...')
            solve_iono(MSs, c, '0-sf', ['MODEL_DATA'], 6.0, 16*base_solint, sf_phaseSolMode)
        ### DONE

        with w.if_todo('c%02i_subfield_corr_tecCS0' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if sf_phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec0-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if sf_phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec0-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        # merge solutions into one h5parms for plotting and apply
        with w.if_todo('c%02i_merge_h5_subfield' % c):
            sol_dir = 'self/solutions'
            lib_util.check_rm(f'{sol_dir}/cal-tec-sf-merged-c{c}.h5')
            # reference, unflag and reset the added stations
            s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-sf-merged-c{c}.h5 --h5_tables {sol_dir}/cal-tec0-sf-c{c}.h5 {sol_dir}/cal-tec1-sf-c{c}.h5 {sol_dir}/cal-tec2-sf-c{c}.h5 \
                   --h5_time_freq {sol_dir}/cal-tec2-sf-c{c}.h5 --no_antenna_crash --no_pol', log='h5_merger.log',
                commandType='python')
            s.run(check=True)
            lib_util.run_losoto(s, f'tec-sf-merged-c{c}', f'{sol_dir}/cal-tec-sf-merged-c{c}.h5',
                                [f'{parset_dir}/losoto-plot-scalar.parset'],
                                plots_dir=f'self/plots/plots-tec-sf-merged-c{c}', h5_dir='self/solutions')
        ### DONE

        # Do a quick debug image...
        with w.if_todo('c%02i_image-subfield' % c):
            logger.info('Test image subfield...')
            lib_util.run_wsclean(s, 'wscleanSF-c'+str(c)+'.log', MSs.getStrWsclean(), name=f'img/subfield-{c}', data_column='SUBFIELD_DATA', size=3000, scale='4arcsec',
                                 weight='briggs -0.3', niter=100000, gridder='wgridder',  parallel_gridding=6, shift=f'{field_center[0].to(u.hourangle).to_string()} {field_center[1].to_string()}',
                                 no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                                 join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6), deconvolution_channels=3, baseline_averaging='',
                                 multiscale='',  multiscale_scale_bias=0.7, pol='i')
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

                # subtract internal region from CORRECTED_DATA_FR to create SUBFIELD_DATA
                logger.info('Subtract main-lobe (SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA"',
                        log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
                clean_empty(MSs, 'only_sidelobe', 'SUBFIELD_DATA', size=10000)
            ### DONE

            # Do a rough correction of the sidelobe data using the subfield solutions
            with w.if_todo(f'correct-sidelobe-c{c}'): # just for testing/debug
                logger.info('Correct sidelobe data with subfield iono solutions...')
                # merged h5parm is always phase no matter the soltype!
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            # DONE

            imagename_lr = 'img/wide-lr'
            # Image the sidelobe data
            with w.if_todo('image_sidelobe'):
                logger.info('Cleaning low-res...')
                lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=True, data_column='SUBFIELD_DATA',
                                     parallel_gridding=4, size=imgsizepix_lr, scale='30arcsec',
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
                # blank within main-libe to not subtract anything from there (maybe extended sources not properly sobtracted at the beginning)
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
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000 cor.invert=False',
                        log='$nameMS_sidelobe_corrupt.log', commandType='DP3')

                logger.info('Subtract corrupted sidelobe model (CORRECTED_DATA_FR = CORRECTED_DATA_FR - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set CORRECTED_DATA_FR = CORRECTED_DATA_FR - MODEL_DATA"', log='$nameMS_taql-c' + str(c) + '.log', commandType='general')

            if c < 2 :
                # Only after the 0th and 1st iteration: apply the subfield solutions to the data.
                with w.if_todo('c%02i_corr_sf_sols' % c):
                    logger.info('Correct subfield ionosphere (CORRECTED_DATA_FR -> CORRECTED_DATA)...')
                    # Correct MSs:CORRECTED_DATA_FR -> CORRECTED_DATA
                    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA_FR  \
                                cor.parmdb=self/solutions/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000',
                            log='$nameMS_sf-correct.log', commandType='DP3')
                    # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
                    if c == 1:
                        logger.info('Correct DI amplitude (CORRECTED_DATA -> CORRECTED_DATA)...')
                        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                                cor.parmdb=self/solutions/cal-amp-di.h5 cor.correction=amplitudeSmooth cor.invert=True',
                                log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
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
# debugging images -> can be removed in production
[ os.system('mv img/subfield-'+str(c)+'-MFS-image*.fits self/images') for c in range(maxIter) ]
[ os.system('mv img/only*image.fits self/images') for c in range(maxIter) ]
[ os.system('mv img/empty*image.fits self/images') for c in range(maxIter) ]
# os.system('mv img/wideP-MFS-*-image.fits self/images')
# os.system('mv img/wide-lr-MFS-image.fits self/images')

# Copy model
os.system(f'mv img/wideM-{maxIter-1}-*-model.fits self/skymodel')
os.system(f'mv img/wideM-{maxIter-1}-*-model-pb.fits self/skymodel')

logger.info("Done.")
