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

# TODO subtraction of sidelobe for RS smaring
# TODO add timesmearing
# TODO the subfield algorithm should not cut any sources ... how to best implement that? Something with mask islands?

# Waiting for bug fixes in other software
# TODO add LoTSS query for statring model once bug is fixed! (Don't use for now, it crashes the VO server)
# TODO add BDA

import sys, os, glob, random, json
import numpy as np
from regions import Regions
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs
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
remove3c = parset.getboolean('LOFAR_ddparallel', 'remove3c') # get rid of 3c sources in the sidelobes
fulljones = parset.getboolean('LOFAR_ddparallel', 'fulljones') # do fulljones DIE amp correction
min_facets = parset.get('LOFAR_ddparallel', 'min_facets') # ''=default (differs for SPARSE and OUTER), otherwise provide comma seperated list [2,3,6..]
min_flux_factor = parset.getfloat('LOFAR_ddparallel', 'min_flux_factor') # min facet flux factor, default = 1. Higher value -> less facets.
develop = parset.getboolean('LOFAR_ddparallel', 'develop') # for development, make more output/images
sf_phaseSolMode = 'phase' #'tec'
start_sourcedb = parset.get('model','sourcedb')
userReg = parset.get('model','userReg')

#############################################################################

def clean_empty(MSs, name, col='CORRECTED_DATA', size=5000, shift=None):
    """ For testing/debugging only"""
    if shift is None:
        lib_util.run_wsclean(s, 'wsclean-empty.log', MSs.getStrWsclean(), name=f'img/{name}',
                            data_column=col, size=size, scale=f'{int(pixscale*2)}arcsec', niter=1, nmiter=0,
                            weight='briggs 0.0', gridder='wgridder', parallel_gridding=1,
                            no_update_model_required='', apply_primary_beam='')
    else:
        lib_util.run_wsclean(s, 'wsclean-empty.log', MSs.getStrWsclean(), name=f'img/{name}',
                            data_column=col, size=size, scale=f'{int(pixscale*2)}arcsec', niter=1, nmiter=0,
                            weight='briggs 0.0', gridder='wgridder', parallel_gridding=1, shift= f'{shift[0]} {shift[1]}',
                            no_update_model_required='', apply_primary_beam='')

def corrupt_model_dirs(MSs, c, tc, model_columns, solmode='phase'):
    """ CORRUPT the MODEL_DATA columns for model_columns
    Parameters
    ----------
    MSs: MSs object
    c: int, cycle
    tc: tec/phase step suffix
    model_columns: list of patch names
    solmode: which solution mode to corrupt for (phase, tec, tecandphase)
    """
    logger.info(f'Corrupt models: {solmode}{tc}...')
    for model_column in model_columns:
        if solmode in ['tec', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}]  \
                    cor.parmdb={sol_dir}/cal-tec{tc}-c{c}.h5 cor.correction=tec000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        elif solmode in ['phase', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                    cor.parmdb={sol_dir}/cal-tec{tc}-c{c}.h5 cor.correction=phase000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        elif solmode in ['amplitude']:
            if not model_column.startswith('patch'):    
                MSs.run(
                    f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                        cor.parmdb={sol_dir}/cal-amp-3C.h5 cor.correction=amplitude000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')

def solve_iono(MSs, c, tc, model_columns, smMHz, solint, solmode, resetant=None, constrainant=None, model_column_fluxes=None, variable_solint=False):
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
    variable_solint: bool, optional. Default = False. Use twice as long solints for patches < 8 Jy and for times as long for patches <3 Jy
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
            losoto_parsets = [parset_dir+'/losoto-refph.parset', resetant_parset, parset_dir+'/losoto-plot-scalarph.parset']
        else:
            losoto_parsets = [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-plot-scalarph.parset']
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
    if variable_solint and solmode == 'tec':
        raise ValueError
    elif variable_solint: # if actived, use twice the solint for fainter directions
        solint *= 4 # use twice as long solint
        # get two solutions per solint (i.e. one per time step) for bright directions
        solutions_per_direction[model_column_fluxes > 4] = 2
        solutions_per_direction[model_column_fluxes > 8] = 4
    
    maxThreads = 1 if len(model_columns) > 30 else None
    if solmode == 'phase':
        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec{tc}.h5 sol.solint={solint} \
                  sol.mode=scalarphase sol.smoothnessconstraint={smMHz}e6 sol.smoothnessreffrequency=54e6 sol.nchan=1 {antennaconstraint} \
                  sol.modeldatacolumns="[{",".join(model_columns)}]" sol.solutions_per_direction="{np.array2string(solutions_per_direction, separator=",")}"',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3', maxProcs=maxThreads)
    else:
        MSs.run(f'DP3 {parset_dir}/DP3-solTEC.parset msin=$pathMS sol.h5parm=$pathMS/tec{tc}.h5 sol.solint={solint} {antennaconstraint} \
                  sol.modeldatacolumns="[{",".join(model_columns)}]" sol.mode={solmode}',
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3', maxProcs=maxThreads)

    lib_util.run_losoto(s, f'tec{tc}-c{c}', [ms+f'/tec{tc}.h5' for ms in MSs.getListStr()], losoto_parsets, 
                        plots_dir=f'{plot_dir}/plots-tec{tc}-c{c}', h5_dir=sol_dir)
    

def corr_die_amp(MSs, col='CORRECTED_DATA', fulljones=fulljones, invert=True):
    """
    invert: set to False to corrupt the chosen column
    """
    if fulljones:
        logger.info('Correct amp-di (fulljones)...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={col} msout.datacolumn={col} \
                cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=fulljones cor.soltab=[amplitudeSmooth,phaseSmooth] \
                cor.updateweights=False cor.invert={invert}',
                log='$nameMS_diampcor.log', commandType='DP3')
                
        logger.info('Correct amp-di normalisation...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={col} msout.datacolumn={col} \
                cor.parmdb={sol_dir}/cal-amp-dinorm.h5 cor.correction=amplitude000 \
                cor.updateweights=False cor.invert={invert}',
                log='$nameMS_diampcor.log', commandType='DP3')
    else:
        logger.info('Correct amp-di...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={col} msout.datacolumn={col} \
                cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=amplitudeSmooth \
                cor.updateweights=False cor.invert={invert}',
                log='$nameMS_diampcor.log', commandType='DP3')
        

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

def add_3c_models(sm, phasecentre, null_mid_freq, max_sep=50., threshold=0):
    
    with open(parset_dir+"/3C_coordinates.json", "r") as file:
        all_3c = json.load(file)
    
    phasecentre = SkyCoord(phasecentre[0], phasecentre[1], unit=(u.deg, u.deg))
    beam_hdu = fits.open(beamMask)[0]
    beam_wcs = wcs.WCS(beam_hdu.header)
        
    logger.info('Adding 3C models...')
    for source, coord in all_3c.items():
        
        pos = SkyCoord(ra=coord[0], dec=coord[1], unit=(u.hourangle, u.deg))
        sep = phasecentre.separation(pos).deg
        
        if sep < 20:
            pix_pos = np.round(np.asarray(wcs.utils.skycoord_to_pixel(pos, beam_wcs, origin=1))).astype(int)
            if (0 <= pix_pos[0] < beam_hdu.data.shape[-2]) and (0 <= pix_pos[1] < beam_hdu.data.shape[-1]):
                within_beam = bool(beam_hdu.data[0,0,pix_pos[1], pix_pos[0]])
        elif sep > 20 and source == "3C 274":
            logger.info(f'3C source {source} (seperation: {sep:.2f} deg) too far from center to reliably subtract: ignore.')
            continue
        else:
            within_beam = False
            
        if within_beam:
            logger.info(f'3C source {source} is within primary beam. Not Adding model for subtraction.')
            continue
        elif phasecentre.separation(pos).deg > max_sep:
            continue
        
        if threshold == 0:
            #determining a linear threshold based on the distance from the null and beam corrected flux
            hnmf = null_mid_freq/2
            a = (19 - 1)/(50 - hnmf)
            threshold = a * (sep - hnmf) + 1
        
        if source in ["3C 196", "3C 380", "3C 295"]: # take pre-existing model for calibrators
            sourcedb = os.path.dirname(__file__) + f'/../models/calib-simple.skymodel'
            sm_3c = lsmtool.load(sourcedb, beamMS=sm.beamMS)
            sm_3c.select(f'patch=={source.replace(" ","")}')
        else:
            sourcedb = os.path.dirname(__file__) + f'/../models/3CRR/{source.replace(" ","")}.txt'
            if not os.path.exists(sourcedb):
                logger.warning(f'No model found for {source} (seperation {sep:.2f} deg)')
                continue
            sm_3c = lsmtool.load(sourcedb, beamMS=sm.beamMS)
        
        sm_3c.setColValues("Patch", ["source_"+source.replace(" ","")]*len(sm_3c.getColValues("I")))
        flux_3c =  sm_3c.getColValues("I", aggregate="sum", applyBeam=True)[0]
        if flux_3c > threshold:
            sm_3c.setColValues("Patch", ["source_"+source.replace(" ","")]*len(sm_3c.getColValues("I")))
            logger.info(f'3C source {source} (seperation: {sep:.2f} deg) app. flux {flux_3c:.2f} Jy is above threshold {threshold:.2f} Jy: keep.')
            sm.concatenate(sm_3c)
            sm.setPatchPositions(method='wmean', applyBeam=True)
        else:
            logger.debug(f'3C source {source} (seperation: {sep:.2f} deg) app. flux {flux_3c:.2f} Jy is below threshold {threshold:.2f} Jy: ignore.')
        
        del sm_3c
        threshold = 0
    return sm

def make_source_regions(sm, c):
    lib_util.check_rm(f'ddparallel/skymodel/regions_c{c}')
    lib_util.check_rm(f'ddparallel/skymodel/sources_c{c}.reg')
    os.makedirs(f'ddparallel/skymodel/regions_c{c}')
    for p in sm.getPatchNames():
        sm_p = sm.copy()
        sm_p.select(f'patch=={p}')
        sm_p.write(f'ddparallel/skymodel/regions_c{c}/{p}.reg', format='ds9', clobber=True)
        regs = Regions.read(f'ddparallel/skymodel/regions_c{c}/{p}.reg')
        col = '#{:06x}'.format(random.randint(0, 256 ** 3))
        for reg in regs:
            reg.visual['facecolor'] = col
            reg.visual['edgecolor'] = col
        regs.write(f'ddparallel/skymodel/regions_c{c}/{p}.reg',overwrite=True)
        os.system(f'cat ddparallel/skymodel/regions_c{c}/*.reg >> ddparallel/skymodel/sources_c{c}.reg')
    

#############################################################################
# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')

    # here images, models, solutions for each group will be saved
    lib_util.check_rm('ddparallel')
    if not os.path.exists('ddparallel/images'): os.makedirs('ddparallel/images')
    if not os.path.exists('ddparallel/solutions'): os.makedirs('ddparallel/solutions')
    if not os.path.exists('ddparallel/plots'): os.makedirs('ddparallel/plots')
    if not os.path.exists('ddparallel/skymodel'): os.makedirs('ddparallel/skymodel')
### DONE

sol_dir = 'ddparallel/solutions'
plot_dir = 'ddparallel/plots'

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_flags=True)
MSs.print_HAcov()

# make beam to the first mid null - outside of that do a rough subtraction and/or 3C peeling. Use sources inside for calibration
phasecentre = MSs.getListObj()[0].getPhaseCentre()
null_mid_freq = max(MSs.getListObj()[0].getFWHM(freq='mid', elliptical=True)) * 1.8 # FWHM to null

# set image size - this should be a bit more than the beam region used for calibration
pixscale = MSs.getListObj()[0].getPixelScale()
imgsizepix_wide = int(1.85*max(MSs.getListObj()[0].getFWHM(freq='mid', elliptical=True))*3600/pixscale) # roughly to null
if imgsizepix_wide > 10000: imgsizepix_wide = 10000
if imgsizepix_wide % 2 != 0: imgsizepix_wide += 1  # prevent odd img sizes
imgsizepix_lr = int(5*max(MSs.getListObj()[0].getFWHM(freq='mid', elliptical=True))*3600/(pixscale*8))
if imgsizepix_lr % 2 != 0: imgsizepix_lr += 1  # prevent odd img sizes

logger.info(f'Setting wide-field image size: {imgsizepix_wide}pix; scale:  {pixscale:.2f}arcsec.')

# set clean componet fit order (use 5 for large BW)
if MSs.getChout(4.e6) >= 7:  # Bandwidth of 28 MHz or more
    cc_fit_order = 5
else: cc_fit_order = 3

fullband = MSs.getBandwidth()
nchan = MSs.mssListObj[0].getNchan()
tint = MSs.mssListObj[0].getTimeInt()
if int(np.rint(fullband / nchan < 195.3e3/4)):
    base_nchan = int(np.rint((195.3e3/4)/(fullband/nchan))) # this is 1 for dutch observations, and larger (2,4) for IS observations
else: base_nchan = 1
if tint < 4:
    base_solint = int(np.rint(4/tint)) # this is already 4 for dutch observations
else: base_solint = 1

mask_threshold = [5.0,4.5,4.0,4.0,4.0,4.0] # sigma values for beizorro mask in cycle c

# define list of facet fluxes per iteration -> this can go into the config
if 'OUTER' in MSs.getListObj()[0].getAntennaSet():
    facet_fluxes = min_flux_factor*np.array([4,2.0, 1.2, 1.0, 0.9, 0.8])*(54e6/np.mean(MSs.getFreqs()))**0.7 # this is not the total flux, but the flux of bright sources used to construct the facets. still needs to be tuned, maybe also depends on the field
elif 'SPARSE' in MSs.getListObj()[0].getAntennaSet():
    facet_fluxes = min_flux_factor*np.array([4,2.4, 1.3, 1.1, 1.0, 0.9])*(54e6/np.mean(MSs.getFreqs()))**0.7 # this is not the total flux, but the flux of bright sources used to construct the facets. still needs to be tuned, maybe also depends on the field

if min_facets: # if manually provided
    if not isinstance(min_facets, list):
        min_facets = min_facets.replace('[', '').replace(']', '').split(',')
        min_facets = np.array(min_facets).astype(int)
else: #default settings
    # use more facets for SPARSE (larger FoV)
    if 'SPARSE' in MSs.getListObj()[0].getAntennaSet():
        min_facets = [3, 6, 18, 24, 24, 24]
    elif 'OUTER' in MSs.getListObj()[0].getAntennaSet():
        min_facets = [2, 4, 12, 20, 20, 20]
    else:
        raise ValueError(f'{MSs.getListObj()[0].getAntennaSet()} not recognized.')

smMHz2 = [2.0,5.0,5.0,5.0,5.0,5.0]
smMHz1 = [8.0,12.0,12.0,12.0,12.0,12.0]
# smMHz0 = [6.0,10.0,10.0,10.0,10.0,10.0]

# Make beam mask/reg
beamMask = 'ddparallel/beam.fits'
beamReg = 'ddparallel/beam.reg'
MSs.getListObj()[0].makeBeamReg(beamReg, freq='mid', to_pbval=0)
if not os.path.exists(beamMask):
    logger.info('Making mask of primary beam...')
    lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name=beamMask.replace('.fits',''), size=imgsizepix_lr, scale='30arcsec')
    os.system(f'mv {beamMask.replace(".fits","-image.fits")} {beamMask}') # beam-image.fits -> beam.fits
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 1.)
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 0., inverse=True)
subfield_path = 'ddparallel/skymodel/subfield.reg'

#################################################################################################

with w.if_todo('solve_fr'):
    # Probably we do not need BLsmoothing since we have long time intervals and smoothnessconstraint?
    logger.info('Converting to circular (DATA -> CIRC_PHASEDIFF_DATA)...')
    # lin2circ conversion TCxx.MS:DATA -> CIRC_PHASEDIFF_DATA # use no dysco here!
    MSs.run(f'DP3 {parset_dir}/DP3-lin2circ.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CIRC_PHASEDIFF_DATA', log='$nameMS_lin2circ.log', commandType="DP3")

    # Get circular phase diff TC.MS: CIRC_PHASEDIFF_DATA -> CIRC_PHASEDIFF_DATA
    logger.info('Get circular phase difference...')
    MSs.run('taql "UPDATE $pathMS SET\
         CIRC_PHASEDIFF_DATA[,0]=0.5*EXP(1.0i*(PHASE(CIRC_PHASEDIFF_DATA[,0])-PHASE(CIRC_PHASEDIFF_DATA[,3]))), \
         CIRC_PHASEDIFF_DATA[,3]=CIRC_PHASEDIFF_DATA[,0], \
         CIRC_PHASEDIFF_DATA[,1]=0+0i, \
         CIRC_PHASEDIFF_DATA[,2]=0+0i"', log='$nameMS_taql.log', commandType='general')

    logger.info('Creating MODEL_DATA_FR...')  # take from MODEL_DATA but overwrite
    MSs.addcol('MODEL_DATA_FR', 'DATA', usedysco=False)
    MSs.run('taql "UPDATE $pathMS SET MODEL_DATA_FR[,0]=0.5+0i, MODEL_DATA_FR[,1]=0.0+0i, MODEL_DATA_FR[,2]=0.0+0i, \
         MODEL_DATA_FR[,3]=0.5+0i"', log='$nameMS_taql.log', commandType='general')

    # Solve TC.MS: CIRC_PHASEDIFF_DATA against MODEL_DATA_FR (only solve - solint=2m nchan=0 as it has the smoothnessconstrain)
    logger.info('Solving circ phase difference ...')
    MSs.run('DP3 ' + parset_dir + '/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.solint=' + str(30 * base_solint),
            log='$nameMS_solFR.log', commandType="DP3")
    lib_util.run_losoto(s, f'fr', [ms + '/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
    os.system(f'mv cal-fr.h5 {sol_dir}')
    os.system(f'mv plots-fr {plot_dir}')

    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, MODEL_DATA_FR"',
            log='$nameMS_taql.log', commandType='general')

with w.if_todo('cor_fr'):
    # Correct FR with results of solve - TC.MS: DATA -> CORRECTED_DATA_FR
    logger.info('Correcting FR (DATA -> CORRECTED_DATA_FR...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA_FR \
            cor.parmdb={sol_dir}/cal-fr.h5 cor.correction=rotationmeasure000',
            log='$nameMS_corFR.log', commandType='DP3')
    
    # Copy TC.MS: CORRECTED_DATA_FR -> CORRECTED_DATA
    logger.info('Creating CORRECTED_DATA = CORRECTED_DATA_FR...')
    MSs.addcol('CORRECTED_DATA', 'CORRECTED_DATA_FR')
### DONE

#####################################################################################################
# Self-cal cycle
for c in range(maxIter):
    logger.info('Start selfcal cycle: '+str(c))

    # get sourcedb
    sourcedb = f'ddparallel/skymodel/tgts-c{c}.skymodel'
    beamMS = MSs.getListStr()[int(len(MSs.getListStr()) / 2)] # use central MS, should not make a big difference
    if not os.path.exists(sourcedb):
        logger.info(f'Creating skymodel {sourcedb}...')
        if c == 0:
            # if provided, use manual model
            if start_sourcedb:
                logger.info(f'Using input skymodel {start_sourcedb}')
                sm = lsmtool.load(start_sourcedb, beamMS=beamMS)
            else:
                # Create initial sourcedb from GSM
                logger.info('Get skymodel from GSM.')
                s.add(f'wget -O {sourcedb} "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord={phasecentre[0]},{phasecentre[1]}&radius={null_mid_freq/2}&unit=deg"',
                      log='wget.log', commandType='general')
                s.run(check=True)
                sm = lsmtool.load(sourcedb, beamMS=beamMS)
        else:
            # get wsclean skymodel of previous iteration
            wsc_src = f'img/wideM-{c-1}-sources-pb.txt'
            sm = lsmtool.load(wsc_src, beamMS=beamMS)
        
        # if using e.g. LoTSS, adjust for the frequency
        logger.debug(f'Extrapolating input skymodel fluxes from {sm.getDefaultValues()["ReferenceFrequency"]/1e6:.0f}MHz to {np.mean(MSs.getFreqs())/1e6:.0f}MHz assuming si=-0.7')
        si_factor = (np.mean(MSs.getFreqs())/sm.getDefaultValues()['ReferenceFrequency'])**0.7 # S144 = si_factor * S54
        sm.select(f'{beamMask}==True')  # remove outside of FoV (should be subtracted (c>0) or not present (c==0)!)
        sm.group('threshold', FWHM=5/60, root='Src') # group nearby components to single source patch
        sm.setPatchPositions(method='wmean', applyBeam=True)
        sm = lib_dd_parallel.merge_nearby_bright_facets(sm, 1/60, 0.5, applyBeam=True)
        # TODO we need some logic here to avoid picking up very extended sources. Also case no bright sources in a field.
        patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=True)
        if sum(patch_fluxes/si_factor > facet_fluxes[c]) < min_facets[c]:
            bright_sources_flux = np.sort(patch_fluxes)[-min_facets[c]] / si_factor
            logger.warning(f'Less than {min_facets[c]} bright sources above minimum flux {facet_fluxes[c]:.2f} Jy! Using sources above {bright_sources_flux:.2f} Jy')
        else:
            bright_sources_flux = facet_fluxes[c]
        bright_names = sm.getPatchNames()[patch_fluxes > bright_sources_flux*si_factor]
        bright_pos = sm.getPatchPositions(bright_names)
        sm.group('voronoi', targetFlux=bright_sources_flux*si_factor, applyBeam=True, root='', byPatch=True)
        sm.setPatchPositions(bright_pos)
        lib_dd_parallel.rename_skymodel_patches(sm, applyBeam=True)

        if c == 0 and remove3c:
            # Add models of bright 3c sources to the sky model. Model will be subtracted from data before imaging.
            sm = add_3c_models(sm, phasecentre=phasecentre, null_mid_freq=null_mid_freq)
            sm.setColValues("Q", np.zeros(len(sm.getColValues("I")))) # force non I Stokes to zero
            sm.setColValues("U", np.zeros(len(sm.getColValues("I"))))
            sm.setColValues("V", np.zeros(len(sm.getColValues("I"))))
        
        make_source_regions(sm, c)
        sm.write(sourcedb, clobber=True)
        logger.info(f'Using {len(sm.getPatchNames())} patches.')
    else:
        logger.info(f'Load existing skymodel {sourcedb}')
        sm = lsmtool.load(sourcedb, beamMS=beamMS)
    
    sm.plot(f'ddparallel/skymodel/patches-c{c}.png', 'patch')

    # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
    sourcedb_basename = sourcedb.split('/')[-1]
    for MS in MSs.getListStr():
        lib_util.check_rm(MS + '/' + sourcedb_basename)
        logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
        os.system('cp -r ' + sourcedb + ' ' + MS)

    patches = sm.getPatchNames()
    patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=True)
    for patch in patches[np.argsort(patch_fluxes)[::-1]]:
        logger.info(f'{patch}: {patch_fluxes[patches==patch][0]:.1f} Jy')
  
    with w.if_todo('c%02i_init_model' % c):
        for patch in patches:
            # Add model to MODEL_DATA
            # TODO: add time smearing in the predict parset
            logger.info(f'Add model to {patch}...')
            MSs.run(f'DP3 {parset_dir}/DP3-predict-beam.parset msin=$pathMS pre.sourcedb=$pathMS/{sourcedb_basename} pre.sources={patch} msout.datacolumn={patch}',
                    log='$nameMS_pre.log', commandType='DP3')
            # pos = sm.getPatchPositions()[patch]
            # size = int((1.1*sm.getPatchSizes()[np.argwhere(sm.getPatchNames()==patch)]) // 4)
            # logger.info(f'Test image MODEL_DATA...')
            # lib_util.run_wsclean(s, 'wscleanMODEL-c' + str(c) + '.log', MSs.getStrWsclean(), name=f'ddparallel/skymodel/{patch}',
            #                      data_column=patch, size=size, scale='8arcsec', shift=f'{pos[0].to(u.hourangle).to_string()} {pos[1].to_string()}',
            #                      weight='briggs -0.5', niter=10000, gridder='wgridder', parallel_gridding=6, no_update_model_required='', minuv_l=30, mgain=0.9,
            #                      parallel_deconvolution=512, beam_size=15, join_channels='', fit_spectral_pol=3,
            #                      channels_out=MSs.getChout(4.e6), deconvolution_channels=3, pol='i', nmiter=5 )
        #MSs = lib_ms.AllMSs(glob.glob("*-smearing.MS"), s, check_flags=True)
    ### DONE

    ### solve ionosphere phase - ms:SMOOTHED_DATA - > reset for all BUT most distant RS!
    with w.if_todo('c%02i_solve_tecRS' % c):
        # Smooth MSs:CORRECTED_DATA -> SMOOTHED_DATA
        MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
        logger.info('Solving TEC (RS)...')
        solve_iono(MSs, c, '-RS', patches, smMHz2[c], 2*base_solint, 'phase', resetant='CS', model_column_fluxes=patch_fluxes, variable_solint=True)
    ### DONE

    ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('c%02i_corrupt_tecRS' % c):
        corrupt_model_dirs(MSs, c, '-RS', patches)
    ### DONE

    ### solve ionosphere phase - ms:SMOOTHED_DATA - > reset for central CS
    with w.if_todo('c%02i_solve_tecCS' % c):
        logger.info('Solving TEC (CS)...')
        solve_iono(MSs, c, '-CS', patches, smMHz1[c], 30*base_solint, 'phase', constrainant=None,  model_column_fluxes=patch_fluxes, variable_solint=True) # 'RS'
    ### DONE

    ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('c%02i_corrupt_tecCS' % c):
        corrupt_model_dirs(MSs, c, '-CS', patches, 'phase')
    ### DONE

    ########################### 3C-subtract PART ####################################
    ### CORRUPT the Amplitude of MODEL_DATA columns for all 3CRR patches
    if c == 0 and remove3c:
        _3c_patches = [p for p in patches if not p.startswith('patch_')]
        if len(_3c_patches) > 0:
            with w.if_todo('3c_solve_amp'):
                logger.info('Solving amplitude for 3C...')
                # Solve diagonal amplitude MSs:SMOOTHED_DATA
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.model_weighted_constraints=true\
                          sol.mode=diagonalamplitude sol.nchan=1 sol.smoothnessconstraint=4e6 sol.smoothnessreffrequency=54e6 sol.h5parm=$pathMS/amp-3C.h5 sol.datause=full \
                          sol.modeldatacolumns="[{",".join(patches)}]" sol.solint=60', log=f'$nameMS_solamp_3c_c{c}.log', commandType="DP3")

                losoto_parsets = [parset_dir + '/losoto-clip.parset', parset_dir + '/losoto-plot-amp.parset']
                lib_util.run_losoto(s, f'amp-3C', [ms + f'/amp-3C.h5' for ms in MSs.getListStr()], losoto_parsets,
                                    plots_dir=f'{plot_dir}/plots-amp-3C', h5_dir=sol_dir)
                ### DONE

            with w.if_todo('3C_corrupt_subtract'):
                MSs.run('addcol2ms.py -m $pathMS -c FLAG_BKP -i FLAG', log='$nameMS_addcol.log', commandType='python')
                MSs.run('taql "update $pathMS set FLAG_BKP = FLAG"', log='$nameMS_taql.log', commandType='general')
                
                debug_3c_sub = False
                if debug_3c_sub:
                    with open(parset_dir+"/3C_coordinates.json", "r") as file:
                        import json
                        all_3c = json.load(file)
                
                for patch in _3c_patches:
                    logger.info(f'Subtracting {patch}...')
                    # Corrupt MODEL_DATA with amplitude, set MODEL_DATA = 0 where data are flagged, then unflag everything
                    if debug_3c_sub:
                        coords = all_3c[patch.replace('source_','').replace("C","C ")]
                        coords[0] = coords[0].replace(" ", "h", 1).replace(" ", "m", 1) + "s"
                        coords[1] = coords[1].replace(" ", "d", 1).replace(" ", "m", 1) + "s"
                        clean_empty(MSs,f'{patch}_data_zoom_pre', 'CORRECTED_DATA_FR', shift=coords, size=2000)
                        clean_empty(MSs,f'{patch}_model', f'{patch}', shift=coords, size=2000)
                        
                        # TEST: comment out amps
                        corrupt_model_dirs(MSs, c, 1, [patch], solmode='amplitude')
                        MSs.run(f'taql "update $pathMS set {patch}[FLAG] = 0"', log='$nameMS_taql.log', commandType='general')
                        clean_empty(MSs,f'{patch}_model_amp', f'{patch}', shift=coords, size=2000)
                    
                    MSs.run(
                        f"taql 'UPDATE $pathMS SET CORRECTED_DATA_FR = CORRECTED_DATA_FR - {patch}'",
                        log = f'$nameMS_subtract_{patch}.log', 
                        commandType = 'general'
                    )
                    if debug_3c_sub: clean_empty(MSs,f'{patch}_data_zoom_after', 'CORRECTED_DATA_FR', shift=coords, size=2000)

                    MSs.run(
                        f"taql 'UPDATE $pathMS SET CORRECTED_DATA = CORRECTED_DATA - {patch}'",
                        log = f'$nameMS_subtract_{patch}.log', 
                        commandType = 'general'
                    )
                    
                    MSs.deletecol(patch)
                    sm = lsmtool.load(sourcedb, beamMS=beamMS)
                    sm.select(f'Patch != {patch}')
                    sm.write(sourcedb, clobber=True)
                    patches = np.delete(patches, np.where(patches == patch)[0])
                    MSs.run(f'taql "update $pathMS set FLAG = FLAG_BKP"', log='$nameMS_taql.log', commandType='general')

                MSs.deletecol('FLAG_BKP')  
            ### DONE

    ########################### AMP-CAL PART ####################################
    # Only once in cycle 1: do di amp to capture element beam 2nd order effect
    # TODO add updateweights in production
    if c == 1:
        with w.if_todo('amp_di_solve'):
            # no smoothing, CORRECTED_DATA is smoothed above
            logger.info('Setting MODEL_DATA to sum of corrupted patch models...')
            MSs.run(f'taql "UPDATE $pathMS SET MODEL_DATA={"+".join(patches)}"', log='$nameMS_taql_phdiff.log', commandType='general')
            
            if fulljones:
                logger.info('Solving amp-di (fulljones)...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=full sol.nchan=12 sol.modeldatacolumns=[MODEL_DATA] \
                     sol.mode=fulljones sol.h5parm=$pathMS/amp-di.h5 sol.solint={150*base_solint} sol.minvisratio=0.5 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]',
                     log='$nameMS_diampsol.log', commandType='DP3')

                lib_util.run_losoto(s, f'amp-di', [ms + f'/amp-di.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-fj.parset', f'{parset_dir}/losoto-amp-difj.parset'],
                                plots_dir=f'{plot_dir}/plots-amp-di', h5_dir=sol_dir)
            
                # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
                logger.info('Correct amp-di (CORRECTED_DATA -> CORRECTED_DATA)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=fulljones cor.soltab=[amplitudeSmooth,phaseSmooth] \
                    cor.updateweights=False',
                    log='$nameMS_diampcor.log', commandType='DP3')
                
                # TODO: add a contrant to all antennas
                #sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]',
                logger.info('Solving amp-di for normalisation...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.datause=full sol.nchan=0 sol.solint=0 \
                    sol.modeldatacolumns=[MODEL_DATA] sol.mode=diagonal sol.h5parm=$pathMS/amp-dinorm.h5',
                    log='$nameMS_diampsol.log', commandType='DP3')
                
                lib_util.run_losoto(s, f'amp-dinorm', [ms + f'/amp-dinorm.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-amp.parset'], plots_dir=f'{plot_dir}/plots-amp-di', h5_dir=sol_dir)
            
                logger.info('Correct amp-di normalisation (CORRECTED_DATA -> CORRECTED_DATA)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb={sol_dir}/cal-amp-dinorm.h5 cor.correction=amplitude000 cor.updateweights=False',
                    log='$nameMS_diampcor.log', commandType='DP3')

            else:
                logger.info('Solving amp-di...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=full sol.nchan=12 sol.modeldatacolumns=[MODEL_DATA] \
                     sol.mode=diagonal sol.h5parm=$pathMS/amp-di.h5 sol.solint={150*base_solint} sol.minvisratio=0.5 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]',
                     log='$nameMS_diampsol.log', commandType='DP3')

                lib_util.run_losoto(s, f'amp-di', [ms + f'/amp-di.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-amp.parset', f'{parset_dir}/losoto-plot-ph.parset', f'{parset_dir}/losoto-amp-di.parset'],
                                plots_dir=f'{plot_dir}/plots-amp-di', h5_dir=sol_dir)

                # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
                logger.info('Correct amp-di (CORRECTED_DATA -> CORRECTED_DATA)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=amplitudeSmooth cor.updateweights=False',
                    log='$nameMS_diampcor.log', commandType='DP3')
    
    # just correct for it after 1st cycle (not used)
    elif c >= 2:
        corr_die_amp(MSs, col='CORRECTED_DATA', fulljones=fulljones)

    #if c > 1:
    #    # Starting from cycle 2: do DD-slow amp solutions
    #    with w.if_todo(f'c%02i_solve_amp_dd' % c):
    #        # no smoothing, cycle 2 CORRECTED_DATA should be the same
    #        # MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
    #        solutions_per_direction = np.ones(len(patches), dtype=int)
    #        solutions_per_direction[patch_fluxes > 4] = 2
    #        solutions_per_direction[patch_fluxes > 8] = 4
    #        logger.info('Solving amp (slow)...')
    #        # this might be very memory hungry since we solve on long timescales
    #        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset memory_logging=True msin=$pathMS sol.nchan=24 \
    #                      sol.mode=diagonal sol.datause=full sol.h5parm=$pathMS/amp-dd-c{c}.h5 sol.solint={300*base_solint} \
    #                      sol.modeldatacolumns="[{",".join(patches)}]" sol.solutions_per_direction="{np.array2string(solutions_per_direction, separator=",")}"',
    #                log='$nameMS_solamp-c' + str(c) + '.log', commandType='DP3', maxProcs=1)
    #
    #        lib_util.run_losoto(s, f'amp-dd-c{c}', [ms + f'/amp-dd-c{c}.h5' for ms in MSs.getListStr()],
    #                            [f'{parset_dir}/losoto-norm.parset',f'{parset_dir}/losoto-plot-scalaramp.parset'],
    #                            plots_dir=f'{plot_dir}/plots-amp-dd-c{c}', h5_dir=sol_dir)

    # merge solutions into one h5parms for large scale image
    with w.if_todo('c%02i_merge_h5' % c):
        lib_util.check_rm(f'{sol_dir}/cal-tec-merged-c{c}.h5')
        # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec-CS-c{c}.h5', sourcedb)
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec-RS-c{c}.h5', sourcedb)

        # by construction, the 3c sources are always at the last indices
        filter_directions = f"--filter_directions '["+', '.join([str(i) for i in range(len(patches))])+"]'" if (c == 0 and remove3c) else ''
        # if we have dd-amplitudes, add polarization
        usepol = ' --no_pol ' if c <2 else ''
        # reference, unflag and reset the added stations f'{sol_dir}/cal-tec0-c{c}.h5',
        logger.info('Merge solutions...')
        s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-RSm-c{c}.h5 --h5_tables {sol_dir}/cal-tec-RS-c{c}.h5 --h5_time_freq {sol_dir}/cal-tec-RS-c{c}.h5 \
              --no_antenna_crash {usepol} --propagate_flags {filter_directions}' , log='h5_merger.log', commandType='python')
        s.run(check=True)
        lib_util.run_losoto(s, f'tec-RSm-c{c}', f'{sol_dir}/cal-tec-RSm-c{c}.h5', [f'{parset_dir}/losoto-plot-scalarph.parset'],
                            plots_dir=f'{plot_dir}/plots-tec-RSm-c{c}', h5_dir=sol_dir)
        s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-CSm-c{c}.h5 --h5_tables {sol_dir}/cal-tec-CS-c{c}.h5 --h5_time_freq {sol_dir}/cal-tec-RS-c{c}.h5 \
              --no_antenna_crash {usepol} --propagate_flags {filter_directions}' , log='h5_merger.log', commandType='python')
        s.run(check=True)
        lib_util.run_losoto(s, f'tec-CSm-c{c}', f'{sol_dir}/cal-tec-CSm-c{c}.h5', [f'{parset_dir}/losoto-plot-scalarph.parset'],
                            plots_dir=f'{plot_dir}/plots-tec-CSm-c{c}', h5_dir=sol_dir)
        if c > 1:
            lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-amp-dd-c{c}.h5', sourcedb)
            s.add(f'h5_merger.py --h5_out {sol_dir}/cal-amp-dd-merged-c{c}.h5 --h5_tables {sol_dir}/cal-amp-dd-c{c}.h5 --h5_time_freq {sol_dir}/cal-tec-RS-c{c}.h5 \
              --no_antenna_crash {usepol} --propagate_flags {filter_directions}' , log='h5_merger.log', commandType='python')
            s.run(check=True)
            lib_util.run_losoto(s, f'amp-dd-merged-c{c}', f'{sol_dir}/cal-amp-dd-merged-c{c}.h5', [f'{parset_dir}/losoto-plot-scalaramp.parset'],
                                plots_dir=f'{plot_dir}/plots-amp-dd-c{c}', h5_dir=sol_dir)
            s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-merged-c{c}.h5 --h5_tables {sol_dir}/cal-tec-RSm-c{c}.h5 {sol_dir}/cal-tec-CSm-c{c}.h5 {sol_dir}/cal-amp-dd-merged-c{c}.h5 \
                   --h5_time_freq {sol_dir}/cal-tec-RS-c{c}.h5 --no_antenna_crash {usepol} --propagate_flags' , log='h5_merger.log', commandType='python')
            s.run(check=True)
            lib_util.run_losoto(s, f'tec-merged-c{c}', f'{sol_dir}/cal-tec-merged-c{c}.h5',
                                [f'{parset_dir}/losoto-plot-scalarph.parset', f'{parset_dir}/losoto-plot-scalaramp.parset'], plots_dir=f'{plot_dir}/plots-tec-merged-c{c}',
                                h5_dir=sol_dir)
        else:
            s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-merged-c{c}.h5 --h5_tables {sol_dir}/cal-tec-RSm-c{c}.h5 {sol_dir}/cal-tec-CSm-c{c}.h5 \
                   --h5_time_freq {sol_dir}/cal-tec-RS-c{c}.h5 --no_antenna_crash {usepol} --propagate_flags', log='h5_merger.log', commandType='python')
            s.run(check=True)
            lib_util.run_losoto(s, f'tec-merged-c{c}', f'{sol_dir}/cal-tec-merged-c{c}.h5',
                                [f'{parset_dir}/losoto-plot-scalarph.parset'], plots_dir=f'{plot_dir}/plots-tec-merged-c{c}',
                                h5_dir=sol_dir)

    facetregname = f'{sol_dir}/facets-c{c}.reg'

    channels_out = MSs.getChout(4.e6) if MSs.getChout(4.e6) > 1 else 2
    with w.if_todo('c%02i-imaging' % c):
        logger.info('Preparing region file...')
        s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+f' --h5 {sol_dir}/cal-tec-merged-c{c}.h5 --imsize {imgsizepix_wide} \
            --pixelscale {pixscale} --writevoronoipoints --output {facetregname}', log='facet_generator.log', commandType='python')
        s.run()

        imagename = 'img/wide-' + str(c)
        imagenameM = 'img/wideM-' + str(c)
        # common imaging arguments used by all of the following wsclean calls
        widefield_kwargs = dict(data_column='CORRECTED_DATA', size=imgsizepix_wide, scale=f'{pixscale}arcsec', weight='briggs -0.5', niter=1000000,
                                gridder='wgridder',  parallel_gridding=32, minuv_l=30, mgain=0.85, parallel_deconvolution=1024,
                                join_channels='', fit_spectral_pol=3, channels_out=channels_out, deconvolution_channels=3, multiscale='',
                                multiscale_scale_bias=0.65, pol='i', facet_regions=facetregname)
        # for low-freq data, allow the beam to be fitted, otherwise (survey) force 15"
        if not(np.mean(MSs.getFreqs()) < 50e6):
            widefield_kwargs['beam_size'] = 15

        if c < 2: # cylce 0 and 1 only dd-phase
            widefield_kwargs['apply_facet_solutions'] = f'{sol_dir}/cal-tec-merged-c{c}.h5 phase000'
        else:
            widefield_kwargs['apply_facet_solutions'] = f'{sol_dir}/cal-tec-merged-c{c}.h5 phase000,amplitude000'

        # TODO make this faster - experiment with increased parallel-gridding as well as shared facet reads option
        # c0: make quick initial image to get a mask
        if c==0:
            logger.info('Making wide-field image for clean mask...')
            lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename,  no_update_model_required='',
                                 auto_threshold=5.0, auto_mask=8.0, multiscale_max_scales=3, nmiter=6, **widefield_kwargs)
            # make initial mask
            current_best_mask = make_current_best_mask(imagename, mask_threshold[c], userReg)
            # safe a bit of time by reusing psf and dirty in first iteration
            reuse_kwargs = {'reuse_psf':imagename, 'reuse_dirty':imagename}
        else:
            current_best_mask = f'img/wideM-{c-1}-mask.fits' # is this already set by the make_current_best_mask() below?
            reuse_kwargs = {}

        # main wsclean call, with mask now
        logger.info('Making wide field image ...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM,  fits_mask=current_best_mask,
                             save_source_list='', update_model_required='',  nmiter=12,  auto_threshold=2.0, auto_mask=4.0,
                             apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='',
                             local_rms='', local_rms_window=50, local_rms_strength=0.5, **widefield_kwargs, **reuse_kwargs)

        current_best_mask = make_current_best_mask(imagenameM, mask_threshold[c]-0.5, userReg)

        # reset NaNs if present
        im = lib_img.Image(f'{imagenameM}-MFS-image.fits')
        im.nantozeroModel()

    #####################################################################################################
    # Find calibration solutions for subfield (does not touch CORRECTED_DATA or CORRECTED_DATA_FR)
    if c < 2:
        # User provided subfield
        if subfield:
            if len(Regions.read(subfield)) > 1:
                raise ValueError(f'Manual subfield region {subfield} contains more than one region.')
            os.system(f'cp {subfield} {subfield_path}')

        # Generate a subdiels
        if not os.path.exists(subfield_path):
            sm = lsmtool.load(f'img/wideM-{c}-sources.txt')
            # sm.remove('img/wide-lr-mask.fits=1')  # remove sidelobe sources that were subtracted
            sm.remove('MajorAxis > 80')  # remove largest scales
            field_center1, field_size1 = lib_dd.make_subfield_region(subfield_path, MSs.getListObj()[0], sm,
                                                                         subfield_min_flux, debug_dir='img/')
            
        subfield_reg = Regions.read(subfield_path)[0]
        subfield_center = [subfield_reg.center.ra, subfield_reg.center.dec]
        subfield_size = np.max([subfield_reg.width.to_value('deg'), subfield_reg.height.to_value('deg')])
        
        with w.if_todo('c%02i_extreg_prepare' % c):
            # prepare model of central/external regions combining the subfield_reg and the mask of
            # sources so not to cut anything at the edges
            logger.info('Blanking central region of model files and reverse...')
            current_best_mask = f'img/wideM-{c}-mask.fits'
            subfield_intmask='ddparallel/images/subfieldint-mask.fits'
            os.system(f'cp {current_best_mask} {subfield_intmask}')
            lib_img.blank_image_reg(subfield_intmask, subfield_path, blankval = 1.)
            lib_img.select_connected_island(subfield_intmask, subfield_center)
            subfield_extmask='ddparallel/images/subfieldext-mask.fits'
            os.system(f'cp {current_best_mask} {subfield_extmask}')
            lib_img.blank_image_fits(subfield_extmask, subfield_intmask, blankval = 0.)

            for im in glob.glob(f'img/wideM-{c}*model*.fits'):
                wideMint = im.replace('wideM','wideMint')
                os.system('cp %s %s' % (im, wideMint))
                lib_img.blank_image_fits(wideMint, subfield_intmask, blankval = 0., inverse=True)
                wideMext = im.replace('wideM','wideMext')
                os.system('cp %s %s' % (im, wideMext))
                lib_img.blank_image_fits(wideMext, subfield_extmask, blankval = 0., inverse=True)
        # DONE

        with w.if_todo('c%02i_xtreg_subtract' % c):
            # Predict (with dd-sol) external region for subtraction - MSs: MODEL_DATA
            logger.info('Predict corrupted model of external region (wsclean)...')
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMext-{c} -j {s.max_cpucores} -channels-out {channels_out} \
                    -facet-regions {facetregname}  -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
                    -apply-facet-solutions {sol_dir}/cal-tec-merged-c{c}.h5 phase000 {MSs.getStrWsclean()}',
                log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
            s.run(check=True)

            # cycle > 0: need to add DI-corruption on top (previous iteration sub-field) - MSs: MODEL_DATA -> MODEL_DATA
            if c > 0:
                logger.info('Add previous iteration sub-field corruption on top of DD-corruption...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb={sol_dir}/cal-tec-sf-merged-c' + str(c-1) + '.h5 cor.correction=phase000 cor.invert=False',
                        log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
                # Corrupt MSs:MODEL_DATA -> MODEL_DATA (sf corr, dieamp corr)
                corr_die_amp(MSs, col='SUBFIELD_DATA', fulljones=fulljones, invert=False)

            # subtract external region MSs: CORRECTED_DATA_FR - MODEL_DATA -> SUBFIELD_DATA
            MSs.addcol('SUBFIELD_DATA','CORRECTED_DATA_FR')
            logger.info('Subtracting external region model (SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        ### DONE

        with w.if_todo('c%02i_intreg_predict' % c):
            # Predict internal region - MSs: MODEL_DATA            
            logger.info('Predict model of internal region...')
            s.add(f'wsclean -predict -padding 1.8 -name img/wideMint-{c} -j {s.max_cpucores} -channels-out {channels_out} \
                    {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
            s.run(check=True)
        ### DONE

        if c>0:
            with w.if_todo('c%02i_subfierld_corr_diamp' % c):
                # Correct MSs:SUBFIELD_DATA -> SUBFIELD_DATA (sf corr, dieamp corr)
                corr_die_amp(MSs, col='SUBFIELD_DATA', fulljones=fulljones)
            ### DONE

        with w.if_todo('c%02i_subfield_solve_tecRS' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info(f'Solving {sf_phaseSolMode} (fast RS)...')
            solve_iono(MSs, c, '2-sf', ['MODEL_DATA'], 1.5, base_solint, sf_phaseSolMode, resetant='intermediate')
        ### DONE

        with w.if_todo('c%02i_subfield_corr_tecRS' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if sf_phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec2-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if sf_phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec2-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        with w.if_todo('c%02i_subfield_solve_tecCS' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info(f'Solving {sf_phaseSolMode} (mid RS)...')
            solve_iono(MSs, c, '1-sf', ['MODEL_DATA'], 5.0, 4*base_solint, sf_phaseSolMode, resetant='inner')
        ### DONE

        with w.if_todo('c%02i_subfield_corr_tecCS' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if sf_phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec1-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if sf_phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec1-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        with w.if_todo('c%02i_subfield_solve_tecCS0' % c):
            # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
            MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
            # solve ionosphere phase - ms:SMOOTHED_DATA
            logger.info(f'Solving {sf_phaseSolMode} (slowCS)...')
            solve_iono(MSs, c, '0-sf', ['MODEL_DATA'], 10.0, 16*base_solint, sf_phaseSolMode)
        ### DONE

        with w.if_todo('c%02i_subfield_corr_tecCS0' % c):
            # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct subfield iono...')
            if sf_phaseSolMode in ['tec', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec0-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            if sf_phaseSolMode in ['phase', 'tecandphase']:
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec0-sf-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE

        # merge solutions into one h5parms for plotting and apply
        with w.if_todo('c%02i_merge_h5_subfield' % c):
            lib_util.check_rm(f'{sol_dir}/cal-tec-sf-merged-c{c}.h5')
            # reference, unflag and reset the added stations
            s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-sf-merged-c{c}.h5 --h5_tables {sol_dir}/cal-tec0-sf-c{c}.h5 {sol_dir}/cal-tec1-sf-c{c}.h5 {sol_dir}/cal-tec2-sf-c{c}.h5 \
                   --h5_time_freq {sol_dir}/cal-tec2-sf-c{c}.h5 --no_antenna_crash --no_pol', log='h5_merger.log',
                commandType='python')
            s.run(check=True)
            lib_util.run_losoto(s, f'tec-sf-merged-c{c}', f'{sol_dir}/cal-tec-sf-merged-c{c}.h5',
                                [f'{parset_dir}/losoto-plot-scalarph.parset'],
                                plots_dir=f'{plot_dir}/plots-tec-sf-merged-c{c}', h5_dir=sol_dir)
        ### DONE

        # Do a quick image of the subfield, not strictly necessary but good to have...
        with w.if_todo('c%02i_image-subfield' % c):
            logger.info('Test image subfield...')
            lib_util.run_wsclean(s, 'wscleanSF-c'+str(c)+'.log', MSs.getStrWsclean(), name=f'img/subfield-{c}', data_column='SUBFIELD_DATA', size=int(1.2*subfield_size*3600/pixscale), scale=f'{pixscale}arcsec',
                                 weight='briggs -0.5', niter=100000, gridder='wgridder',  parallel_gridding=6, shift=f'{subfield_center[0].to(u.hourangle).to_string()} {subfield_center[1].to_string()}',
                                 no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                                 join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6), deconvolution_channels=3, baseline_averaging='',
                                 multiscale='',  multiscale_scale_bias=0.7, pol='i')
        ### DONE

        #####################################################################################################
        # Subtract side-lobe sources (modifies CORRECTED_DATA_FR removing sidelobe sources)
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
                # Recreate MODEL_DATA of internal region for subtraction
                logger.info('Predict corrupted model of internal region...')
                s.add(f'wsclean -predict -padding 1.8 -name img/wideMintpb-0 -j {s.max_cpucores} -channels-out {channels_out} \
                       -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam -facet-regions {facetregname} \
                       -apply-facet-solutions {sol_dir}/cal-tec-merged-c{c}.h5 phase000 {MSs.getStrWsclean()}',
                    log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
                s.run(check=True)

                # subtract internal region from MSs: CORRECTED_DATA_FR - MODEL_DATA -> SUBFIELD_DATA
                logger.info('Subtract main-lobe (SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set SUBFIELD_DATA = CORRECTED_DATA_FR - MODEL_DATA"',
                        log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
                if develop: clean_empty(MSs, 'only_sidelobe', 'SUBFIELD_DATA', size=10000)
            ### DONE

            # Do a rough correction of the sidelobe data using the subfield solutions
            # MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            with w.if_todo(f'correct-sidelobe'): # just for testing/debug
                logger.info('Correct sidelobe data with subfield iono solutions...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                        cor.parmdb={sol_dir}/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000',
                        log='$nameMS_sf-correct.log', commandType='DP3')
            ### DONE

            imagename_lr = 'img/wide-lr'
            # Image the sidelobe data
            # MSs: create MODEL_DATA (with just the sidelobe flux)
            with w.if_todo('image_lr'):
                logger.info('Cleaning sidelobe low-res...')
                # TODO test if IDG beam correction is fast enough, otherwise ignore beam corr or use facet beam [for idg remove parallel gridding (parallel_gridding=4) and BLavg (baseline_averaging='')]
                channels_out_lr = MSs.getChout(2.e6) if MSs.getChout(2.e6) > 1 else 2
                lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, do_predict=True, data_column='SUBFIELD_DATA',
                                     size=imgsizepix_lr, scale='30arcsec', save_source_list='',  parallel_gridding=4, baseline_averaging='',
                                     #use_idg='', idg_mode='cpu', grid_with_beam='', beam_aterm_update=900, use_differential_lofar_beam='',
                                     weight='briggs -0.5', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
                                     taper_gaussian='200arcsec', mgain=0.85, channels_out=channels_out_lr, parallel_deconvolution=512,
                                     local_rms='', auto_mask=3, auto_threshold=1.5, join_channels='', fit_spectral_pol=5)
            ### DONE

            # Subtract full low-resolution field (including possible large-scale emission within primary beam)
            # to get empty data set for flagging
            # MSs: SUBFIELD_DATA - MODEL_DATA -> SUBFIELD_DATA
            with w.if_todo('subtract_lr'):
                logger.info('Subtract low-resolution to get empty data set (SUBFIELD_DATA = SUBFIELD_DATA - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set SUBFIELD_DATA = SUBFIELD_DATA - MODEL_DATA"',
                    log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
                if develop: clean_empty(MSs, 'empty', 'SUBFIELD_DATA')
            ### DONE

            # Flag on residuals (SUBFIELD_DATA)
            with w.if_todo('flag_residuals_lr'):
                logger.info('Flagging residuals (SUBFIELD_DATA)...')
                MSs.run('DP3 ' + parset_dir + '/DP3-flag.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA \
                    aoflagger.strategy=' + parset_dir + '/LBAdefaultwideband.lua',
                    log='$nameMS_flag-c' + str(c) + '.log', commandType='DP3')
            ### DONE

            # Now subtract the sidelobe
            with w.if_todo('subtract_sidelobe'):
                sidelobe_predict_mode = 'wsclean'
                if sidelobe_predict_mode == 'wsclean':
                    # blank within main-libe to not subtract anything from there (maybe extended sources not properly subtracted at the beginning)
                    for im in glob.glob(f'{imagename_lr}*model*fits'):
                        wideLRext = im.replace(imagename_lr, f'{imagename_lr}-blank')
                        os.system('cp %s %s' % (im, wideLRext))
                        lib_img.blank_image_reg(wideLRext, beamReg , blankval=0.)

                    logger.info('Predict model of sidelobe region (wsclean)...')
                    # -use-idg -idg-mode cpu -grid-with-beam -beam-aterm-update 900 -use-differential-lofar-beam
                    s.add(f'wsclean -predict -padding 1.8 -name {imagename_lr}-blank -j {s.max_cpucores} \
                        -channels-out {channels_out_lr} {MSs.getStrWsclean()}',
                        log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
                    s.run(check=True)
                elif sidelobe_predict_mode=='DP3':
                    logger.info('Predict model of sidelobe region (DP3)...')
                    sidelobe_sky = lsmtool.load(f'{imagename_lr}-sources.txt')
                    sidelobe_sky.remove(f'{beamMask}==True')
                    MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={imagename_lr}-sources.txt msout.datacolumn=MODEL_DATA',
                        log='$nameMS_pre.log', commandType='DP3')
                
                logger.info('Corrupt sidelobe model with subfield solutions...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.parmdb={sol_dir}/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000 cor.invert=False',
                        log='$nameMS_sidelobe_corrupt.log', commandType='DP3')

                logger.info('Subtract corrupted sidelobe model (CORRECTED_DATA_FR = CORRECTED_DATA_FR - MODEL_DATA)...')
                MSs.run('taql "update $pathMS set CORRECTED_DATA_FR = CORRECTED_DATA_FR - MODEL_DATA"', log='$nameMS_taql-c' + str(c) + '.log', commandType='general')
            ### DONE

        # apply the subfield solutions to the data.
        with w.if_todo('c%02i_corr_sf_sols' % c):
            logger.info('Correct subfield ionosphere (CORRECTED_DATA_FR -> CORRECTED_DATA)...')
            # Correct MSs:CORRECTED_DATA_FR -> CORRECTED_DATA (sf corr)
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA_FR  \
                    cor.parmdb={sol_dir}/cal-tec-sf-merged-c' + str(c) + '.h5 cor.correction=phase000',
                    log='$nameMS_sf-correct.log', commandType='DP3')
        ### DONE
        
        # finally re-correct for die amp on the newly created CORRECTED_DATA
        if c >= 1:
            with w.if_todo('c%02i_corr_dieamp' % c):
                # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA (sf corr, dieamp corr)
                corr_die_amp(MSs, col='CORRECTED_DATA', fulljones=fulljones)
            ### DONE

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image*.fits ddparallel/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-MFS-residual*.fits ddparallel/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources*.txt ddparallel/images') for c in range(maxIter) ]
# debugging images -> can be removed in production
[ os.system('mv img/subfield-'+str(c)+'-MFS-image*.fits ddparallel/images') for c in range(maxIter) ]
# os.system('mv img/wideP-MFS-*-image.fits ddparallel/images')
# os.system('mv img/wide-lr-MFS-image.fits ddparallel/images')

# debug images
if develop:
    os.system('mv img/only*image.fits ddparallel/images')
    os.system('mv img/empty*image.fits ddparallel/images')

# Copy model
os.system(f'mv img/wideM-{maxIter-1}-*-model.fits ddparallel/skymodel')
os.system(f'mv img/wideM-{maxIter-1}-*-model-fpb.fits ddparallel/skymodel')
os.system(f'mv img/wideM-{maxIter-1}-*-model-pb.fits ddparallel/skymodel')

w.alldone()
