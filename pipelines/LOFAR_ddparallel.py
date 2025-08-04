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

# Waiting for bug fixes in other software
# TODO add BDA
# TODO add timesmearing in dp3 predict
# TODO model_weighted_constraint

import os, sys, glob, random, json
import numpy as np
from regions import Regions
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs
from losoto.h5parm import h5parm
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, lib_h5, lib_dd_parallel, lib_cat
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
fulljones = parset.getboolean('LOFAR_ddparallel', 'fulljones') # do fulljones DIE amp correction instead of diagonal (default: False)
min_facets = parset.get('LOFAR_ddparallel', 'min_facets') # ''=default (differs for SPARSE and OUTER), otherwise provide comma seperated list [2,3,6..]
max_facets = parset.get('LOFAR_ddparallel', 'max_facets') # ''=default (differs for SPARSE and OUTER), otherwise provide comma seperated list [5,10,20..]
develop = parset.getboolean('LOFAR_ddparallel', 'develop') # for development, make more output/images
use_shm = parset.getboolean('LOFAR_ddparallel', 'use_shm') # use shared memory for wsclean
data_dir = parset.get('LOFAR_ddparallel','data_dir')
start_sourcedb = parset.get('model','sourcedb')
userReg = parset.get('model','userReg')
sf_phaseSolMode = 'phase' #'tec'

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

def corrupt_model_dirs(MSs, c, model_columns, solmode='phase'):
    """ CORRUPT the MODEL_DATA columns for model_columns
    Parameters
    ----------
    MSs: MSs object
    c: int, cycle
    model_columns: list of patch names
    solmode: which solution mode to corrupt for (phase, tec, tecandphase)
    """
    logger.info(f'Corrupt models: {solmode}...')
    for model_column in model_columns:
        if solmode in ['tec', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}]  \
                    cor.parmdb={sol_dir}/cal-tec-c{c}.h5 cor.correction=tec000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        elif solmode in ['phase', 'tecandphase']:
            MSs.run(
                f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                    cor.parmdb={sol_dir}/cal-tec-c{c}.h5 cor.correction=phase000 cor.invert=False',
                log='$nameMS_corrupt.log', commandType='DP3')
        elif solmode in ['amplitude']:
            if not model_column.startswith('patch'):    
                MSs.run(
                    f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={model_column} msout.datacolumn={model_column} cor.direction=[{model_column}] \
                        cor.parmdb={sol_dir}/cal-amp-3C.h5 cor.correction=amplitude000 cor.invert=False',
                    log='$nameMS_corrupt.log', commandType='DP3')

def corr_die_amp(MSs, col='CORRECTED_DATA', fulljones=fulljones, invert=True):
    """
    invert: set to False to corrupt the chosen column
    """
    if fulljones:
        logger.info('Correct amp-di (fulljones)...') if invert else logger.info('Corrupt amp-di (fulljones)...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={col} msout.datacolumn={col} \
                cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=fulljones cor.soltab=[amplitudeSmooth,phaseSmooth] \
                cor.updateweights=False cor.invert={invert}',
                log='$nameMS_diampcor.log', commandType='DP3')
                
        logger.info('Perform amp-di normalisation...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn={col} msout.datacolumn={col} \
                cor.parmdb={sol_dir}/cal-amp-dinorm.h5 cor.correction=amplitude000 \
                cor.updateweights=False cor.invert={invert}',
                log='$nameMS_diampcor.log', commandType='DP3')
    else:
        logger.info('Correct amp-di...') if invert else logger.info('Corrupt amp-di...')
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

def add_3c_models(sm, phasecentre, null_mid_freq, beamMask, max_sep=50., threshold=0):
    
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
            else:
                within_beam = False
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
            sourcedb = os.path.dirname(__file__) + '/../models/calib-simple.skymodel'
            sm_3c = lsmtool.load(sourcedb, beamMS=sm.beamMS)
            sm_3c.select(f'patch=={source.replace(" ","")}')
            sm_3c.setColValues('LogarithmicSI', ['True']*len(sm_3c.getColValues('I'))) # add standard spidx
        else:
            sourcedb = os.path.dirname(__file__) + f'/../models/3CRR/{source.replace(" ","")}.txt'
            if not os.path.exists(sourcedb):
                logger.warning(f'No model found for {source} (seperation {sep:.2f} deg)')
                continue
            sm_3c = lsmtool.load(sourcedb, beamMS=sm.beamMS)
        
        sm_3c.setColValues("Patch", ["source_"+source.replace(" ","")]*len(sm_3c.getColValues("I")))
        flux_3c =  sm_3c.getColValues("I", aggregate="sum", applyBeam=True)[0]
        if flux_3c > threshold:
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
    lib_util.check_rm(f'ddparallel/skymodel/patches_c{c}.reg')
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
        os.system(f'cat ddparallel/skymodel/regions_c{c}/*.reg >> ddparallel/skymodel/patches_c{c}.reg')
    lib_util.check_rm(f'ddparallel/skymodel/regions_c{c}')

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

MSs = lib_ms.AllMSs( glob.glob(data_dir + 'mss/TC*[0-9].MS'), s, check_flags=True, check_consistency=True)
MSs.print_HAcov()
[MS.print_ateam_demix() for MS in MSs.getListObj()]
for ateam in ['CasA', 'CygA', 'TauA', 'VirA']:
    dist = MSs.getListObj()[0].distBrightSource(ateam)
    logger.info('Distance from %s: %.0f deg' % (ateam, dist))

# print fractional flags
for MS in MSs.getListObj():
    logger.info(f'{MS.nameMS}: Fractional flags: {MS.fractionalFlag()*100:.1f}%.')

# make beam to the first mid null - outside of that do a rough subtraction and/or 3C peeling. Use sources inside for calibration
phasecentre = MSs.getListObj()[0].getPhaseCentre()
null_mid_freq = max(MSs.getListObj()[0].getFWHM(freq='mid', elliptical=True)) * 1.8 # FWHM to null

# set image size - this should be a bit more than the beam region used for calibration
pixscale = MSs.getListObj()[0].getPixelScale()
imgsizepix_wide = int(1.85*max(MSs.getListObj()[0].getFWHM(freq='min', elliptical=True))*3600/pixscale) # roughly to biggest null
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
nchan_ph = round(0.195312e6 / MSs.getListObj()[0].getChanband())  # number of channels in 1 SBs
nchan_amp = 2* round(0.192e6 / MSs.getListObj()[0].getChanband()) # number of channels in 2 SBs
tint = MSs.mssListObj[0].getTimeInt()
if np.round(tint) != 4:
    raise ValueError('Input data should be at 4s time resolution.')

mask_threshold = [5.0,4.5,4.0,4.0,4.0,4.0] # sigma values for beizorro mask in cycle c

# define list of facet fluxes per iteration -> this can go into the config
# if we have LOTSS-DR3 or a custom sky model, we will start from 3Jy not 4Jy sources!
if 'OUTER' in MSs.getListObj()[0].getAntennaSet():
    facet_fluxes = np.array([4, 2.2, 1.2, 1.0, 0.9, 0.8])*(54e6/np.mean(MSs.getFreqs()))**0.7 # this is not the total flux, but the flux of bright sources used to construct the facets. still needs to be tuned, maybe also depends on the field
elif 'SPARSE' in MSs.getListObj()[0].getAntennaSet():
    facet_fluxes = np.array([4, 2.6, 1.3, 1.1, 1.0, 0.9])*(54e6/np.mean(MSs.getFreqs()))**0.7 # this is not the total flux, but the flux of bright sources used to construct the facets. still needs to be tuned, maybe also depends on the field

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

if max_facets: # if manually provided
    if not isinstance(max_facets, list):
        max_facets = max_facets.replace('[', '').replace(']', '').split(',')
        max_facets = np.array(max_facets).astype(int)
else: #default settings
    # use more facets for SPARSE (larger FoV)
    if 'SPARSE' in MSs.getListObj()[0].getAntennaSet():
        max_facets = [10, 22, 35, 35, 35, 35]
    elif 'OUTER' in MSs.getListObj()[0].getAntennaSet():
        max_facets = [8, 16, 25, 25, 25, 25]
    else:
        raise ValueError(f'{MSs.getListObj()[0].getAntennaSet()} not recognized.')

if (min_facets[0] > max_facets[0]) or (min_facets[1] > max_facets[1]):
    raise ValueError(f'min_facets {min_facets} and max_facets {max_facets} are not compatible.')


# Make beam mask/reg
beamMask = 'ddparallel/beam.fits'
beamReg = 'ddparallel/beam.reg'
MSs.getListObj()[0].makeBeamReg(beamReg, freq='min', to_pbval=0)
if not os.path.exists(beamMask):
    logger.info('Making mask of primary beam...')
    lib_util.run_wsclean(s, 'wscleanLRmask.log', MSs.getStrWsclean(), name=beamMask.replace('.fits',''), data_column='DATA', \
                         size=imgsizepix_lr, scale='30arcsec')
    os.system(f'mv {beamMask.replace(".fits","-image.fits")} {beamMask}') # beam-image.fits -> beam.fits
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 1.)
    lib_img.blank_image_reg(beamMask, beamReg, blankval = 0., inverse=True)
subfield_path = 'ddparallel/skymodel/subfield.reg'

with w.if_todo('addcol'):
    MSs.addcol('CORRECTED_DATA', 'DATA') # carries on varies correction
    MSs.addcol('PREPARED_DATA', 'DATA') # used to store data with annoying sources (3c, sidelobe) subtracted

#################################################################################################
# Find FR, all it does is creating the cal-fr.h5
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

    MSs.run(f'DP3 {parset_dir}/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5',
            log='$nameMS_solFR.log', commandType="DP3")
    lib_util.run_losoto(s, 'fr', [ms + '/fr.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-fr.parset'])
    os.system(f'mv cal-fr.h5 {sol_dir}')
    os.system(f'mv plots-fr {plot_dir}')

    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, MODEL_DATA_FR"',
            log='$nameMS_taql.log', commandType='general')

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
            if start_sourcedb == '':  # if not provided, use LOTSS-DR3 as default, if the field is not fully covered, resort to GSM
                if lib_dd_parallel.check_lotss_coverage(phasecentre, null_mid_freq/2):
                    logger.info('Target fully in LoTSS-DR3 - start from LoTSS + use more initial facets.')
                    start_sourcedb = 'LOTSS-DR3'
                    facet_fluxes[0] = 3.5
                else:
                    logger.info('Target not fully in LoTSS-DR3 - start from GSM.')
                    start_sourcedb = 'GSM'
            # case use survey to start
            if start_sourcedb.upper() in ['GSM','LOTSS','TGSS','VLSSR','NVSS','WENSS']:
                logger.info(f'Get skymodel from {start_sourcedb}...')
                sm = lsmtool.load(start_sourcedb, VOPosition=phasecentre, VORadius=null_mid_freq/2, beamMS=beamMS)
                if start_sourcedb.upper() == 'LOTSS':
                    sm.setColValues('I', sm.getColValues('I')/1000) # convert mJy to Jy TODO fix in LSMtool
                    sm.select('I>0.05', applyBeam=True)
                    sm.setColValues('SpectralIndex', [[-0.7]]*len(sm.getColValues('I'))) # add standard spidx
                    sm.setColValues('LogarithmicSI', ['True']*len(sm.getColValues('I'))) # add standard spidx
            # load LoTSS DR3 model and select decrease size in DEC
            # using components downloaded from https://www.lofar-surveys.org/downloads/DR3/catalogues/LoTSS_DR3_v0.1_gaus.fits
            # componentlist prepared on 11-03-2025 (total flux in Jy)
            elif start_sourcedb.upper() == 'LOTSS-DR3':
                sm = lib_cat.get_LOTSS_DR3_cone_as_skymodel(phasecentre, null_mid_freq/2, 'ddparallel/skymodel/starting.skymodel', beamMS=beamMS)
            # otherwise if provided, use manual model
            else:
                logger.info(f'Using input skymodel {start_sourcedb}')
                sm = lsmtool.load(start_sourcedb, beamMS=beamMS)
        else:
                # get wsclean skymodel of previous iteration
                wsc_src = f'img/wideM-{c-1}-sources-pb.txt'
                sm = lsmtool.load(wsc_src, beamMS=beamMS)
        # if using e.g. LoTSS, adjust for the frequency
        logger.debug(f'Extrapolating input skymodel fluxes from {np.mean(sm.getColValues("ReferenceFrequency"))/1e6:.0f}MHz to {np.mean(MSs.getFreqs())/1e6:.0f}MHz assuming si=-0.7')
        si_factor = (np.mean(MSs.getFreqs())/np.mean(sm.getColValues('ReferenceFrequency')))**0.7 # S60 = si_factor * S54
        sm.select(f'{beamMask}==True')  # remove outside of FoV (should be subtracted (c>0) or not present (c==0)!)
        sm.group('threshold', FWHM=2/60, root='Src') # group nearby components to single source patch
        sm.setPatchPositions(method='wmean', applyBeam=True)
        sm = lib_dd_parallel.merge_nearby_bright_facets(sm, 1/60, 0.5, applyBeam=True)
        # TODO we need some logic here to avoid picking up very extended sources.
        patch_fluxes = sm.getColValues('I', aggregate='sum', applyBeam=True)
        # disbale 3cremoval if bright source in the field and set max_facets to 1
        if c==0 and (patch_fluxes/si_factor > 100).any():
            logger.warning(f'Found patch with flux > 100 Jy: {patch_fluxes[patch_fluxes > 100]}, turning off 3c removal and forcing max_facets=1.')
            remove3c = False
            min_facets[0] = 1
            max_facets[0] = 1
            mask_threshold = [7.0,6.5,6.0,6.0,6.0,6.0] # sigma values for beizorro mask in cycle c           
        # check if there are less than the minimum requested bright sources to form the facets
        if sum(patch_fluxes/si_factor > facet_fluxes[c]) < min_facets[c]: # convert skymodel fluxes to MS central freq
            bright_sources_flux = np.sort(patch_fluxes)[-min_facets[c]] / si_factor # bright sources flux is at MSs central freq
            logger.warning(f'Less than {min_facets[c]} bright sources above minimum flux {facet_fluxes[c]:.2f} Jy! Using sources above {bright_sources_flux:.2f} Jy')
        # check if there are more than the maximum requested bright sources to form the facets
        elif sum(patch_fluxes/si_factor > facet_fluxes[c]) > max_facets[c]:
            bright_sources_flux = np.sort(patch_fluxes)[-max_facets[c]] / si_factor
            logger.warning(f'More than {max_facets[c]} bright sources above minimum flux {facet_fluxes[c]:.2f} Jy! Using sources above {bright_sources_flux:.2f} Jy')
        else:
            bright_sources_flux = facet_fluxes[c]
        bright_names = sm.getPatchNames()[patch_fluxes >= bright_sources_flux*si_factor]
        bright_pos = sm.getPatchPositions(bright_names)
        sm.group('voronoi', targetFlux=bright_sources_flux*si_factor, applyBeam=True, root='', byPatch=True)
        sm.setPatchPositions(bright_pos)
        lib_dd_parallel.rename_skymodel_patches(sm, applyBeam=True)

        if c == 0 and remove3c:
            # Add models of bright 3c sources to the sky model. Model will be subtracted from data before imaging.
            sm = add_3c_models(sm, phasecentre=phasecentre, beamMask=beamMask, null_mid_freq=null_mid_freq)
            for stokes in ['Q','U','V']:
                sm.setColValues(stokes, np.zeros(len(sm.getColValues("I")))) # force non I Stokes to zero
                sm.table[stokes].unit = u.Jy

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
            # Add model to MODEL_DATA and do FR corruption
            # TODO add time smearing in the predict parset?
            # note: beammode=full applies array_beam and element_beam, but the element_beam at the centre is already corrected, so MODEL_DATA needs to be
            # corrected for element_beam at the phase centre. However, this is very slow, so we just corrupt for the array_beam and assume the element beam
            # correction at the phase centre is ok everywhere.
            logger.info(f'Add model to {patch}...')
            correctfreqsmearing = c == 0 # only in cycle zero correct freq smearing
            MSs.run(f'DP3 {parset_dir}/DP3-predict-beam.parset msin=$pathMS pre.sourcedb=$pathMS/{sourcedb_basename} pre.sources={patch} msout.datacolumn={patch} pre.correctfreqsmearing={correctfreqsmearing}',
                    log='$nameMS_pre.log', commandType='DP3')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn={patch} msout.datacolumn={patch}\
                       cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")
            # pos = sm.getPatchPositions()[patch]
            # size = int((1.1*sm.getPatchSizes()[np.argwhere(sm.getPatchNames()==patch)]) // 4)
            # logger.info(f'Test image MODEL_DATA...')
            # lib_util.run_wsclean(s, 'wscleanMODEL-c' + str(c) + '.log', MSs.getStrWsclean(), name=f'ddparallel/skymodel/{patch}',
            #                      data_column=patch, size=size, scale='8arcsec', shift=f'{pos[0].to(u.hourangle).to_string()} {pos[1].to_string()}',
            #                      weight='briggs -0.5', niter=10000, gridder='wgridder', parallel_gridding=MSs.getChout(4.e6), no_update_model_required='', minuv_l=30, mgain=0.9,
            #                      parallel_deconvolution=512, beam_size=15, join_channels='', fit_spectral_pol=3,
            #                      channels_out=MSs.getChout(4.e6), deconvolution_channels=3, pol='i', nmiter=5 )
        #MSs = lib_ms.AllMSs(glob.glob("*-smearing.MS"), s, check_flags=True)
    ### DONE

    ### solve ionosphere phase - ms:SMOOTHED_DATA - > reset for all BUT most distant RS!
    with w.if_todo('c%02i_solve_iono' % c):
        # Smooth MSs:CORRECTED_DATA -> SMOOTHED_DATA
        MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
        logger.info('Solving ionosphere (DD)...')
        smMHz = np.array([[2.5,4.0,10.0,15.0],[6.0,10.0,15.0,25.0]]) # [cycle0, cycle1]
        smMHz_factors = [smMHz[0]/np.max(smMHz[0]), smMHz[1]/np.max(smMHz[1])] # factors should be <1 otherwise trimming of kernel is off
        solutions_per_direction = 15*np.ones(len(patches), dtype=int)
        # get twice as many solutions for brighter directions solint (i.e. one per time step) for bright directions
        solutions_per_direction[patch_fluxes > 4] *= 2
        # solint = int(solint / np.min(solutions_per_direction))
        # solutions_per_direction = (solutions_per_direction / np.min(solutions_per_direction)).astype(int)

        # if len(patches) >= 30 and solint > 60:
        #     logger.warning('Detected many directions - limit number of parallel DP3 processes to 1.')
        #     maxProcs = 1
        # elif len(patches) >= 18 and solint > 60:
        #     logger.warning('Detected many directions - limit number of parallel DP3 processes to 3.')
        #     maxProcs = 3
        # else:
        #     maxProcs = 8
        # TODO use smoothness_dd_factors
        nchan_ph = round(0.195312e6 / MSs.getListObj()[0].getChanband())  # number of channels in 1 SBs
        avg_factors = [15,5,2,1]
        ant_avg_factors = f"[CS*:{avg_factors[0]},[RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LBA]:{avg_factors[1]},[RS208LBA,RS307LBA,RS406LBA,RS407LBA]:{avg_factors[2]},[RS210LBA,RS310LBA,RS409LBA,RS508LBA,RS509LBA]:{avg_factors[3]}]"
        ant_smooth_factors = f"[CS*:{smMHz_factors[c][3]},[RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LBA]:{smMHz_factors[c][2]},[RS208LBA,RS307LBA,RS406LBA,RS407LBA]:{smMHz_factors[c][1]},[RS210LBA,RS310LBA,RS409LBA,RS508LBA,RS509LBA]:{smMHz_factors[c][0]}]"

        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec.h5 sol.solint=60 \
                  sol.mode=scalarphase sol.smoothnessconstraint={max(smMHz[c])}e6 sol.smoothnessreffrequency=54e6 sol.nchan={nchan_ph}  \
                  sol.modeldatacolumns="[{",".join(patches)}]" sol.solutions_per_direction="{np.array2string(solutions_per_direction, separator=",")}" \
                  sol.antenna_averaging_factors={ant_avg_factors} sol.antenna_smoothness_factors={ant_smooth_factors} ',
                  log='$nameMS_solTEC-c' + str(c) + '.log', commandType='DP3', maxProcs=8)

        lib_util.run_losoto(s, f'tec-c{c}', [ms + f'/tec.h5' for ms in MSs.getListStr()],
                            [parset_dir + '/losoto-refph.parset', parset_dir + '/losoto-plot-scalarph.parset'],
                            plots_dir=f'{plot_dir}/plots-tec-c{c}', h5_dir=sol_dir)

        # make sure the h5parm directions are correctly set - this should actually work automatically with DP3 -> eventually fix this in the DP3 solve call
        lib_h5.point_h5dirs_to_skymodel(f'{sol_dir}/cal-tec-c{c}.h5', sourcedb)
    ### DONE

    ### CORRUPT the MODEL_DATA columns for all patches
    with w.if_todo('c%02i_corrupt_iono' % c):
        corrupt_model_dirs(MSs, c, patches)
    ### DONE

    # merge all models into a single columns (FR corrupted)
    with w.if_todo('c%02i_add_patches' % c):
        logger.info('Setting MODEL_DATA to sum of corrupted patch models...')
        MSs.addcol('MODEL_DATA', 'DATA', usedysco=False)
        non_3c_patches = [p for p in patches if p.startswith('patch_')]
        MSs.run(f'taql "UPDATE $pathMS SET MODEL_DATA={"+".join(non_3c_patches)}"', log='$nameMS_taql.log', commandType='general')

    ########################### 3C-subtract PART ####################################
    # PREPARED_DATA and CORRECTED_DATA will have 3c models (corrupted for tec, amp and FR) removed
    if c == 0 and remove3c:
        _3c_patches = [p for p in patches if not p.startswith('patch_')]
        if len(_3c_patches) > 0:
            with w.if_todo('3c_solve_amp'):
                logger.info('Solving amplitude for 3C...')
                # Solve diagonal amplitude MSs:SMOOTHED_DATA
                # 225 is 15 minutes
                # sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS201LBA,CS301LBA,CS401LBA,CS501LBA,CS103LBA,CS302LBA]]',
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.model_weighted_constraints=true\
                          sol.mode=diagonalamplitude sol.nchan={nchan_amp} sol.smoothnessconstraint=4e6 sol.smoothnessreffrequency=54e6 sol.h5parm=$pathMS/amp-3C.h5 sol.datause=full \
                          sol.modeldatacolumns="[MODEL_DATA,{",".join(_3c_patches)}]" sol.solint='+str(224),
                          log=f'$nameMS_solamp_3c_c{c}.log', commandType="DP3", maxProcs=4)

                #losoto_parsets = [parset_dir + '/losoto-clip.parset', parset_dir + '/losoto-plot-amp.parset']
                losoto_parsets = [parset_dir + '/losoto-plot-amp.parset']
                lib_util.run_losoto(s, 'amp-3C', [ms + '/amp-3C.h5' for ms in MSs.getListStr()], losoto_parsets,
                                    plots_dir=f'{plot_dir}/plots-amp-3C', h5_dir=sol_dir)
                ### DONE

            with w.if_todo('3C_corrupt_subtract'):
                MSs.run('addcol2ms.py -m $pathMS -c FLAG_BKP -i FLAG', log='$nameMS_addcol.log', commandType='python')
                MSs.run('taql "update $pathMS set FLAG_BKP = FLAG"', log='$nameMS_taql.log', commandType='general')
                
                debug_3c_sub = False
                if debug_3c_sub:
                    MSs.addcol('DATA_SUB', 'PREPARED_DATA')
                    with open(parset_dir+"/3C_coordinates.json", "r") as file:
                        import json
                        all_3c = json.load(file)
                        
                for patch in _3c_patches:
                    
                    if debug_3c_sub:
                        MSs.run('taql "update $pathMS set DATA_SUB = PREPARED_DATA"', log='$nameMS_taql.log', commandType='general')
                        coords = all_3c[patch.replace('source_','').replace("C","C ")]
                        coords[0] = coords[0].replace(" ", "h", 1).replace(" ", "m", 1) + "s"
                        coords[1] = coords[1].replace(" ", "d", 1).replace(" ", "m", 1) + "s"
                        clean_empty(MSs,f'{patch}_data_zoom_pre', 'PREPARED_DATA', shift=coords, size=2000)
                        clean_empty(MSs,f'{patch}_model', f'{patch}', shift=coords, size=2000)
                        MSs.run(
                            f"taql 'UPDATE $pathMS SET DATA_SUB = DATA_SUB - {patch}'",
                            log = f'$nameMS_taql.log',
                            commandType = 'general'
                        )
                        clean_empty(MSs,f'{patch}_data_after_phase', 'DATA_SUB', shift=coords, size=2000)

                    logger.info(f'Subtracting {patch}...')
                    # Corrupt MODEL_DATA with amplitude, set MODEL_DATA = 0 where data are flagged, then unflag everything
                    corrupt_model_dirs(MSs, c, [patch], solmode='amplitude')
                    MSs.run(f'taql "update $pathMS set {patch}[FLAG] = 0"', log='$nameMS_taql.log', commandType='general')
                    
                    if debug_3c_sub:
                        clean_empty(MSs,f'{patch}_model_amp', f'{patch}', shift=coords, size=2000)
                        MSs.run('taql "update $pathMS set DATA_SUB = PREPARED_DATA"', log='$nameMS_taql.log', commandType='general')
                        MSs.run(
                            f"taql 'UPDATE $pathMS SET DATA_SUB = DATA_SUB - {patch}'",
                            log = f'$nameMS_subtract_{patch}.log',
                            commandType = 'general'
                        )
                        clean_empty(MSs,f'{patch}_data_after_amp', 'DATA_SUB', shift=coords, size=2000)
                        
                    MSs.run(
                        f"taql 'UPDATE $pathMS SET PREPARED_DATA = PREPARED_DATA - {patch}'",
                        log = f'$nameMS_taql.log', 
                        commandType = 'general'
                    )

                    MSs.run(
                        f"taql 'UPDATE $pathMS SET CORRECTED_DATA = CORRECTED_DATA - {patch}'",
                        log = f'$nameMS_taql.log', 
                        commandType = 'general'
                    )
                    
                    if debug_3c_sub: MSs.deletecol('DATA_SUB')
                    MSs.deletecol(patch)
                    sm = lsmtool.load(sourcedb, beamMS=beamMS)
                    sm.select(f'Patch != {patch}')
                    sm.write(sourcedb, clobber=True)
                    patches = np.delete(patches, np.where(patches == patch)[0])
                
                MSs.run('taql "update $pathMS set FLAG = FLAG_BKP"', log='$nameMS_taql.log', commandType='general')
                MSs.deletecol('FLAG_BKP')

                # remove subtracted 3c direction from h5 for wide-field imaging
                lib_util.check_rm(f'{sol_dir}/cal-tec-no3c-c{c}.h5')

                # by construction, the 3c sources are always at the last indices
                filter_directions = "--filter_directions '[" + ', '.join([str(i) for i in range(len(patches))]) + "]'"
                logger.info('Merge solutions...')
                s.add(f'h5_merger.py --h5_out {sol_dir}/cal-tec-no3c-c{c}.h5 --h5_tables {sol_dir}/cal-tec-c{c}.h5 --h5_time_freq {sol_dir}/cal-tec-c{c}.h5 \
                      --no_antenna_crash --no_pol --propagate_flags {filter_directions}', log='h5_merger.log', commandType='python')
                s.run(check=True)
                lib_util.run_losoto(s, f'tec-no3c-c{c}', f'{sol_dir}/cal-tec-no3c-c{c}.h5',
                                    [f'{parset_dir}/losoto-plot-scalarph.parset'],
                                    plots_dir=f'{plot_dir}/plots-tec-no3c-c{c}',
                                    h5_dir=sol_dir)
        else:
            remove3c = False
        ### DONE

    ########################### AMP-CAL PART ####################################
    # Only once in cycle 1: do di amp to capture element beam 2nd order effect
    # TODO add updateweights in production
    if c == 1:
        with w.if_todo('amp_di_solve'):

            if fulljones:
                logger.info('Solving amp-di (fulljones)...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=full sol.nchan={nchan_amp} sol.modeldatacolumns=[MODEL_DATA] \
                     sol.mode=fulljones sol.h5parm=$pathMS/amp-di.h5 sol.solint=150 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]',
                     log='$nameMS_diampsol.log', commandType='DP3')

                lib_util.run_losoto(s, 'amp-di', [ms + '/amp-di.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-fj.parset', f'{parset_dir}/losoto-amp-difj.parset'],
                                plots_dir=f'{plot_dir}/plots-amp-di', h5_dir=sol_dir)
            
                # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
                logger.info('Correct amp-di (CORRECTED_DATA -> CORRECTED_DATA)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=fulljones cor.soltab=[amplitudeSmooth,phaseSmooth] \
                    cor.updateweights=False',
                    log='$nameMS_diampcor.log', commandType='DP3')
                
                #sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS406LBA,RS407LBA,RS409LBA,RS310LBA,RS503LB,RS508LBA,RS509LBA]]',
                logger.info('Solving amp-di for normalisation...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.datause=full sol.nchan=0 sol.solint=0 \
                        sol.modeldatacolumns=[MODEL_DATA] sol.mode=scalaramplitude sol.h5parm=$pathMS/amp-dinorm.h5',
                        log='$nameMS_diampsol.log', commandType='DP3')

                lib_util.run_losoto(s, 'amp-dinorm', [ms + '/amp-dinorm.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-amp2d.parset'], plots_dir=f'{plot_dir}/plots-amp-di', h5_dir=sol_dir)

                with h5parm(f'{sol_dir}/cal-amp-dinorm.h5') as h5:
                    meanamp = np.nanmean(h5.getSolset('sol000').getSoltab('amplitude000').val)
                    logger.info(f'Amplitude normalization: {meanamp:.3f} ')
            
                logger.info('Correct amp-di normalisation (CORRECTED_DATA -> CORRECTED_DATA)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb={sol_dir}/cal-amp-dinorm.h5 cor.correction=amplitude000 cor.updateweights=False',
                    log='$nameMS_diampcor.log', commandType='DP3')

            else:
                logger.info('Solving amp-di (diagonal)...')
                MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.datause=full sol.nchan={nchan_amp} sol.modeldatacolumns=[MODEL_DATA] \
                     sol.mode=diagonal sol.h5parm=$pathMS/amp-di.h5 sol.solint=150 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]',
                     log='$nameMS_diampsol.log', commandType='DP3')

                lib_util.run_losoto(s, 'amp-di', [ms + '/amp-di.h5' for ms in MSs.getListStr()],
                                [f'{parset_dir}/losoto-plot-amp.parset', f'{parset_dir}/losoto-plot-ph.parset', f'{parset_dir}/losoto-amp-di.parset'],
                                plots_dir=f'{plot_dir}/plots-amp-di', h5_dir=sol_dir)

                # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
                logger.info('Correct amp-di (CORRECTED_DATA -> CORRECTED_DATA)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb={sol_dir}/cal-amp-di.h5 cor.correction=amplitudeSmooth cor.updateweights=False',
                    log='$nameMS_diampcor.log', commandType='DP3')

    ########################### IMAGING ####################################
    with w.if_todo('c%02i-corrFR' % c):
        # Correct for FR CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Correcting for FR...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA\
                cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")

    facetregname = f'{sol_dir}/facets-c{c}.reg'
    wide_h5 = f'{sol_dir}/cal-tec-no3c-c{c}.h5' if os.path.exists(f'{sol_dir}/cal-tec-no3c-c{c}.h5') else f'{sol_dir}/cal-tec-c{c}.h5'

    channels_out = MSs.getChout(4.e6) if MSs.getChout(4.e6) > 1 else 2
    with w.if_todo('c%02i-imaging' % c):
        logger.info('Preparing region file...')
        s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+f' --h5 {wide_h5} --imsize {int(1.1*imgsizepix_wide)} \
            --pixelscale {pixscale} --writevoronoipoints --output {facetregname}', log='facet_generator.log', commandType='python')
        s.run(check = True)

        imagename = 'img/wide-' + str(c)
        imagenameM = 'img/wideM-' + str(c)
        # common imaging arguments used by all of the following wsclean calls
        widefield_kwargs = dict(data_column='CORRECTED_DATA', use_shm=use_shm, size=imgsizepix_wide, scale=f'{pixscale}arcsec', weight='briggs -0.5', niter=1000000,
                                gridder='wgridder',  parallel_gridding=len(sm.getPatchNames())*channels_out, minuv_l=30, mgain=0.85, parallel_deconvolution=1024,
                                join_channels='', fit_spectral_pol=3, channels_out=channels_out, deconvolution_channels=3, multiscale='',
                                multiscale_scale_bias=0.65, pol='i', facet_regions=facetregname, apply_facet_solutions=f'{wide_h5} phase000', concat_mss=False)

        # for low-freq or low-dec data, allow the beam to be fitted, otherwise (survey) force 15"
        if not(np.mean(MSs.getFreqs()) < 50e6) and not (phasecentre[1] < 23):
            widefield_kwargs['beam_size'] = 15

        # c0: make quick initial image to get a mask
        if c==0:
            logger.info('Making wide-field image for clean mask...')
            lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename,  no_update_model_required='',
                                 auto_threshold=5.0, auto_mask=8.0, multiscale_max_scales=3, nmiter=6, keep_concat=True, **widefield_kwargs)
            # make initial mask
            current_best_mask = make_current_best_mask(imagename, mask_threshold[c], userReg)
            # safe a bit of time by reusing psf and dirty in first iteration
            reuse_kwargs = {'reuse_concat':True, 'reuse_psf':imagename, 'reuse_dirty':imagename}
        else:
            current_best_mask = f'img/wideM-{c-1}-mask.fits' # is this already set by the make_current_best_mask() below? (not if we restart)
            reuse_kwargs = {}

        # main wsclean call, with mask now
        logger.info('Making wide field image ...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagenameM, fits_mask=current_best_mask,
                             save_source_list='', no_update_model_required='',  nmiter=12,  auto_threshold=2.0, auto_mask=4.0,
                             apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='',
                             local_rms='', local_rms_window=50, local_rms_strength=0.5, **widefield_kwargs, **reuse_kwargs)

        current_best_mask = make_current_best_mask(imagenameM, mask_threshold[c]-0.5, userReg)

        # reset NaNs if present
        im = lib_img.Image(f'{imagenameM}-MFS-image.fits')
        im.nantozeroModel()

    #####################################################################################################
    # Find calibration solutions for subfield - recreate SUBFIELD_DATA

    # User provided subfield
    if subfield:
        if len(Regions.read(subfield)) > 1:
            raise ValueError(f'Manual subfield region {subfield} contains more than one region.')
        os.system(f'cp {subfield} {subfield_path}')

    # Generate a subfield
    if not os.path.exists(subfield_path):
        sm = lsmtool.load(f'img/wideM-{c}-sources.txt')
        # sm.remove('img/wide-lr-mask.fits=1')  # remove sidelobe sources that were subtracted
        sm.remove('MajorAxis > 80')  # remove largest scales
        field_center1, field_size1 = lib_dd.make_subfield_region(subfield_path, MSs.getListObj()[0], sm,
                                                                     subfield_min_flux, pixscale, imgsizepix_wide, debug_dir='img/')

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
                -apply-facet-solutions {wide_h5} phase000 {MSs.getStrWsclean()}',
            log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
        s.run(check=True)

        # Corrupt for FR - MSs: MODEL_DATA -> MODEL_DATA
        logger.info('Corrupting external region model with FR...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA\
                       cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")

        # cycle > 0: need to add DI-corruption on top (previous iteration sub-field) and DI-amp - MSs: MODEL_DATA -> MODEL_DATA
        if c == 0:
            MSs.addcol('SUBFIELD_DATA','PREPARED_DATA')
        else:
            logger.info('Add previous iteration sub-field corruption...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                    cor.parmdb={sol_dir}/cal-tec-sf-c' + str(c-1) + '.h5 cor.correction=phase000 cor.invert=False',
                    log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
            # Corrupt MSs:MODEL_DATA -> MODEL_DATA (sf corr, dieamp corr)
            corr_die_amp(MSs, col='MODEL_DATA', fulljones=fulljones, invert=False)

        # subtract external region MSs: PREPARED_DATA - MODEL_DATA -> SUBFIELD_DATA
        logger.info('Subtracting external region model (SUBFIELD_DATA = PREPARED_DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set SUBFIELD_DATA = PREPARED_DATA - MODEL_DATA"', log='$nameMS_taql.log', commandType='general')
        if develop: clean_empty(MSs, 'empty_sf', 'SUBFIELD_DATA', size=10000)
    ### DONE

    with w.if_todo('c%02i_intreg_predict' % c):
        # Predict internal region - MSs: MODEL_DATA
        logger.info('Predict model of internal region...')
        s.add(f'wsclean -predict -padding 1.8 -name img/wideMint-{c} -j {s.max_cpucores} -channels-out {channels_out} \
                -facet-regions {facetregname}  -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
                {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
        s.run(check=True)

        # Corrupt for FR internal region - MSs: MODEL_DATA -> MODEL_DATA
        logger.info('Corrupting internal region model with FR...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA\
                       cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")
    ### DONE

    if c>0:
        with w.if_todo('c%02i_subfierld_corr_diamp' % c):
            # Correct amp MSs:SUBFIELD_DATA -> SUBFIELD_DATA
            corr_die_amp(MSs, col='SUBFIELD_DATA', fulljones=fulljones)
        ### DONE

    with w.if_todo('c%02i_subfield_solve_tec' % c):
        # Smooth MSs: SUBFIELD_DATA -> SMOOTHED_DATA
        MSs.run_Blsmooth('SUBFIELD_DATA', logstr=f'smooth-c{c}')
        # solve ionosphere phase - ms:SMOOTHED_DATA
        logger.info('Solving ionosphere (subfield)...')
        smMHz_sf = np.array([2.5,4.0,10.0,15.0])
        smMHz_factors_sf = smMHz_sf/np.max(smMHz_sf) # factors should be <1 otherwise trimming of kernel is off, so normalize
        ant_avg_factors = f'"[CS*:15,[RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LBA]:4,[RS208LBA,RS307LBA,RS406LBA,RS407LBA]:2,[RS210LBA,RS310LBA,RS409LBA,RS508LBA,RS509LBA]:1]"'
        ant_smooth_factors = f'"[CS*:{smMHz_factors_sf[3]},[RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LBA]:{smMHz_factors_sf[2]},[RS208LBA,RS307LBA,RS406LBA,RS407LBA]:{smMHz_factors_sf[1]},[RS210LBA,RS310LBA,RS409LBA,RS508LBA,RS509LBA]:{smMHz_factors_sf[0]}]"'

        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS sol.h5parm=$pathMS/tec-sf.h5 sol.solint=15 \
                  sol.mode=scalarphase sol.smoothnessconstraint={max(smMHz_sf)}e6 sol.smoothnessreffrequency=54e6 sol.nchan=1  \
                  sol.modeldatacolumns="[MODEL_DATA]" sol.solutions_per_direction="[15]" \
                  sol.antenna_averaging_factors={ant_avg_factors} sol.antenna_smoothness_factors={ant_smooth_factors}',
                log='$nameMS_solTEC-sf-c' + str(c) + '.log', commandType='DP3', maxProcs=8)

        lib_util.run_losoto(s, f'tec-sf-c{c}',[ms + f'/tec-sf.h5' for ms in MSs.getListStr()],
                            [parset_dir + '/losoto-refph.parset', f'{parset_dir}/losoto-plot-scalarph.parset'],
                            plots_dir=f'{plot_dir}/plots-tec-sf-c{c}', h5_dir=sol_dir)
    ### DONE

    with w.if_todo('c%02i_subfield_corr_tec' % c):
        # Correct MSs: SUBFIELD_DATA -> SUBFIELD_DATA
        logger.info('Correct subfield iono...')
        if sf_phaseSolMode in ['tec', 'tecandphase']:
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                    cor.parmdb={sol_dir}/cal-tec-sf-c' + str(c) + '.h5 cor.correction=tec000 ',
                    log='$nameMS_sf-correct.log', commandType='DP3')
        if sf_phaseSolMode in ['phase', 'tecandphase']:
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                    cor.parmdb={sol_dir}/cal-tec-sf-c' + str(c) + '.h5 cor.correction=phase000',
                    log='$nameMS_sf-correct.log', commandType='DP3')
    ### DONE

    # Do a quick image of the subfield, not strictly necessary but good to have...
    with w.if_todo('c%02i_image-subfield' % c):
        logger.info('Correcting subfield region with FR...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA\
                       cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
        logger.info('Test image subfield...')
        # Note that here the beam is not applied, so the image is suboptimal
        lib_util.run_wsclean(s, 'wscleanSF-c'+str(c)+'.log', MSs.getStrWsclean(), name=f'img/subfield-{c}', data_column='SUBFIELD_DATA', size=int(1.2*subfield_size*3600/pixscale), scale=f'{pixscale}arcsec',
                             weight='briggs -0.5', niter=100000, gridder='wgridder',  parallel_gridding=MSs.getChout(4.e6), shift=f'{subfield_center[0].to(u.hourangle).to_string()} {subfield_center[1].to_string()}',
                             no_update_model_required='', minuv_l=30, beam_size=15, mgain=0.85, nmiter=12, parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                             join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6), deconvolution_channels=3, baseline_averaging='',
                             multiscale='',  multiscale_scale_bias=0.7, pol='i')
    ### DONE

    #####################################################################################################
    # Subtract side-lobe sources (modifies PREPARED_DATA removing sidelobe sources) - recreate SUBFIELD_DATA
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
                   -apply-facet-solutions {wide_h5} phase000 {MSs.getStrWsclean()}',
                log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
            s.run(check=True)
            # Corrupt for FR - MSs: MODEL_DATA -> MODEL_DATA
            logger.info('Corrupting mainlobe region model with FR...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA\
                       cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")

            # subtract internal region from MSs: PREPARED_DATA - MODEL_DATA -> SUBFIELD_DATA
            logger.info('Subtract main-lobe (SUBFIELD_DATA = PREPARED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = PREPARED_DATA - MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')
            if develop: clean_empty(MSs, 'only_sidelobe', 'SUBFIELD_DATA', size=10000)
        ### DONE

        # Do a rough correction of the sidelobe data using the subfield solutions and FR
        # MSs: SUBFIELD_DATA -> SUBFIELD_DATA
        with w.if_todo('correct-sidelobe'): # just for testing/debug
            logger.info('Correct sidelobe region with subfield solutions...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA \
                    cor.parmdb={sol_dir}/cal-tec-sf-c' + str(c) + '.h5 cor.correction=phase000',
                    log='$nameMS_sf-correct.log', commandType='DP3')
            # Corrected for FR - MSs: SUBFIELD_DATA -> SUBFIELD_DATA
            logger.info('Correct sidelobe region model with FR...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=SUBFIELD_DATA msout.datacolumn=SUBFIELD_DATA\
                       cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")
        ### DONE

        imagename_lr = 'img/wide-lr'
        channels_out_lr = MSs.getChout(2.e6) if MSs.getChout(2.e6) > 1 else 2
        # Image the sidelobe data and predict model
        # MSs: create MODEL_DATA (with just the sidelobe flux)
        # TODO CHECK: should we do the clean + predict with beam?
        facetregname_lr = f'{sol_dir}/facets-lr-c{c}.reg'
        with w.if_todo('image_lr'):
            logger.info('Preparing region file...')
            # this is too slow
            #s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+f' --grid 10 --imsize {int(1.1*imgsizepix_lr)} \
            #--pixelscale 30 --writevoronoipoints --output {facetregname_lr}', log='facet_generator.log', commandType='python')
            #s.run(check = True)
            #facet_regions=facetregname_lr, apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='', 
            # --facet-regions {facetregname_lr} -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam
            logger.info('Cleaning sidelobe low-res...')
            lib_util.run_wsclean(s, 'wscleanLR.log', MSs.getStrWsclean(), name=imagename_lr, data_column='SUBFIELD_DATA',
                                size=imgsizepix_lr, scale='30arcsec', save_source_list='',  parallel_gridding=channels_out_lr, baseline_averaging='',
                                weight='briggs -0.5', niter=50000, no_update_model_required='', minuv_l=30, maxuvw_m=6000,
                                taper_gaussian='200arcsec', mgain=0.85, channels_out=channels_out_lr, parallel_deconvolution=512,
                                local_rms='', auto_mask=3, auto_threshold=1.5, join_channels='', fit_spectral_pol=5,
                                do_predict=True)
            #s.add(f'wsclean -predict -padding 1.8 -name {imagename_lr} -j {s.max_cpucores} -channels-out {channels_out_lr} \
            #          {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
            #s.run(check=True)
        ### DONE

        # Subtract full low-resolution field (including possible large-scale emission within primary beam)
        # to get empty data set for flagging
        # MSs: SUBFIELD_DATA - MODEL_DATA -> SUBFIELD_DATA
        with w.if_todo('subtract_lr'):
            logger.info('Subtract low-resolution to get empty data set (SUBFIELD_DATA = SUBFIELD_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set SUBFIELD_DATA = SUBFIELD_DATA - MODEL_DATA"',
                log='$nameMS_taql.log', commandType='general')
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
                for im in glob.glob(f'{imagename_lr}-0*model*fits'):
                    wideLRext = im.replace(imagename_lr, f'{imagename_lr}-blank')
                    os.system('cp %s %s' % (im, wideLRext))
                    lib_img.blank_image_reg(wideLRext, beamReg , blankval=0.)

                logger.info('Predict model of sidelobe region (wsclean)...')
                # --facet-regions {facetregname_lr} -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
                s.add(f'wsclean -predict -padding 1.8 -name {imagename_lr}-blank -j {s.max_cpucores} -channels-out {channels_out_lr} \
                      {MSs.getStrWsclean()}', log='wscleanPRE-c' + str(c) + '.log', commandType='wsclean')
                s.run(check=True)

            elif sidelobe_predict_mode=='DP3':
                logger.info('Predict model of sidelobe region (DP3)...')
                sidelobe_sky = lsmtool.load(f'{imagename_lr}-sources.txt')
                sidelobe_sky.remove(f'{beamMask}==True')
                MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={imagename_lr}-sources.txt msout.datacolumn=MODEL_DATA',
                    log='$nameMS_pre.log', commandType='DP3')
                
            logger.info('Corrupt sidelobe model with subfield solutions...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                    cor.parmdb={sol_dir}/cal-tec-sf-c' + str(c) + '.h5 cor.correction=phase000 cor.invert=False',
                    log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
            logger.info('Corrupt sidelobe  model with FR...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA\
                       cor.correction=rotationmeasure000 cor.invert=False', log='$nameMS_corFR.log', commandType="DP3")

            logger.info('Subtract corrupted sidelobe model (PREPARED_DATA = PREPARED_DATA - MODEL_DATA)...')
            MSs.run('taql "update $pathMS set PREPARED_DATA = PREPARED_DATA - MODEL_DATA"', log='$nameMS_taql.log', commandType='general')
        ### DONE

    # apply the subfield solutions to the data.
    with w.if_todo('c%02i_corr_sf_sols' % c):
        logger.info('Correct subfield ionosphere (PREPARED_DATA -> CORRECTED_DATA)...')
        # Correct MSs:PREPARED_DATA -> CORRECTED_DATA (sf corr)
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=PREPARED_DATA msout.datacolumn=CORRECTED_DATA  \
                cor.parmdb={sol_dir}/cal-tec-sf-c' + str(c) + '.h5 cor.correction=phase000',
                log='$nameMS_sf-correct.log', commandType='DP3')
    ### DONE

    # finally re-correct for die amp on the newly created CORRECTED_DATA
    if c >= 1:
        with w.if_todo('c%02i_corr_dieamp' % c):
            # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA (sf corr, dieamp corr)
            corr_die_amp(MSs, col='CORRECTED_DATA', fulljones=fulljones)
        ### DONE

with w.if_todo('final_fr_corr'):
    # Before leaving, apply also FR - MSs: CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Final FR correction...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={sol_dir}/cal-fr.h5 msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA\
            cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType="DP3")

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image*.fits ddparallel/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-MFS-residual*.fits ddparallel/images') for c in range(maxIter) ]
[ os.system('mv img/wideM-'+str(c)+'-sources*.txt ddparallel/images') for c in range(maxIter) ]
[ os.system('mv img/subfield-'+str(c)+'-MFS-image*.fits ddparallel/images') for c in range(maxIter) ]
os.system('mv img/wide-lr-MFS-image.fits ddparallel/images')

# debug images
if develop:
    os.system('mv img/only*image.fits ddparallel/images')
    os.system('mv img/empty*image.fits ddparallel/images')
else:
    # remove MODEL_DATA cols...
    logger.info('Removing unwanted columns...')
    for patch in patches:
        MSs.deletecol(patch)

# Copy model
os.system(f'mv img/wideM-{maxIter-1}-*-model.fits ddparallel/skymodel')
os.system(f'mv img/wideM-{maxIter-1}-*-model-fpb.fits ddparallel/skymodel')
os.system(f'mv img/wideM-{maxIter-1}-*-model-pb.fits ddparallel/skymodel')

w.alldone()
