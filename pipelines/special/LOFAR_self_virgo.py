#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# they need to be in "./mss/"

import sys, os, glob
import numpy as np
from regions import Regions
from astropy.io import fits
from astropy.wcs import WCS
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
subfield_min_flux = parset.getfloat('LOFAR_self','subfield_min_flux') # default 40 Jy
subfield = parset.get('LOFAR_self','subfield') # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
maxIter = parset.getint('LOFAR_self','maxIter') # default = 2 (try also 3)
phaseSolMode = parset.get('LOFAR_self', 'ph_sol_mode') # tecandphase, tec, phase
sourcedb = parset.get('model','sourcedb')
apparent = parset.getboolean('model','apparent')
userReg = parset.get('model','userReg')

spec = lambda nu: 10**(2.4466 - 0.8116 * ((np.log10(nu/1.e9))**1) - 0.0483 * ((np.log10(nu/1.e9))**2) ) # PB17

#############################################################################

def testimage(c, MSs):
    imagename = 'img/test-' + str(c)
    lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(),
                         name=imagename, size=2000, scale='4arcsec', baseline_averaging='',
                         weight='briggs -0.3', niter=1000000, no_update_model_required='', minuv_l=30,
                         parallel_gridding=6,  mgain=0.8, parallel_deconvolution=1024,
                         auto_mask=5, auto_threshold=3.5, join_channels='',
                         channels_out=6, multiscale='')


def clean_self(c, MSs):
    imagenameM = 'img/wideM-' + str(c)
    lib_util.run_wsclean(s, 'wscleanM-c' + str(c) + '.log', MSs.getStrWsclean(), concat_mss=True, fits_mask='mask-image.fits',
                         name=imagenameM, reuse_psf=imagenameM, reuse_dirty=imagenameM, save_source_list='', size=1500,
                         scale='4arcsec', weight='briggs -0.3', niter=1000000, no_update_model_required='', minuv_l=30,
                         parallel_gridding=6,  mgain=0.7,  auto_threshold=2., join_channels='', fit_spectral_pol=5, channels_out=MSs.getChout(4.e6),
                         multiscale='')

    os.system('cat ' + logger_obj.log_dir + '/wscleanM-c' + str(c) + '.log | grep "background noise"')

    img = lib_img.Image(imagenameM)
    img.rescaleModel(spec)

    logger.info('Predict model...')
    pred_str = f'wsclean -predict -padding 1.8 -name img/wideM-{c} -j {s.max_processors} -channels-out {MSs.getChout(4e6)}'
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

if not os.path.exists('template-image.fits'):
    logger.info('Create mask...')
    # dummy clean to get image -> mask
    lib_util.run_wsclean(s, 'wsclean-template.log', MSs.getStrWsclean(), niter=0, channel_range='0 1',
                         no_reorder='',
                         interval='0 10', name='template', scale='4asec', size=1500,
                         nmiter=0)
    # create basemask
    lib_img.blank_image_reg('template-image.fits', 'm87_flux.reg', inverse=True, blankval=0.)
    lib_img.blank_image_reg('template-image.fits', 'm87_flux.reg', inverse=False, blankval=1.)
#################################################################
# Get online model
# if sourcedb == '':
#     if not os.path.exists('tgts.skymodel'):
#         fwhm = MSs.getListObj()[0].getFWHM(freq='min')
#         radeg = phasecentre[0]
#         decdeg = phasecentre[1]
#         # get model the size of the image (radius=fwhm/2)
#         os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
#         lsm = lsmtool.load('tgts.skymodel', beamMS=MSs.getListStr()[0])#, beamMS=MSs.getListObj()[0])
#         lsm.remove('I<1')
#         lsm.write('tgts-beam.skymodel', applyBeam=True, clobber=True)
#         lsm.write('tgts.skymodel', applyBeam=False, clobber=True)
#         apparent = False
#     sourcedb = 'tgts.skymodel'
#################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS + '/' + sourcedb_basename)
    logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
    os.system('cp -r ' + sourcedb + ' ' + MS)

# Here the model is added only to CS+RS, IS used only for FR and model is not needed
with w.if_todo('init_model'):
    # model = '/beegfs/p1uy068/virgo/models/LBA/fits-models-freqcut/Vir'
    # n = len(glob.glob(f'{model}-[0-9]*-model.fits'))
    # logger.info('Predict (wsclean: %s - chan: %i)...' % (model, n))
    # s.add(
    #     f'wsclean -predict -no-reorder -name {model} -j {s.max_processors} -use-wgridder -channels-out {n} {MSs.getStrWsclean()}',
    #     log='wscleanPRE-init.log', commandType='wsclean', processors='max')
    # s.run(check=True)
    logger.info('Add model to MODEL_DATA...')
    if apparent:
        MSs.run('DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=false pre.sourcedb=$pathMS/' + sourcedb_basename,
            log='$nameMS_pre.log', commandType='DP3')
    else:
        MSs.run('DP3 ' + parset_dir + '/DP3-predict.parset msin=$pathMS pre.usebeammodel=true pre.usechannelfreq=True \
                 pre.beammode=array_factor pre.onebeamperpatch=True pre.sourcedb=$pathMS/' + sourcedb_basename,
                log='$nameMS_pre.log', commandType='DP3')
# ### DONE

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
    MSs.addcol('FR_MODEL_DATA', 'MODEL_DATA', usedysco=False)
    MSs.run('taql "UPDATE $pathMS SET FR_MODEL_DATA[,0]=0.5+0i, FR_MODEL_DATA[,1]=0.0+0i, FR_MODEL_DATA[,2]=0.0+0i, \
         FR_MODEL_DATA[,3]=0.5+0i"', log='$nameMS_taql_frmodel.log', commandType='general')

    # Solve cal_SB.MS:CIRC_PHASEDIFF_DATA against FR_MODEL_DATA (only solve - solint=2m nchan=0 as it has the smoothnessconstrain)
    logger.info('Solving circ phase difference ...')
    MSs.run('DP3 ' + parset_dir + '/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.solint=' + str(30 * base_solint) + ' ',
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
    with w.if_todo('set_corrected_data_c%02i' % c):
        logger.info('Creating CORRECTED_DATA from FR_CORRECTED_DATA...')
        MSs.addcol('CORRECTED_DATA', 'FR_CORRECTED_DATA')
        MSs.run('taql UPDATE $pathMS SET CORRECTED_DATA=FR_CORRECTED_DATA', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

    with w.if_todo('solve_tec1_c%02i' % c):
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')

        # solve ionosphere phase - ms:SMOOTHED_DATA (1m 2SB)
        logger.info('Solving TEC1...')
        if phaseSolMode == 'phase': #phase
            solver_params = f'sol.mode=scalarphase sol.smoothnessconstraint=0.5e6 sol.smoothnessreffrequency=54e6'
            losoto_parsets = [parset_dir+'/losoto-plot-scalar.parset']
        else: # TEC or TecAndPhase
            solver_params = f'sol.mode={phaseSolMode} sol.approximatetec=True sol.maxapproxiter=250 sol.approxtolerance=1e-3'
            losoto_parsets = [parset_dir+'/losoto-plot-tec.parset']

        MSs.run(f"DP3 {parset_dir}/DP3-solTEC.parset msin=$pathMS sol.h5parm=$pathMS/tec1.h5 sol.solint={base_solint} {solver_params}",
                log='$nameMS_solTEC-c'+str(c)+'.log', commandType='DP3')

        lib_util.run_losoto(s, 'tec1-c'+str(c), [ms+'/tec1.h5' for ms in MSs.getListStr()], losoto_parsets)
        os.system('mv cal-tec1-c'+str(c)+'.h5 self/solutions/')
        os.system('mv plots-tec1-c'+str(c)+' self/plots/')
    ### DONE

    with w.if_todo('cor_tec1_c%02i' % c):
        # correct TEC - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        if phaseSolMode in ['tec', 'tecandphase']:
            logger.info('Correcting TEC1...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA\
                    cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=tec000',
                    log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DP3')
        if phaseSolMode in ['phase', 'tecandphase']:
            logger.info('Correcting ph1...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA\
                    cor.parmdb=self/solutions/cal-tec1-c'+str(c)+'.h5 cor.correction=phase000',
                    log='$nameMS_corTEC-c'+str(c)+'.log', commandType='DP3')
        # testimage('iono', MSs)
    ### DONE

    with w.if_todo('solve_fast_g_c%02i' % c):
        # DIE Calibration - ms:CORRECTED_DATA (8m, 4SB)
        logger.info('Solving fast G ...')
        # MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
        MSs.run(f'DP3 {parset_dir}/DP3-solGfj.parset msin=$pathMS sol.h5parm=$pathMS/g1.h5 sol.solint={2*base_solint} sol.nchan=1 sol.smoothnessconstraint=4e6 sol.mode=scalarcomplexgain',
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DP3')
        lib_util.run_losoto(s, 'g1-c'+str(c), [MS+'/g1.h5' for MS in MSs.getListStr()],
                            [parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-scalar.parset'])#, parset_dir+'/losoto-bp.parset'])
        os.system('mv plots-g1-c'+str(c)+' self/plots/')
        os.system('mv cal-g1-c'+str(c)+'.h5 self/solutions/')
    ### DONE

    with w.if_todo('cor_fast_g_c%02i' % c):
        # correct G - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting G...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                cor.parmdb=self/solutions/cal-g1-c'+str(c)+'.h5 cor.correction=amplitude000',
                log='$nameMS_corG-c'+str(c)+'.log', commandType='DP3')
        # testimage('amp', MSs)
        logger.info('Correcting gph...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                cor.parmdb=self/solutions/cal-g1-c'+str(c)+'.h5 cor.correction=phase000',
                log='$nameMS_corG-c'+str(c)+'.log', commandType='DP3')
        # testimage('amp+ph', MSs)


    with w.if_todo('solve_slow_g_c%02i' % c):
        # DIE Calibration - ms:CORRECTED_DATA (8m, 4SB)
        logger.info('Solving slow G (full jones)...')
        # MSs.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-c{c}')
        MSs.run('DP3 '+parset_dir+'/DP3-solGfj.parset msin=$pathMS sol.h5parm=$pathMS/g2.h5 sol.solint='+str(120*base_solint)+' sol.nchan='+str(16*base_nchan),
                log='$nameMS_solG-c'+str(c)+'.log', commandType='DP3')
        lib_util.run_losoto(s, 'g2-c'+str(c), [MS+'/g2.h5' for MS in MSs.getListStr()],
                [parset_dir+'/losoto-plot-fullj.parset'])#, parset_dir+'/losoto-bp.parset'])
        os.system('mv plots-g2-c'+str(c)+' self/plots/')
        os.system('mv cal-g2-c'+str(c)+'.h5 self/solutions/')
    ### DONE

    with w.if_todo('cor_g_c%02i' % c):
        # correct G - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        logger.info('Correcting G...')
        MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                cor.parmdb=self/solutions/cal-g2-c'+str(c)+'.h5 cor.correction=fulljones cor.soltab=[amplitude000,phase000]',
                log='$nameMS_corG-c'+str(c)+'.log', commandType='DP3')
        # testimage('fj', MSs)

    with w.if_todo(f'clean_c{c}'):
        clean_self(c, MSs)

    with w.if_todo(f'flag_residuals'):
        logger.info('SET CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA...')
        MSs.run('taql UPDATE $pathMS SET CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        logger.info('Flagging residuals...')
        MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.lua',
                log='$nameMS_flag-c'+str(c)+'.log', commandType='DP3')

    sys.exit()