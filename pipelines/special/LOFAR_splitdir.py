#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# HE: this is an experimental script to split out individual directions from LBA observations with IS present
#     use cases are e.g. the preparation of in-field calibrators

import sys, os, glob, argparse
import numpy as np
import lsmtool as lsm
import astropy.wcs
import warnings
from losoto.h5parm import h5parm
from LiLF import lib_ms, lib_img, lib_util, lib_log

#####################################################
def test_image_dutch(MSs, imgname):
    """ Create a quick debug image..."""
    lib_util.run_wsclean(s, 'wsclean-test.log', MSs.getStrWsclean(), name=f'img/{imgname}',
                         data_column='SUBTRACTED_DATA', size=3000, scale=f'4arcsec',
                         weight='briggs -0.3', niter=100000, gridder='wgridder', parallel_gridding=6,
                         no_update_model_required='', minuv_l=30, maxuvw_m=max_uvw_m_dutch, mgain=0.85, nmiter=12,
                         parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                         join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6),
                         deconvolution_channels=3, baseline_averaging='',
                         multiscale='', multiscale_scale_bias=0.7, pol='i')

parser = argparse.ArgumentParser(description='Split out a single direction by subtracting the rest field and correcting the stations.')
parser.add_argument('--dirreg',action='store', default=None, type=str, help='Provide a region for the source to be split of (ds9 circle or square).')
parser.add_argument('--name', dest='name', type=str, default='target_extracted', help='Name of the direction.')
parser.add_argument('--dutchdir', type=str, default='dutchdir', help='Directory of the dutch processing.')
parser.add_argument('--freqres', type=float, default=0.195312, help='Freq. resolution of the split-off MSs in Mhz. Default=0.195312MHz (1 subband)')
parser.add_argument('--timeres', type=int, default=16, help='Time resolution of the split-off MSs in s. Default: 16s.')

args = parser.parse_args()
regfile, name, dutchdir = args.dirreg, args.name, args.dutchdir
time_resolution, freq_resolution = args.timeres, args.freqres

logger_obj = lib_log.Logger(f'pipeline-splitdir')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker(f'pipeline-splitdir.walker')
warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_splitdir','parset_dir')

if not regfile:
    logger.error('No target direction region file provided (--dirreg.')
    sys.exit()

if not os.path.exists(f'mss-{name}'):
    os.makedirs(f'mss-{name}')

if not os.path.exists(f'img'):
    os.makedirs(f'img')

dir_reg = lib_util.Region_helper(regfile)
center = dir_reg.get_center() # center of the extract region

MSs = lib_ms.AllMSs(glob.glob('mss-IS/*MS'), s, check_flags=False, check_sun=False)
max_uvw_m_dutch = MSs.getMaxBL(check_flags=True, dutch_only=True) # this is important, it is the maximum uvw value in meters of any dutch-dutch baseline. Everything above this value is certainly IS data

# 1. Apply dutch direction-independent solutions (FR, amp and phase)
with w.if_todo('correct_dutch_di'):
    logger.info('Correcting FR (Dutch stations) DATA -> CORRECTED_DATA...')
    # Correct FR with results of solve - TC.MS: DATA -> CORRECTED_DATA
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={dutchdir}/ddparallel/solutions/cal-fr.h5 \
            cor.correction=rotationmeasure000', log='$nameMS_corFR.log', commandType='DP3')
    logger.info('Correcting subfield phase (Dutch stations) CORRECTED_DATA -> CORRECTED_DATA...')
    # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
            cor.parmdb={dutchdir}/ddparallel/solutions/cal-tec-sf-merged-c1.h5 cor.correction=phase000',
            log='$nameMS_sf-correct.log', commandType='DP3')
    # Correct MSs:CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Correcting DI amplitude (Dutch stations) CORRECTED_DATA -> CORRECTED_DATA...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                cor.parmdb={dutchdir}/ddparallel/solutions/cal-amp-di.h5 cor.correction=amplitudeSmooth ',
            log='$nameMS_sidelobe_corrupt.log', commandType='DP3')
    logger.info('Test image (corrected)...')
    lib_util.run_wsclean(s, 'wsclean-test0.log', MSs.getStrWsclean(), name=f'img/dutchdicorr',
                         data_column='CORRECTED_DATA', size=3000, scale=f'4arcsec',
                         weight='briggs -0.3', niter=100000, gridder='wgridder', parallel_gridding=6,
                         no_update_model_required='', minuv_l=30, maxuvw_m=max_uvw_m_dutch, mgain=0.85, nmiter=12,
                         parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                         join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6),
                         deconvolution_channels=3, baseline_averaging='',
                         multiscale='', multiscale_scale_bias=0.7, pol='i')
### DONE

# 2. Prepare the h5parm solutions from dutch processing so they can be used for IS data (reorder dirs and add IS stations with unit solutions)
with w.if_todo('interph5'):
    # first we need to reorder the soltab dir axis to have the same order as the solset.getSou() dict, otherwise h5_merger creates a mess (best would be to fix this in h5_merger)
    os.system(f'cp {dutchdir}/ddserial/c00/solutions/interp.h5 interp.h5')
    with h5parm('interp.h5', readonly=False) as h5:
        solset = h5.getSolset('sol000')
        soltab_ph = solset.getSoltab('phase000')
        soltab_amp = solset.getSoltab('amplitude000')

        sou = solset.getSou()
        order_ph = []
        order_amp = []
        for src in sou:
            order_ph.append(np.argwhere(soltab_ph.dir == src)[0])
            order_amp.append(np.argwhere(soltab_amp.dir == src)[0])
        order_ph = np.squeeze(order_ph)
        order_amp = np.squeeze(order_amp)
        h5.getSolset('sol000').getSoltab('phase000').setValues(soltab_ph.getValues()[0][order_ph])
        h5.getSolset('sol000').getSoltab('phase000').setAxisValues('dir', list(sou.keys()))
        h5.getSolset('sol000').getSoltab('amplitude000').setValues(soltab_amp.getValues()[0][order_amp])
        h5.getSolset('sol000').getSoltab('amplitude000').setAxisValues('dir', list(sou.keys()))
    s.add(f'h5_merger.py --h5_out interp_merged.h5 --h5_tables interp.h5 -ms {MSs.getListStr()[0]} --add_ms_stations --no_antenna_crash --propagate_flags')
    s.run(check=True)

# 3. Predict corrupted visibilities for all but the direction to split off - set to zero for all non-dutch baselines!
with w.if_todo('predict'):
    # prepare model of central/external regions
    logger.info('Blanking direction region of model files and reverse...')
    for im in glob.glob(f'{dutchdir}/ddserial/c00/images/wideDD-*model*.fits'):
        wideMext = im.replace('wideDD-c00',f'wideDD-c00-{name}').split('/')[-1]
        os.system('cp %s %s' % (im, wideMext))
        lib_img.blank_image_reg(wideMext, regfile, blankval = 0.)
    # Recreate MODEL_DATA of external region for subtraction
    logger.info('Predict corrupted model of external region (wsclean)...')
    s.add(f'wsclean -predict -padding 1.8 -name wideDD-c00-{name} -j {s.max_cpucores} -channels-out {len(glob.glob(f"wideDD-c00-{name}*fpb.fits"))} \
            -facet-regions {dutchdir}/ddserial/c00/images/wideDD-c00_facets.reg -maxuvw-m {max_uvw_m_dutch} -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
            -apply-facet-solutions interp_merged.h5 phase000,amplitude000 {MSs.getStrWsclean()}',
          log='wscleanPRE.log', commandType='wsclean')
    s.run(check=True)
    # Set to zero for non-dutch baselines
    MSs.run("taql 'update $pathMS set MODEL_DATA=0 WHERE ANTENNA1 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME !~p/[CR]S*/]'", log='$nameMS_resetISmodel.log', commandType='general')
    MSs.run("taql 'update $pathMS set MODEL_DATA=0 WHERE ANTENNA2 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME !~p/[CR]S*/]'", log='$nameMS_resetISmodel.log', commandType='general')
    # MSs.run('taql "update $pathMS set MODEL_DATA=0 WHERE ANTENNA1 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME ! \
    #                      ~p/[CR]S*/] && ANTENNA2 in [SELECT ROWID() FROM ::ANTENNA WHERE NAME ! ~p/[CR]S*/]"', log='$nameMS_resetISmodel.log', commandType='general')

# 4. subtract for dutch baselines
with w.if_todo('subtract'):
    MSs.addcol('SUBTRACTED_DATA', 'CORRECTED_DATA')
    logger.info('Subtracting external region model (SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
            log='$nameMS_taql.log', commandType='general')
    test_image_dutch(MSs, 'dutchsub')
    lib_util.run_wsclean(s, 'wsclean-test0.log', MSs.getStrWsclean(), name=f'img/allsub',
                         data_column='SUBTRACTED_DATA', size=6000, scale=f'0.5arcsec',
                         weight='briggs -0.3', niter=100000, gridder='wgridder', parallel_gridding=6,
                         no_update_model_required='', minuv_l=30, mgain=0.85, nmiter=12,
                         parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                         join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6),
                         deconvolution_channels=3, baseline_averaging='', multiscale='', multiscale_scale_bias=0.7, pol='i')

# 5. apply closest direction solutions for dutch baselines
with w.if_todo('correct_dutch_dd'):
    # apply init - closest DDE sol
    # TODO: this assumes phase000 and optionally, amplitude000
    h5init = h5parm('interp_merged.h5')
    solset_dde = h5init.getSolset('sol000')
    # get closest dir to target reg center
    dirs = np.array([solset_dde.getSou()[k] for k in solset_dde.getSoltab('phase000').dir])
    dir_dist = lib_util.distanceOnSphere(dirs[:, 0], dirs[:, 1], *np.deg2rad(center), rad=True)
    closest = solset_dde.getSoltab('phase000').dir[np.argmin(dir_dist)]
    logger.info('Init apply: correct closest DDE solutions ({})'.format(closest))
    logger.info('Correct init ph...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA '
                                          'msout.datacolumn=SUBTRACTED_DATA cor.parmdb=interp_merged.h5 cor.correction=phase000 cor.direction=' + closest,
                    log='$nameMS_init-correct.log', commandType='DP3')
    if 'amplitude000' in solset_dde.getSoltabNames():
        logger.info('Correct init amp...')
        MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA msout.datacolumn=SUBTRACTED_DATA \
                         cor.parmdb=interp_merged.h5 cor.correction=amplitude000 cor.direction=' + closest,
                        log='$nameMS_init-correct.log', commandType='DP3')
    h5init.close()
    test_image_dutch(MSs, 'dutchsubcorr')
### DONE

# 6. phase shift and average the data -> to 16s,
with w.if_todo('phaseshift'):
    t_avg_factor = int(round(time_resolution/MSs.getListObj()[0].getTimeInt()))
    f_avg_factor = int(round(freq_resolution*1e6/MSs.getListObj()[0].getChanband()))
    logger.info(f'Phase shift and avg to {time_resolution}s, {freq_resolution:.4f}MHz (x{t_avg_factor} in t; x{f_avg_factor} in f)...')
    MSs.run(f'DP3 {parset_dir}/DP3-shiftavg.parset msin=$pathMS msout=mss-{name}/$nameMS.MS msin.datacolumn=SUBTRACTED_DATA '
            f'shift.phasecenter=[{center[0]}deg,{center[1]}deg] avg.freqstep={f_avg_factor} avg.timestep={t_avg_factor}',
            log='$nameMS_shiftavg.log', commandType='DP3')
    test_image_dutch(MSs, 'dutchsubcorrshift')

MSs_extract = lib_ms.AllMSs( glob.glob(f'mss-{name}/*.MS'), s )

with w.if_todo('beamcorr'):
    logger.info('Correcting beam...')
    MSs_extract.run('DP3 ' + parset_dir + '/DP3-beam.parset msin=$pathMS', log='$nameMS_beam.log', commandType='DP3')
    test_image_dutch(MSs, 'dutchsubcorrshiftbeam')

sys.exit()

for p in close_pointings:
    MSs = lib_ms.AllMSs( glob.glob('mss-extract/'+p+'/*MS'), s )
    ch_out = MSs.getChout(4e6)  # chout from dd
    fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
    phase_center = MSs.getListObj()[0].getPhaseCentre()
    # read image, h5parm, make mask
    wideDD_image = lib_img.Image(f'extract/init/{p}/wideDD-{highest_ddcal}-MFS-image-pb.fits')
    dde_h5parm = 'extract/init/'+p+'/interp.h5'
    # make mask for subtraction
    mask_ddcal = wideDD_image.imagename.replace('.fits', '_mask-ddcal.fits')  # this is used to find calibrators
    wideDD_image.makeMask(threshpix=5, atrous_do=True, maskname=mask_ddcal, write_srl=True, write_ds9=True)

    # Delete old columns to avoid dysco issues
    with w.if_todo('remove_columns_' + p):
        logger.info('Removing old MODEL_DATA and SUBTRACTED_DATA columns...')
        datadir = os.listdir(f'mss-extract/{p}/')
        for dir in datadir:
            MSs.run(f'taql "ALTER TABLE mss-extract/{p}/{dir} DELETE COLUMN MODEL_DATA, SUBTRACTED_DATA"', log=f'{dir}_deloldcols.log', commandType='python')

    with w.if_todo('predict_rest_' + p):

        # # Add mock MODEL column to avoid DDFacet overflow
        MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA', log='$nameMS_addmodelcol.log', commandType='python')

        # Predict+corrupt in MODEL_DATA of everything BUT the calibrator
        inmask = f'extract/init/{p}/wideDD-{highest_ddcal}-MFS-image-pb_mask-ddcal.fits'
        outmask = inmask + '.mask'
        lib_img.blank_image_reg(inmask, target_reg_file, outfile=outmask, inverse=False, blankval=0.)

        for im in glob.glob(f'extract/init/{p}/wideDD-{highest_ddcal}-0*model-pb.fits'):
            wideDDext = im.replace('wideDD', 'wideDDext')
            os.system('cp %s %s' % (im, wideDDext))
            lib_img.blank_image_reg(wideDDext, target_reg_file, blankval=0.)

        # # if we have subtract reg, unmask that part again to predict+subtract it.
        if subtract_reg_file:
            os.system(f'cp ../{subtract_reg_file} .')
            logger.info(f"Re-adding sources in subtract-region {subtract_reg_file} to subtraction model.")
            lib_img.blank_image_reg(outmask, subtract_reg_file, inverse=False, blankval=1.)

        h5init = h5parm(dde_h5parm)
        solset_dde = h5init.getSolset('sol000')

        if 'amplitude000' in solset_dde.getSoltabNames():
            correct_for = 'phase000, amplitude000'
        else:
            correct_for = 'phase000'

        facet_path = f'extract/init/{p}/wideDD-{highest_ddcal}_facets.reg'
        #TODO apply beam here when bug is fixed...
        s.add(f'wsclean -predict -padding 1.8 -name extract/init/{p}/wideDDext-{highest_ddcal} -j ' + str(s.max_cpucores) + ' -channels-out ' + str(
            ch_out) + ' -facet-regions ' + facet_path + ' -diagonal-solutions -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
            -apply-facet-solutions ' + dde_h5parm + ' ' + correct_for + ' \
            -reorder -parallel-reordering 4 ' + MSs.getStrWsclean(),
              log='wscleanPRE.log', commandType='wsclean', processors='max')
        s.run(check=True)


if no_selfcal: # finish here
    logger.info('No selfcal option is set in lilf.config. Done.')
    sys.exit(0)


# initial imaging to get the model in the MODEL_DATA (could also be done using the Dico DDFacet model)
do_beam = len(close_pointings) > 1 # if > 1 pointing, correct beam every cycle, otherwise only at the end.
with w.if_todo('image_init'):
    logger.info('Initial imaging...')
    clean('init', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=do_beam, userReg=userReg, datacol='CORRECTED_DATA')

# Smoothing - ms:DATA -> ms:SMOOTHED_DATA
with w.if_todo('smooth'):
    logger.info('BL-based smoothing...')
    MSs_extract.run('BLsmooth.py -c 1 -n 8 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python', maxThreads=1)

# get initial noise and set iterators for timeint solutions
image = lib_img.Image('img/extractM-init-MFS-image.fits', userReg=userReg)
rms_noise_pre, mm_ratio_pre = image.getNoise(), image.getMaxMinRatio()
rms_noise_init, mm_ratio_init = rms_noise_pre, mm_ratio_pre
doamp = False

image.selectCC(checkBeam=False)
sm = lsm.load(image.skymodel_cut)
total_flux = np.sum(sm.getColValues('I', aggregate='sum'))

# Per default we have 32s exposure
# shortest time interval for phase solutions as a function of total flux
ph_int = [8]
if total_flux > 3:
    ph_int.append(4)
if total_flux > 5:
    ph_int.append(2)
if total_flux > 10:
    ph_int.append(1)
iter_ph_solint = lib_util.Sol_iterator(ph_int)
iter_amp_solint = lib_util.Sol_iterator([60, 30, 15])
iter_amp2_solint = lib_util.Sol_iterator([60, 30])
logger.info('RMS noise (init): %f' % (rms_noise_pre))
logger.info('MM ratio (init): %f' % (mm_ratio_pre))
logger.info('Total flux (init): %f Jy' % (total_flux))
rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)

    h5ph = 'extract/cal-ph-c%02i.h5' % c
    solint_ph = next(iter_ph_solint)
    if doamp:
        h5amp1 = 'extract/cal-amp1-c%02i.h5' % c
        solint_amp = next(iter_amp_solint)
        h5amp2 = 'extract/cal-amp2-c%02i.h5' % c
        solint_amp2 = next(iter_amp2_solint)
        h5fj = 'extract/cal-fulljones-c%02i.h5' % c

    if phSolMode == 'phase':
        logger.info('Phase calibration...')
        with w.if_todo('cal-ph-c%02i' % c):
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
                     sol.mode=scalarphase sol.solint=' + str(solint_ph) + ' sol.smoothnessconstraint=5e6 sol.smoothnessreffrequency=54e6 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]',
                     log='$nameMS_solGph-c%02i.log' % c, commandType='DP3')
            lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs_extract.getListStr()],
                                [parset_dir + '/losoto-plot1.parset'],
                                plots_dir='extract/plots-%s' % c)
            os.system('mv cal-ph.h5 %s' % h5ph)

        with w.if_todo('cor-ph-c%02i' % c):
            # correct ph - ms:DATA -> ms:CORRECTED_DATA
            logger.info('Correct ph...')
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=phase000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DP3')
    else:
        logger.info('Tecandphase calibration...')
        with w.if_todo('cal-tecandph-c%02i' % c):
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-soltecandphase.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
                     sol.mode=tecandphase sol.solint=' + str(solint_ph) + ' sol.nchan=1 sol.smoothnessconstraint=5e6 sol.smoothnessreffrequency=54e6 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]',
                     log='$nameMS_soltecandphase-c%02i.log' % c, commandType='DP3')
            lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs_extract.getListStr()],
                                [parset_dir + '/losoto-plottecandphase.parset'],
                                plots_dir='extract/plots-%s' % c)
            os.system('mv cal-ph.h5 %s' % h5ph)
        logger.info('tecandphase calibration...')

        with w.if_todo('cor-tecandph-c%02i' % c):
            # correct ph - ms:DATA -> ms:CORRECTED_DATA
            logger.info('Correct tec...')
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=tec000',
                        log='$nameMS_correct-c%02i.log' % c, commandType='DP3')
            logger.info('Correct ph...')
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=phase000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DP3')

    if doamp:
        if ampSolMode == 'diagonal':
            with w.if_todo('cal-amp1-c%02i' % c):
                logger.info('Gain amp calibration 1 (solint: %i)...' % solint_amp)
                # Calibration - ms:CORRECTED_DATA
                # possible to put nchan=6 if less channels are needed in the h5parm (e.g. for IDG)
                MSs_extract.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp1.h5 \
                    sol.mode=diagonal sol.solint=' + str(solint_amp) + ' sol.nchan=1 sol.smoothnessconstraint=4e6 sol.minvisratio=0.5 \
                    sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]]', \
                            log='$nameMS_solGamp1-c%02i.log' % c, commandType='DP3')

                losoto_parsets = [parset_dir + '/losoto-clip.parset', parset_dir + '/losoto-norm.parset',
                                      parset_dir + '/losoto-plot2.parset']
                lib_util.run_losoto(s, 'amp1', [ms + '/cal-amp1.h5' for ms in MSs_extract.getListStr()], losoto_parsets,
                                    plots_dir='extract/plots-%s' % c)
                os.system('mv cal-amp1.h5 %s' % h5amp1)

            with w.if_todo('cor-amp1-c%02i' % c):
                logger.info('Correct amp 1...')
                # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb=' + h5amp1 + ' cor.correction=amplitude000',
                            log='$nameMS_correct-c%02i.log' % c, commandType='DP3')

            with w.if_todo('cal-amp2-c%02i' % c):
                logger.info('Gain amp calibration 2 (solint: %i)...' % solint_amp2)
                # Calibration - ms:SMOOTHED_DATA
                MSs_extract.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp2.h5 \
                    sol.mode=diagonal sol.solint=' + str(
                    solint_amp2) + ' sol.nchan=1  sol.smoothnessconstraint=10e6 sol.minvisratio=0.5', \
                            log='$nameMS_solGamp2-c%02i.log' % c, commandType='DP3')

                losoto_parsets = [parset_dir + '/losoto-clip2.parset', parset_dir + '/losoto-norm.parset',
                                  parset_dir + '/losoto-plot3.parset']
                lib_util.run_losoto(s, 'amp2', [ms + '/cal-amp2.h5' for ms in MSs_extract.getListStr()], losoto_parsets,
                                    plots_dir='extract/plots-%s' % c)
                os.system('mv cal-amp2.h5 %s' % h5amp2)

            with w.if_todo('cor-amp2-c%02i' % c):
                logger.info('Correct amp 2...')
                # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                    cor.parmdb=' + h5amp2 + ' cor.correction=amplitude000',
                            log='$nameMS_correct-c%02i.log' % c, commandType='DP3')
        ### DONE

        elif ampSolMode == 'fulljones':
            # Smoothing - ms:DATA -> ms:SMOOTHED_DATA
            with w.if_todo('smooth-c%02i' % c):
                logger.info('BL-based smoothing...')
                MSs_extract.run('BLsmooth.py -c 1 -n 8 -r -i CORRECTED_DATA -o SMOOTHED_CORRECTED_DATA $pathMS',
                                log='$nameMS_smooth.log',
                                commandType='python', maxThreads=1)
                ### DONE

            with w.if_todo('cal_fulljones_%02i' % c):
                logger.info('Solving full-Jones...')
                MSs_extract.run(f'DP3 {parset_dir}/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_CORRECTED_DATA '
                        f'sol.h5parm=$pathMS/cal-fulljones.h5 sol.mode=fulljones sol.smoothnessconstraint=5e6 sol.nchan=1 sol.solint={solint_amp2}',
                        log=f'$nameMS_solFulljones-c{c}.log', commandType="DP3")

                lib_util.run_losoto(s, f'fulljones', [ms + '/cal-fulljones.h5' for ms in MSs_extract.getListStr()], \
                                    [parset_dir + '/losoto-fulljones.parset'], plots_dir='extract/plots-%s' % c)
                os.system('mv cal-fulljones.h5 %s' % h5fj)

            # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
            with w.if_todo('cor_fulljones_c%02i' % c):
                logger.info('Full-Jones correction...')
                MSs_extract.run(f'DP3 {parset_dir}/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                        f'cor.correction=fulljones cor.parmdb={h5fj} cor.soltab=[amplitude000,phase000]',
                        log=f'$nameMS_cor_gain-c{c}.log', commandType='DP3')

    with w.if_todo('image-c%02i' % c):
        logger.info('Imaging...')
        # if we have more than one close pointing, need to apply idg beam each iteration
        clean('c%02i' % c, MSs_extract, size=(1.1*target_reg.get_width(),1.1*target_reg.get_height()), apply_beam=do_beam, userReg=userReg, datacol='CORRECTED_DATA') # size 2 times radius  , apply_beam = c==maxniter

    # get noise, if larger than 98% of prev cycle: break
    extract_image = lib_img.Image('img/extractM-c%02i-MFS-image.fits' % c, userReg=userReg)
    rms_noise, mm_ratio = extract_image.getNoise(), extract_image.getMaxMinRatio()

    extract_image.selectCC(checkBeam=False)
    sm = lsm.load(extract_image.skymodel_cut)
    total_flux = np.sum(sm.getColValues('I', aggregate='sum'))
    logger.info('RMS noise (c:%02i): %f' % (c, rms_noise))
    logger.info('MM ratio (c:%02i): %f' % (c, mm_ratio))
    logger.info('Total flux (c:%02i): %f Jy' % (c, total_flux))

    if rms_noise < rms_noise_pre:
        best_iter = c
    else: best_iter = c - 1

    if ampcal.lower =='true':
        if (rms_noise > 0.98 * rms_noise_pre and mm_ratio < 1.01 * mm_ratio_pre) or rms_noise > 1.2 * rms_noise_pre:
            if (mm_ratio < 10 and c >= 2) or (mm_ratio < 20 and c >= 3) or (c >= 5):
                break
    elif ampcal.lower == 'false':
        pass
    else:
        if (rms_noise > 0.98 * rms_noise_pre and mm_ratio < 1.01 * mm_ratio_pre) or rms_noise > 1.2 * rms_noise_pre:
            if (mm_ratio < 10 and c >= 2) or (mm_ratio < 20 and c >= 3) or (c >= 4):
                break

    if c >= 3 and mm_ratio >= 20:
        if ampcal.lower == 'true':
            logger.info('Starting amplitude calibration in next cycle...')
            doamp = True
        elif ampcal.lower == 'false':
            logger.info('Amplitude calibration set to false. Just using phase...')
            doamp = False
        else:
            logger.info('Starting amplitude calibration in next cycle...')
            doamp = True

    rms_noise_pre = rms_noise
    mm_ratio_pre = mm_ratio

# Finally:
with w.if_todo('final_apply'):
    if best_iter != c: # If last iteration was NOT the best iteration, apply best iteration.
        logger.info('Best iteration: second to last cycle ({})'.format(best_iter))
        h5ph = 'extract/cal-ph-c%02i.h5' % best_iter
        h5amp1 = 'extract/cal-amp1-c%02i.h5' % best_iter
        h5amp2 = 'extract/cal-amp2-c%02i.h5' % best_iter
        # correct ph - ms:DATA -> ms:CORRECTED_DATA
        logger.info('Correct ph...')
        MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                     cor.parmdb=' + h5ph + ' cor.correction=phase000',
                log='$nameMS_correct-final.log', commandType='DP3')

        if phSolMode == 'tecandphase':
            logger.info('Correct tec...')
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + h5ph + ' cor.correction=tec000',
                    log='$nameMS_correct-c%02i.log' % c, commandType='DP3')

        if os.path.exists(h5amp1):
            logger.info('Correct amp 1...')
            # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp1 + ' cor.correction=amplitude000',
                    log='$nameMS_correct-final.log', commandType='DP3')

        if os.path.exists(h5amp2):
            logger.info('Correct amp 2...')
            # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                cor.parmdb=' + h5amp2 + ' cor.correction=amplitude000',
                    log='$nameMS_correct-final.log', commandType='DP3')
    else:
        logger.info('Best iteration: last cycle ({})'.format(best_iter))


with w.if_todo('imaging_final'):
    logger.info('Final imaging w/ beam correction...')
    if use_idg == 'True':
        clean('final', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=True, userReg=userReg, datacol='CORRECTED_DATA')# size 2 times radius  , apply_beam = c==maxniter
    else:
        clean('final', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), userReg=userReg, datacol='CORRECTED_DATA')  # size 2 times radius  , apply_beam = c==maxniter

with w.if_todo('imaging_highres'):
     logger.info('Producing high resolution image...')
     if use_idg == 'True':
        clean('highres', MSs_extract, res='high', size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=True, userReg=userReg, datacol='CORRECTED_DATA', update_model=False)
     else:
         clean('highres', MSs_extract, res='high', size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), userReg=userReg, datacol='CORRECTED_DATA', update_model=False)

if sourcesub == True:
    logger.info('Doing compact sources subtraction + lowres imaging')
    with w.if_todo('find_compact_sources'):

        clean('sub-highres', MSs_extract, res='ultrahigh', size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), userReg=userReg, minuv=minuv_forsub,
              datacol='CORRECTED_DATA', update_model=False)

    with w.if_todo('produce_mask'):
        logger.info('Subtracting compact sources...')
        #os.system(f'MakeMask.py --RestoredIm img/{highimagename} --Th 3')
        highres_image = lib_img.Image('img/extractM-sub-highres-MFS-image.fits')
        mask_highres = highres_image.imagename.replace('.fits', '_mask-highres.fits')
        highres_image.makeMask(threshpix=3, atrous_do=True, maskname=mask_highres, write_srl=True, write_ds9=True)

        fitsmask = 'extractM-sub-highres-MFS-image_mask-highres.fits'
        clean('compactmask', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), fits_mask=f'img/{fitsmask}',
              minuv=minuv_forsub, res='ultrahigh', datacol='CORRECTED_DATA')

    with w.if_todo('source_subtraction'):
        logger.info('Adding DIFFUSE_SUB column to datasets...')
        MSs_extract.addcol('DIFFUSE_SUB', 'DATA', usedysco=True, log='$nameMS_adddiffsubcol.log')
        MSs_extract.run('taql "update $pathMS set DIFFUSE_SUB=CORRECTED_DATA-MODEL_DATA"', log='$nameMS_hressubtract.log', commandType='general')

        logger.info('Final imaging with compact sources subtracted...')

        clean('sourcesubtracted', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), datacol='DIFFUSE_SUB',
              res='low', update_model=False)