#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# HE: this is an experimental script to split out individual directions from LBA observations with IS present
#     use cases are e.g. the preparation of in-field calibrators

# TODO: possible improvement: interpolate IS DI phase solutions

import os, glob, argparse
import numpy as np
import astropy.wcs
import warnings
from losoto.h5parm import h5parm
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd


#####################################################
def test_image_dutch(MSs, imgname, data_col='SUBTRACTED_DATA'):
    if False:
        logger.warning('skip image')
    else:
        """ Create a quick debug image..."""
        lib_util.run_wsclean(s, 'wsclean-test.log', MSs.getStrWsclean(), name=f'img/{imgname}',
                             data_column=data_col, size=3000, scale=f'4arcsec',
                             weight='briggs -0.5', niter=100000, gridder='wgridder', parallel_gridding=6,
                             no_update_model_required='', minuv_l=30, maxuvw_m=max_uvw_m_dutch, mgain=0.88, nmiter=10,
                             parallel_deconvolution=512, auto_threshold=3.0, auto_mask=5.0,
                             join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6),
                             deconvolution_channels=3, baseline_averaging='',
                             multiscale='', multiscale_scale_bias=0.7, pol='i')

def test_image_is(MSs, imgname, data_col='DATA'):
    if False:
        logger.warning('skip image')
    else:
        """ Create a quick debug image..."""
        lib_util.run_wsclean(s, 'wsclean-test.log', MSs.getStrWsclean(), name=f'img/{imgname}',
                             data_column=data_col, size=3000, scale=f'0.15arcsec',
                             weight='briggs -1.1', niter=10000, gridder='wgridder', parallel_gridding=6,
                             no_update_model_required='', minuv_l=30, mgain=0.75, nmiter=10,
                             auto_threshold=3.0, auto_mask=5.0,
                             join_channels='', fit_spectral_pol=3, multiscale_max_scales=5, channels_out=MSs.getChout(4.e6),
                             multiscale='', multiscale_scale_bias=0.7, pol='i')
#####################################################
parser = argparse.ArgumentParser(description='Split out a single direction by subtracting the rest field and correcting the stations.')
### Options for all modes
parser.add_argument('--mode', type=str, help='Either infield, ddcal or widefield.')
parser.add_argument('--infieldreg', default=None, type=str, help='Provide a region for the infield calibrator (ds9 circle or square).')
parser.add_argument('--dutchdir', type=str, default=None, help='Directory of the dutch processing.')
parser.add_argument('--ddserialcycle', type=int, default=0, help='cycle to use.')
parser.add_argument('--mss', type=str, default=None, help='Directory containing the IS MSs (after timesplit).')
### Options for splitting of infield or ddcal
parser.add_argument('--freqres', type=float, default=0.195312, help='Freq. resolution of the split-off MSs in Mhz. Default=0.195312MHz (1 subband)')
parser.add_argument('--timeres', type=int, default=16, help='Time resolution of the split-off MSs in s. Default: 16s.')
### Options for ddcal or widefield imaging
parser.add_argument('--infieldh5', type=str, default=None, help='Path towards IS di solutions that should be applied to the MSs (infield delay calibrator).')
### Options for ddcal only
parser.add_argument('--dirreg', nargs='+', default=None, help='Provide a (list) of ddcal regions to be split of (ds9 circle or square).')
args = parser.parse_args()

mode = args.mode
dirregfile_list, infieldregfile, dutchdir, mss_path = args.dirreg, args.infieldreg, args.dutchdir, args.mss
ddserialcycle = args.ddserialcycle
time_resolution, freq_resolution = args.timeres, args.freqres
infield_h5 = args.infieldh5

logger_obj = lib_log.Logger(f'pipeline-splitdir')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker(f'pipeline-splitdir.walker')
warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_splitdir','parset_dir')

# check input
# general input
if mode not in ['infield', 'ddcal', 'widefield']:
    raise ValueError('mode must be either infield or ddcal or widefield.')
if not infieldregfile:
    raise ValueError('No infield direction region file provided (--infieldreg.')
if not dutchdir:
    raise ValueError('No dutch processing directory provided.')

if mode in ['widefield','ddcal']:
    if not infield_h5: raise ValueError('No infield direction region file provided (--infieldreg.')
if mode == 'ddcal':
    if not dirregfile_list: raise ValueError('No ddcal direction region file provided (--dirreg.')


if not os.path.exists(f'img'):
    os.makedirs(f'img')

if not os.path.exists(f'splitdir'):
    os.makedirs(f'splitdir')

infield_reg = lib_util.Region_helper(infieldregfile)
infield_center = infield_reg.get_center()  # center of the infield cal region

#####################################################
# 1. Get MSs with IS from LiLF/LOFAR_timesplit output and apply dutch direction-independent solutions (FR, amp and phase) to column DATA
with w.if_todo('correct_dutch_di'):
    lib_util.check_rm('mss-hires')
    os.makedirs(f'mss-hires')
    # TODO change order of FR and di-amp for newer data
    # load original MSs - those will NOT be manipulated
    MSs_orig = lib_ms.AllMSs(glob.glob(mss_path + '/*MS'), s, check_flags=False, check_sun=False)
    logger.info('Correcting DI amplitude (Dutch stations) DATA -> DATA...')
    # Correct FR- MSs-orig/TC.MS:DATA -> MSs/TC.MS:DATA
    MSs_avg = lib_ms.AllMSs(glob.glob(dutchdir + '/mss-avg/*MS'), s, check_flags=False, check_sun=False)
    # we cut some channels from Dutch MSs in ddserial - make sure to cut the same amount of channels here
    freqres_dutch = MSs_avg.getListObj()[0].getChanband()
    freqres_is = MSs_orig.getListObj()[0].getChanband()
    msin_nchan = int(MSs_avg.getListObj()[0].getNchan()*freqres_dutch/freqres_is) # 1920 for A2255 data
    MSs_orig.run(f'DP3 {parset_dir}/DP3-cor.parset msin.nchan={msin_nchan} msin=$pathMS msout=mss-hires/$nameMS.MS msout.datacolumn=DATA cor.parmdb={dutchdir}/ddparallel/solutions/cal-amp-di.h5 \
            cor.correction=amplitudeSmooth', log='$nameMS_corFR.log', commandType='DP3')
    MSs = lib_ms.AllMSs(glob.glob('mss-hires/*MS'), s, check_flags=False, check_sun=False)
    logger.info('Correcting subfield phase (Dutch stations) DATA -> DATA...')
    # Correct MSs>DATA -> DATA
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=DATA \
            cor.parmdb={dutchdir}/ddparallel/solutions/cal-tec-sf-c1.h5 cor.correction=phase000',
            log='$nameMS_sf-correct.log', commandType='DP3')
    # Correct MSs:DATA -> DATA
    logger.info('Correcting FR (Dutch stations) DATA -> DATA...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=DATA \
                cor.parmdb={dutchdir}/ddparallel/solutions/cal-fr.h5 cor.correction=rotationmeasure000 ',
            log='$nameMS_sidelobe_corrupt.log', commandType='DP3')


MSs = lib_ms.AllMSs(glob.glob('mss-hires/*MS'), s, check_flags=False, check_sun=False)
orig_center = MSs.getListObj()[0].getPhaseCentre()
max_uvw_m_dutch = 1.05*MSs.getMaxBL(check_flags=True, dutch_only=True) # this is important, it is the maximum uvw value in meters of any dutch-dutch baseline. Everything above this value is certainly IS data

# 2. Prepare the h5parm solutions from dutch processing so they can be used for IS data (reorder dirs and add IS stations with unit solutions)
with w.if_todo('interph5'):
    logger.info('Add IS stations to Dutch dd-solutions..')
    os.system(f'cp {dutchdir}/ddserial/c0{ddserialcycle}/solutions/interp.h5 interp.h5')
    with h5parm('interp.h5') as h5:
        h5_timeres = np.diff(h5.getSolset('sol000').getSoltab('phase000').getAxisValues('time'))[0]
        h5_freqres = np.diff(h5.getSolset('sol000').getSoltab('phase000').getAxisValues('freq'))[0]
    solint = int(np.round(h5_timeres/MSs.getListObj()[0].getTimeInt()))
    nchan = int(round(h5_freqres/np.diff(MSs.getFreqs())[0]))
    # first we need to reorder the soltab dir axis to have the same order as the solset.getSou() dict, otherwise h5_merger creates a mess (best would be to fix this in h5_merger)
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
    s.add(f'h5_merger.py --h5_out interp_merged.h5 --h5_tables interp.h5 -ms "mss-hires/TC*.MS" --freq_av {nchan} --time_av {solint} --add_ms_stations --no_antenna_crash --propagate_flags')
    s.run(check=True)
    # logger.info('Test image (dd-corrected)...')
    # lib_util.run_wsclean(s, 'wsclean-c.log', MSs.getStrWsclean(), name='img/dutchddcorrmerged',
    #                      data_column='CORRECTED_DATA', size=3000, scale='4arcsec',
    #                      weight='briggs -0.3', niter=100000, gridder='wgridder', parallel_gridding=32,
    #                      no_update_model_required='', minuv_l=30, nmiter=20, mgain=0.85,
    #                      parallel_deconvolution=1024,
    #                      auto_threshold=3.0, auto_mask=5.0, join_channels='', fit_spectral_pol=3,
    #                      channels_out=6, deconvolution_channels=3,
    #                      multiscale='', multiscale_scale_bias=0.65, pol='i', beam_size=15,
    #                      apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='',
    #                      facet_regions=f'{dutchdir}/ddserial/c00/images/wideDD-c00_facets.reg', maxuvw_m=max_uvw_m_dutch,
    #                      apply_facet_solutions='interp_merged.h5 phase000,amplitude000')

################################ split our infield or ddcal ################################
if mode in ['infield', 'ddcal']:
    # 3. Predict corrupted visibilities for the full field - set to zero for all non-dutch baselines!
    with w.if_todo('predict'):
        # prepare model of central/external regions
        logger.info('Predict corrupted model of full field (wsclean)...')
        image_channels = len(
            glob.glob(f"{dutchdir}/ddserial/c0{ddserialcycle}/images/wideDD-c0{ddserialcycle}*-fpb.fits"))
        s.add(
            f'wsclean -predict -padding 1.8 -name {dutchdir}/ddserial/c0{ddserialcycle}/images/wideDD-c0{ddserialcycle} -j {s.max_cpucores} -channels-out {image_channels} \
                -facet-regions {dutchdir}/ddserial/c0{ddserialcycle}/solutions/facets-c0{ddserialcycle}.reg -maxuvw-m {max_uvw_m_dutch} -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
                -apply-facet-solutions interp_merged.h5 phase000,amplitude000 {MSs.getStrWsclean()}',
            log='wscleanPRE.log', commandType='wsclean')
        s.run(check=True)
        # Set to zero for non-dutch baselines
        logger.info('Set MODEL_DATA=0 for IS baselines...')
        MSs.run(
            "taql 'update $pathMS set MODEL_DATA=0 WHERE ANTENNA1 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME !~p/[CR]S*/]'",
            log='$nameMS_resetISmodel.log', commandType='general')
        MSs.run(
            "taql 'update $pathMS set MODEL_DATA=0 WHERE ANTENNA2 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME !~p/[CR]S*/]'",
            log='$nameMS_resetISmodel.log', commandType='general')

    # 4.Subtract the full field from all Dutch baselines
    with w.if_todo('subtract'):
        MSs.addcol('SUBTRACTED_DATA', 'DATA')
        logger.info('Subtracting external region model (SUBTRACTED_DATA = DATA - MODEL_DATA)...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"',
                log='$nameMS_taql.log', commandType='general')
        test_image_dutch(MSs, 'dutchsub')

    regfile_list = [infieldregfile] if mode == 'infield' else dirregfile_list
    # iterate over all regions, per region re-add the sources on the region, correct the sols
    for i, regfile in enumerate(regfile_list):
        name = regfile.split('.reg')[0]
        logger.info(f'Splitting out direction {name} ({i+1}/{len(regfile_list)})...')
        dir_reg = lib_util.Region_helper(regfile)
        dir_center = dir_reg.get_center()  # center of the ddcal region
        if not os.path.exists(f'mss-{name}'):
            os.makedirs(f'mss-{name}')
        d = lib_dd.Direction(name)
        center = dir_center
        d.region_file = regfile
        d.set_position(center, orig_center)
        d.set_region_facets(facets_region_file=f'{dutchdir}/ddserial/c0{ddserialcycle}/solutions/facets-c0{ddserialcycle}.reg', loc='splitdir')

        # 5. Predict back the corrupted visibilities for the infield direction -  set to zero for all non-dutch baselines!
        with w.if_todo(f'addback-{name}'):
            # prepare model of central/external regions
            logger.info('Blanking model: all but direction region...')
            for im in glob.glob(f'{dutchdir}/ddserial/c0{ddserialcycle}/images/wideDD-c0{ddserialcycle}*model*fpb.fits'):
                wideMext = 'splitdir/' + im.replace(f'wideDD-c0{ddserialcycle}',f'wideDD-c0{ddserialcycle}-{name}').split('/')[-1]
                os.system('cp %s %s' % (im, wideMext))
                lib_img.blank_image_reg(wideMext, regfile, blankval=0., inverse=True)
            # Recreate MODEL_DATA of external region for re-adding
            logger.info('Predict corrupted model of direction region (wsclean)...')
            s.add(f'wsclean -predict -padding 1.8 -name splitdir/wideDD-c0{ddserialcycle}-{name} -j {s.max_cpucores} -channels-out {len(glob.glob(f"splitdir/wideDD-c0{ddserialcycle}-{name}*fpb.fits"))} \
                    -facet-regions {d.get_region_facets()} -maxuvw-m {max_uvw_m_dutch} -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam -no-solution-directions-check \
                    -apply-facet-solutions interp_merged.h5 phase000,amplitude000 {MSs.getStrWsclean()}',
                  log='wscleanPRE.log', commandType='wsclean')
            s.run(check=True)
            # Set to zero for non-dutch baselines
            MSs.run("taql 'update $pathMS set MODEL_DATA=0 WHERE ANTENNA1 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME !~p/[CR]S*/]'", log='$nameMS_resetISmodel.log', commandType='general')
            MSs.run("taql 'update $pathMS set MODEL_DATA=0 WHERE ANTENNA2 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME !~p/[CR]S*/]'", log='$nameMS_resetISmodel.log', commandType='general')

            logger.info('Add back direction region model (CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA)...')
            MSs.addcol('CORRECTED_DATA', 'DATA')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')

        correct_col = 'CORRECTED_DATA'

        # 6. apply INFIELD direction solutions for dutch baselines
        # TODO for now we correct DUTCH stations for the INFIELD CALIBRATOR direction. Test also using the ddcal direction for the case of splitting off ddcals.
        with w.if_todo(f'correct_dutch_dd-{name}'):
            # apply init - closest DDE sol
            # TODO: this assumes phase000 and optionally, amplitude000
            h5init = h5parm('interp_merged.h5')
            solset_dde = h5init.getSolset('sol000')
            # get closest dir to target reg center
            dirs = np.array([solset_dde.getSou()[k] for k in solset_dde.getSoltab('phase000').dir])
            dir_dist = lib_util.distanceOnSphere(dirs[:, 0], dirs[:, 1], *np.deg2rad(infield_center), rad=True)
            closest = solset_dde.getSoltab('phase000').dir[np.argmin(dir_dist)]
            logger.info('Init apply: correct closest DDE solutions ({})'.format(closest))
            logger.info('Correct init ph...')
            MSs.run('DP3 ' + parset_dir + f'/DP3-cor.parset msin=$pathMS msin.datacolumn={"DATA" if mode=="widefield" else "CORRECTED_DATA"} '
                                                  f'msout.datacolumn={correct_col} cor.parmdb=interp_merged.h5 cor.correction=phase000 cor.direction=' + closest,
                            log='$nameMS_init-correct.log', commandType='DP3')
            if 'amplitude000' in solset_dde.getSoltabNames():
                logger.info('Correct init amp...')
                MSs.run('DP3 ' + parset_dir + f'/DP3-cor.parset msin=$pathMS msin.datacolumn={correct_col} msout.datacolumn={correct_col} \
                                 cor.parmdb=interp_merged.h5 cor.correction=amplitude000 cor.direction=' + closest,
                                log='$nameMS_init-correct.log', commandType='DP3')
            h5init.close()
            # test_image_dutch(MSs, 'dutchsubcorr', data_col=correct_col)
        ### DONE

        with w.if_todo(f'beamcorr-infield-{name}'):
            logger.info('Correcting beam (infield calibrator)...')
            MSs.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn={correct_col} msout.datacolumn={correct_col} \
                    corrbeam.direction=[{infield_center[0]}deg,{infield_center[1]}deg]', log='$nameMS_beam.log', commandType='DP3')
            # test_image_dutch(MSs, 'dutchsubcorrbeam1', data_col=correct_col)

        if mode == 'ddcal':
            # apply infield delay calibrator solutions to full data
            with w.if_todo(f'correct_IS_di-{name}'):
                logger.info('Correcting delay cal full solutions (IS)...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={infield_h5}  \
                        cor.correction=fulljones cor.soltab=[amplitude000,phase000] msin.datacolumn={correct_col} msout.datacolumn={correct_col} ', log='$nameMS_corIS.log', commandType='DP3')
                test_image_dutch(MSs, 'dutchsubcorrbeam1is',data_col=correct_col)

        # 6. phase shift and average the data -> to 16s,
        with w.if_todo(f'phaseshift-{name}'):
            t_avg_factor = int(round(time_resolution / MSs.getListObj()[0].getTimeInt()))
            f_avg_factor = int(round(freq_resolution * 1e6 / MSs.getListObj()[0].getChanband()))
            center = infield_center if mode == 'infield' else dir_center
            logger.info(
                f'Phase shift and avg to {time_resolution}s, {freq_resolution:.4f}MHz (x{t_avg_factor} in t; x{f_avg_factor} in f)...')
            MSs.run(
                f'DP3 {parset_dir}/DP3-shiftavg.parset msin=$pathMS msout=mss-{name}/$nameMS.MS msin.datacolumn=CORRECTED_DATA '
                f'shift.phasecenter=[{center[0]}deg,{center[1]}deg] avg.freqstep={f_avg_factor} avg.timestep={t_avg_factor}',
                log='$nameMS_shiftavg.log', commandType='DP3')
            # test_image_dutch(MSs, 'dutchsubcorrshift')

        MSs_extract = lib_ms.AllMSs(glob.glob(f'mss-{name}/*.MS'), s)

        if mode == 'ddcal':
            with w.if_todo(f'beamcorr-ddcal-{name}'):
                logger.info('Correcting beam (ddcal)...')
                MSs_extract.run('DP3 ' + parset_dir + '/DP3-beam.parset msin=$pathMS', log='$nameMS_beam.log',
                                commandType='DP3')
                # test_image_dutch(MSs_extract, 'dutchsubcorrshiftbeam', data_col='DATA')
                test_image_is(MSs_extract, 'is_shifted', data_col='DATA')
    ### DONE



################################ prepare MSs for widefield imaging ################################
if mode == 'widefield':
    # 5. apply closest direction solutions for dutch baselines
    # TODO for now we correct DUTCH stations for the INFIELD CALIBRATOR direction. Test also using the ddcal direction for the case of splitting off ddcals.
    correct_col = 'CORRECTED_DATA'
    with w.if_todo('correct_dutch_dd'):
        # apply init - closest DDE sol
        # TODO: this assumes phase000 and optionally, amplitude000
        h5init = h5parm('interp_merged.h5')
        solset_dde = h5init.getSolset('sol000')
        # get closest dir to target reg center
        dirs = np.array([solset_dde.getSou()[k] for k in solset_dde.getSoltab('phase000').dir])
        dir_dist = lib_util.distanceOnSphere(dirs[:, 0], dirs[:, 1], *np.deg2rad(infield_center), rad=True)
        closest = solset_dde.getSoltab('phase000').dir[np.argmin(dir_dist)]
        logger.info('Init apply: correct closest DDE solutions ({})'.format(closest))
        logger.info('Correct init ph...')
        MSs.run('DP3 ' + parset_dir + f'/DP3-cor.parset msin=$pathMS msin.datacolumn={"DATA" if mode=="widefield" else "CORRECTED_DATA"} '
                                      f'msout.datacolumn={correct_col} cor.parmdb=interp_merged.h5 cor.correction=phase000 cor.direction=' + closest,
                log='$nameMS_init-correct.log', commandType='DP3')
        if 'amplitude000' in solset_dde.getSoltabNames():
            logger.info('Correct init amp...')
            MSs.run('DP3 ' + parset_dir + f'/DP3-cor.parset msin=$pathMS msin.datacolumn={correct_col} msout.datacolumn={correct_col} \
                             cor.parmdb=interp_merged.h5 cor.correction=amplitude000 cor.direction=' + closest,
                    log='$nameMS_init-correct.log', commandType='DP3')
        h5init.close()
        # test_image_dutch(MSs, 'dutchsubcorr', data_col=correct_col)
    ### DONE

    with w.if_todo('beamcorr-infield'):
        logger.info('Correcting beam (infield calibrator)...')
        MSs.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn={correct_col} msout.datacolumn={correct_col} \
                corrbeam.direction=[{infield_center[0]}deg,{infield_center[1]}deg]', log='$nameMS_beam.log', commandType='DP3')
        # test_image_dutch(MSs, 'dutchsubcorrbeam1', data_col=correct_col)

    # apply infield delay calibrator solutions to full data
    with w.if_todo('correct_IS_di'):
        logger.info('Correcting delay cal full solutions (IS)...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.parmdb={infield_h5}  \
                cor.correction=fulljones cor.soltab=[amplitude000,phase000] msin.datacolumn={correct_col} msout.datacolumn={correct_col} ', log='$nameMS_corIS.log', commandType='DP3')
        test_image_dutch(MSs, 'dutchsubcorrbeam1is',data_col=correct_col)
    ### DONE

    # wide-field imaging - apply the beam again to the phase center...
    with w.if_todo('beamcorr-widefield'):
        logger.info('Correcting beam (original phase center for widefield)...')
        MSs.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                corrbeam.direction=[{orig_center[0]}deg,{orig_center[1]}deg]', log='$nameMS_beam.log',
                commandType='DP3')
        # test_image_dutch(MSs, 'dutchsubcorrbeam1isbeam2', data_col='DATA')
        test_image_is(MSs, 'is_widefield', data_col='DATA')
