#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob, argparse
import numpy as np
import lsmtool as lsm
from astropy.io import fits
from astropy.wcs import WCS
import astropy.wcs
from losoto.h5parm import h5parm
from pathlib import Path
import warnings
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
import astropy.units as u
from LiLF import lib_ms, lib_img, lib_util, lib_log

#####################################################
def clean(p, MSs, res='normal', size=[1, 1], empty=False, userReg=None, apply_beam=False, datacol=None, minuv=30, numiter=100000, fits_mask=None, update_model=True):

    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    if userReg and fits_mask:
        raise ValueError('Cannot provide both fits_mask and UserReg.')
    pixscale = MSs.resolution

    if res == 'normal':
        pixscale = float('%.1f' % (pixscale / 2.5))
    elif res == 'high':
        pixscale = float('%.1f' % (pixscale / 3.5))
    elif res == 'ultrahigh':
        pixscale = float('%.1f' % (pixscale / 3.5))
    elif res == 'low':
        pass  # no change

    imsize = [int(size[0] * 1.5 / (pixscale / 3600.)), int(size[1] * 1.5 / (pixscale / 3600.))]  # add 50%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 256: imsize[0] = 256
    if imsize[1] < 256: imsize[1] = 256

    logger.debug('Image size: ' + str(imsize) + ' - Pixel scale: ' + str(pixscale))

    if res == 'normal':
        weight = 'briggs -0.3'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.6'
        maxuv_l = None
    elif res == 'ultrahigh':
        weight = 'briggs -1'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500
    else:
        logger.error('Wrong "res": %s.' % str(res))
        sys.exit()

    if empty:
        logger.info('Cleaning empty (' + str(p) + ')...')
        imagename = 'img/empty-' + str(p)
        lib_util.run_wsclean(s, 'wscleanE-' + str(p) + '.log', MSs.getStrWsclean(), concat_mss = False, name=imagename,
                             data_column='SUBTRACTED_DATA',
                             size=imsize, scale=str(pixscale) + 'arcsec',
                             weight=weight, niter=0, no_update_model_required='', minuv_l=30, mgain=0,
                             baseline_averaging='')
    else:
        arg_dict = dict()
        if datacol:
            arg_dict['data_column'] = datacol
        # in case userReg is provided -> shallow clean, mask, merge mask with userReg, deep clean with mask
        if userReg:
            # clean 1 -- only if userReg
            logger.info('Cleaning (' + str(p) + ')...')
            imagename = 'img/extract-' + str(p)

            lib_util.run_wsclean(s, 'wscleanA-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename, concat_mss = False,
                                 size=imsize, scale=str(pixscale) + 'arcsec',
                                 weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l,
                                 mgain=0.85, parallel_deconvolution=512, auto_threshold=5, join_channels='',
                                 fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3, baseline_averaging='',**arg_dict)
            # New mask method using Makemask.py
            mask = imagename + '-MFS-image.fits'
            try:
                os.system(f'MakeMask.py --RestoredIm {mask} --Th 5 --Box 150,5')
            except:
                logger.warning('Fail to create mask for %s.' % imagename + '-MFS-image.fits')
                return
            lib_img.blank_image_reg(mask + '.mask.fits', userReg, inverse=False, blankval=1.)

        # clean 2
        # TODO: add deconvolution_channels when bug fixed
        if userReg:
            logger.info('Cleaning w/ mask (' + str(p) + ')...')
        else:
            logger.info('Cleaning (' + str(p) + ')...')

        if targetname:
            imagenameM = f'img/{targetname}-' + str(p)
        else:
            sc = SkyCoord(ra=target_ra * u.deg, dec=target_dec * u.deg, frame='icrs')
            total_ra_min = int(round(sc.ra.hour * 60.0))
            ra_h = (total_ra_min // 60) % 24
            ra_m = total_ra_min % 60
            total_dec_arcmin = int(round(abs(sc.dec.degree) * 60.0))
            dec_d = total_dec_arcmin // 60
            dec_m = total_dec_arcmin % 60
            sign = '+' if sc.dec.degree >= 0 else '-'
            imagenameM = f'img/J{ra_h:02d}{ra_m:02d}{sign}{dec_d:02d}{dec_m:02d}'

        if apply_beam:
            imsize[0] = imsize[1]
            arg_dict['use_idg'] = ''
            arg_dict['idg_mode'] = 'cpu'
            arg_dict['grid_with_beam'] = ''
            arg_dict['beam_aterm_update'] = 800
            if update_model:
                arg_dict['update_model_required'] = ''
            else:
                arg_dict['no_update_model_required'] = ''
        else:
            arg_dict['baseline_averaging'] = ''
            arg_dict['no_update_model_required'] = ''
            if update_model:
                arg_dict['do_predict'] = True
        if userReg:
            arg_dict['reuse_psf'] = imagename
            arg_dict['reuse_dirty'] = imagename
            arg_dict['fits_mask'] = mask + '.mask.fits'
        if fits_mask:
            arg_dict['fits_mask'] = fits_mask

        lib_util.run_wsclean(s, 'wscleanB-' + str(p) + '.log', MSs.getStrWsclean(), concat_mss = False, name=imagenameM,
                             size=imsize, scale=str(pixscale) + 'arcsec', weight=weight, niter=numiter,
                             minuv_l=minuv, maxuv_l=maxuv_l, mgain=0.85, multiscale='',
                             parallel_deconvolution=512, auto_threshold=0.5, auto_mask=3.0, save_source_list='',
                             join_channels='', fit_spectral_pol=3, channels_out=ch_out, **arg_dict)  # , deconvolution_channels=3)

        os.system('cat '+logger_obj.log_dir+'/wscleanB-' + str(p) + '.log | grep "background noise"')


print(r"""
  _      ____     _       ____          _            _____        _                      _    _               
 | |    | __ )   / \     |  _ \   __ _ | |_  __ _   | ____|__  __| |_  _ __  __ _   ___ | |_ (_)  ___   _ __  
 | |    |  _ \  / _ \    | | | | / _` || __|/ _` |  |  _|  \ \/ /| __|| '__|/ _` | / __|| __|| | / _ \ | '_ \ 
 | |___ | |_) |/ ___ \   | |_| || (_| || |_| (_| |  | |___  >  < | |_ | |  | (_| || (__ | |_ | || (_) || | | |
 |_____||____//_/   \_\  |____/  \__,_| \__|\__,_|  |_____|/_/\_\ \__||_|   \__,_| \___| \__||_| \___/ |_| |_|

""")


parser = argparse.ArgumentParser(description='Extraction of targets of interest from LBA observations.')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for observations.')
parser.add_argument('-radec', '--radec', dest='radec', nargs=2, type=float, default=None, help='RA/DEC where to center the extraction in deg. Use if you wish to extract only one target.')
parser.add_argument('--z', dest='redshift', type=float, default=-99, help='Redshift of the target. Not necessary unless one wants to perform compact source subtraction.')
parser.add_argument('--name', dest='name', type=str, default=None, help='Name of the target. Will be used to create the directory containing the extracted data.')
parser.add_argument('--beamcut', dest='beamcut', type=float, default=0.3, help='Beam sensitivity threshold.')
parser.add_argument('--noselfcal', dest='noselfcal', help='Do not perform selfcalibration.', action='store_true')
parser.add_argument('--extreg', dest='extreg', action='store', default=None, type=str, help='Provide an extraction region.')
parser.add_argument('--maskreg', dest='maskreg', action='store', default=None, type=str, help='Provide a user mask for cleaning.')
parser.add_argument('--ampcal', dest='ampcal', action='store', default='auto', type=str, help='Perform amplitude calibration. Can be set to True, False or auto.')
parser.add_argument('--ampsol', dest='ampsol', action='store', default='diagonal', type=str, help='How to solve for amplitudes. Can be set to diagonal or fulljones.')
parser.add_argument('--phsol', dest='phsol', action='store', default='tecandphase', type=str, help='How to solve for phases. Can be set to tecandphase or phase.')
parser.add_argument('--maxniter', dest='maxniter', type=int, default=10, help='Maximum number of selfcalibration cycles to perform.')
parser.add_argument('--subreg', dest='subreg', action='store', default=None, type=str, help='Provide an optional mask for sources that need to be removed.')
parser.add_argument('--idg', dest='idg', action='store', default='True', type=str, help='Use image domain gridding for beam correction. Set to False only in case of memory issues.')

args = parser.parse_args()
coords = args.radec
pathdir = args.path
ztarget = args.redshift
targetname = args.name
beam_cut = args.beamcut
no_selfcal = args.noselfcal
extractreg = args.extreg
userReg = args.maskreg
ampcal = args.ampcal
ampSolMode = args.ampsol
phSolMode = args.phsol
maxniter = args.maxniter
subtract_reg_file = args.subreg
use_idg = args.idg

if not pathdir:
    logger.error('Provide a path (-p) where to look for LBA observations.')
    sys.exit()

if coords is None:
    logger.error('Provide RA and DEC (--radec) in deg where to shift the phase centre.')
    sys.exit()

target_ra = coords[0]
target_dec = coords[1]

target_defined = False
if targetname:
    target_defined = True
else:
    sc = SkyCoord(ra=target_ra * u.deg, dec=target_dec * u.deg, frame='icrs')
    total_ra_min = int(round(sc.ra.hour * 60.0))
    ra_h = (total_ra_min // 60) % 24
    ra_m = total_ra_min % 60
    total_dec_arcmin = int(round(abs(sc.dec.degree) * 60.0))
    dec_d = total_dec_arcmin // 60
    dec_m = total_dec_arcmin % 60
    sign = '+' if sc.dec.degree >= 0 else '-'
    targetname = f'J{ra_h:02d}{ra_m:02d}{sign}{dec_d:02d}{dec_m:02d}'

if not os.path.exists(targetname):
    os.system(f'mkdir {targetname}')
os.chdir(str(targetname))

logger_obj = lib_log.Logger('pipeline-extract')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-extract.walker')
warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_extract','parset_dir')

if target_defined == True:
    logger.info(f'Target name: {targetname}')
else:
    logger.info(f'Target name not provided. Using coordinates to name the target: {targetname}')

logger.info(f'Beam sensitivity cut: {beam_cut}')
if no_selfcal:
    logger.info('No selfcal option is set. No selfcalibration will be performed..')
logger.info(f'Amplitude calibration set to: {ampcal}.')
logger.info(f'Type of amplitude solver: {ampSolMode}.')
logger.info(f'Type of phase solver: {phSolMode}.')
logger.info(f'Max number of selfcalibration cycles: {maxniter}')
logger.info(f'Subtraction region set: {subtract_reg_file}.')

if ampcal.lower() not in ['false', 'true', 'auto']:
    logger.error('ampcal must be true, false or auto.')
    sys.exit()

if phSolMode not in ['tecandphase', 'phase']:
    logger.error('phSolMode {} not supported. Choose tecandphase, phase.')
    sys.exit()

if ampSolMode not in ['diagonal', 'fulljones']:
    logger.error('ampSolMode {} not supported. Choose diagonal, fulljones.')
    sys.exit()

if userReg is not None:
    if not extreg:
        logger.error('To specify a mask region, an extraction region must also be provided.')
        sys.exit()
    else:
        userReg = str(mask_reg)
        logger.info("A mask for cleaning was provided.")

if subtract_reg_file is not None:
    if not extreg:
        logger.error('To specify a mask region, an extraction region must also be provided.')
        sys.exit()
    else:
        ext = 0

ext_region_extent = 0.25 #deg. This is where we start to get pointings, then we can increase the radius depending on the flux density threshold.

if ztarget != -99:
    sourcesub = True
else:
    sourcesub = False

if extractreg:
    target_reg_file = str(extractreg)
    print(target_reg_file)
    logger.info('Extraction region provided. No automatic region will be drawn...')
else:
    logger.info('Extraction region not set by the user. It will be created automatically.')
    target = lib_util.create_extregion(target_ra, target_dec, ext_region_extent)
    if os.path.exists(f'{targetname}.reg'):
        os.remove(f'{targetname}.reg')
    with open(f'{targetname}.reg', 'w') as f:
        f.write(target)
        target_reg_file = f'{targetname}.reg'

target_reg = lib_util.Region_helper(target_reg_file)
center = target_reg.get_center() # center of the extract region

list_dirs = [_d for _d in Path(str(pathdir)).iterdir() if _d.is_dir()]
tocheck = []
for dir in list_dirs:
    if (dir / 'wideDDS-c0-MFS-image-pb.fits').exists() and ((dir / 'interp.h5.gz').exists() or (dir / 'interp.h5').exists()):
        tocheck.append(dir)
        if (dir/'interp.h5.gz').exists():
            s.add(f'gunzip {dir}/interp.h5.gz')
            s.run(check=True)
close_pointings = []


for pointing in tocheck:

    chout_max = len(glob.glob(f'{pointing}/primarybeam.fits'))
    with fits.open(f'{pointing}/primarybeam.fits') as f:
        header, data = lib_img.flatten(f)
        wcs = WCS(header)
        c_pix = np.rint(wcs.wcs_world2pix([center], 0)).astype(int)[0]
        if np.all(c_pix > 0):
            try:
                beam_value = data[c_pix[1]][c_pix[0]]  # Checked -  1 and 0 are correct here.
            except IndexError:
                continue
        else:
            continue
        if beam_value > beam_cut**2: # square since Norm is sqrt(beam response)
            close_pointings.append(str(pointing).split('/')[-1])

if len(close_pointings): # only write if we find something
    with open('pointinglist.txt', 'w') as f:
        f.write('\n'.join(close_pointings))
else:
    with open('pointinglist.txt', 'r') as f:
        close_pointings = f.readlines()
        close_pointings = [line.rstrip() for line in close_pointings]


if not len(close_pointings): # raise error if none are found!
    logger.error(f'Did not find any pointing covering coordinates {target_ra}, {target_dec} with primary beam response > {beam_cut} in {pathdir}.')
    logger.error('If this is somehow unexpected, check the path (-p) and coordinates and try again.')
    logger.error('If you wish to force the extraction, you can lower the beam sensitivity threshold (default = 0.3).')
    sys.exit(1)

print('')
logger.info('The following pointings will be used:')
print('')
for name in close_pointings:
    if name != None:
        logger.info(f'Pointing {name};')
print('')

if sourcesub == True:
    sourceLLS=0.25 #Mpc.
    oneradinmpc = cosmo.angular_diameter_distance(ztarget) / (360. / (2. * np.pi))
    scalebarlengthdeg = sourceLLS / oneradinmpc.value
    minuv_forsub = 1./(scalebarlengthdeg*np.pi/180.)

with w.if_todo('cleaning'):
    logger.info('Cleaning and copying data...')
    lib_util.check_rm('extract')
    lib_util.check_rm('img')
    lib_util.check_rm('mss-extract/shiftavg')
    lib_util.check_rm('extract-images')
    os.makedirs('img')
    os.makedirs('mss-extract/shiftavg')

    for i, p in enumerate(close_pointings):
        os.makedirs('extract-files/init/'+p)
        os.system(f'cp {str(pathdir)}/{p}/wideDDS-c0-MFS-image-pb.fits extract-files/init/{p}')  # copy ddcal images
        os.system(f'cp {str(pathdir)}/{p}/wideDDS-c0-0*-model-fpb.fits extract-files/init/{p}')  # copy models
        os.system(f'cp {str(pathdir)}/{p}/interp.h5 extract-files/init/{p}')  # copy final dde sols
        os.system(f'cp {str(pathdir)}/{p}/facetsS-c0.reg extract-files/init/{p}')  # copy facet file
        lib_util.check_rm('mss-extract/'+p)
        if not os.path.exists('mss-extract/'+p):
            logger.info(f'Uncompressing .MS files of {p}...')
            os.makedirs('mss-extract/' + p)
            tgz_file = Path(pathdir) / p / f'{p}.tgz'
            if tgz_file.exists():
                os.system(f'tar -xzf {tgz_file} -C mss-extract/{p}')

if not extractreg:
    for p in close_pointings:
        image_tocheck = f'extract-files/init/{p}/wideDDS-c0-MFS-image-pb.fits'
        flux_check = lib_img.Image(image_tocheck)
        reg_flux = flux_check.calc_flux(target_reg_file)
        flux_thresh = 5 #Jy. If flux is lower than this, the extent of the extraction region gets increased.
        param=1

        maxthreshold = False
        #LET IT CREATE THE EXT REGION
        while reg_flux < flux_thresh:
            ext_region_extent += 0.084 #We add 5 arcmin every cycle
            if ext_region_extent < 0.75:
                maxthreshold = False
                target = lib_util.create_extregion(target_ra, target_dec, ext_region_extent)
                if os.path.exists(f'{targetname}.reg'):
                    os.remove(f'{targetname}.reg')
                with open(f'{targetname}.reg', 'w') as f:
                    f.write(target)
                target_reg_file = f'{targetname}.reg'
                target_reg = lib_util.Region_helper(target_reg_file)
                center = target_reg.get_center()  # center of the extract region
                reg_flux = flux_check.calc_flux(target_reg_file)
                #logger.info(f'Flux inside region of {round(ext_region_extent * 60)} arcmin radius: {round(reg_flux)}')

            else:
                maxthreshold = True
                ext_region_extent = 0.75
                target = lib_util.create_extregion(target_ra, target_dec, ext_region_extent)
                if os.path.exists(f'{targetname}.reg'):
                    os.remove(f'{targetname}.reg')
                with open(f'{targetname}.reg', 'w') as f:
                    f.write(target)
                target_reg_file = f'{targetname}.reg'
                target_reg = lib_util.Region_helper(target_reg_file)
                center = target_reg.get_center()  # center of the extract region
                break


    if maxthreshold == True:
        logger.info('Low flux (<5 Jy) detected around the target in one or more pointings. Self-calibration could be not optimal.')
        logger.info('A maximum radius of 45 arcmin was chosen for the extraction region.')
    else:
        logger.info(f'The extraction region will have a radius of {int(ext_region_extent * 60)} arcmin.')

for p in close_pointings:
    MSs = lib_ms.AllMSs( glob.glob('mss-extract/'+p+'/mss-avg/*MS'), s )
    ch_out = MSs.getChout(4e6)  # chout from dd
    fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
    phase_center = MSs.getListObj()[0].getPhaseCentre()
    # read image, h5parm, make mask
    wideDD_image = lib_img.Image(f'extract-files/init/{p}/wideDDS-c0-MFS-image-pb.fits')
    dde_h5parm = f'extract-files/init/{p}/interp.h5'
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
        #MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA', log='$nameMS_addmodelcol.log', commandType='python')

        # Predict+corrupt in MODEL_DATA of everything BUT the calibrator
        inmask = f'extract-files/init/{p}/wideDDS-c0-MFS-image-pb_mask-ddcal.fits'
        outmask = inmask + '.mask'
        lib_img.blank_image_reg(inmask, target_reg_file, outfile=outmask, inverse=False, blankval=0.)

        for im in glob.glob(f'extract-files/init/{p}/wideDDS-c0-0*model-fpb.fits'):
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
            correct_for = 'phase000,amplitude000'
        else:
            correct_for = 'phase000'

        facet_path = f'extract-files/init/{p}/facetsS-c0.reg'
        s.add(f'wsclean -predict -padding 1.8 -name extract-files/init/{p}/wideDDextS-c0 -j ' + str(s.max_cpucores) + ' -channels-out ' + str(
            ch_out) + ' -facet-regions ' + facet_path + ' -apply-facet-beam -facet-beam-update 120 -use-differential-lofar-beam \
            -apply-facet-solutions ' + dde_h5parm + ' ' + correct_for + ' \
            -reorder -parallel-reordering 4 ' + MSs.getStrWsclean(),
            log='wscleanPRE.log', commandType='wsclean')
        s.run(check=True)
        h5init.close()

    with w.if_todo('subtract_rest_'+p):

        # Remove corrupted data from DATA
        logger.info('Add columns...')
        #MSs.run('addcol2ms.py -m $pathMS -c SUBTRACTED_DATA -i CORRECTED_DATA', log='$nameMS_addcol.log', commandType='python')
        MSs.addcol('SUBTRACTED_DATA', 'DATA', log='$nameMS_addcol.log')
        logger.info('Set SUBTRACTED_DATA = DATA - MODEL_DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"',
                log='$nameMS_subtract.log', commandType='general')

    # Phase shift in the target location
    with w.if_todo('phaseshift_'+p):
        t_avg_factor = int(round(32/MSs.getListObj()[0].getTimeInt()))
        logger.info('Phase shift and avg...')
        MSs.run(f'DP3 {parset_dir}/DP3-shiftavg.parset msin=$pathMS msout=mss-extract/shiftavg/{p}_$nameMS.MS-extract msin.datacolumn=SUBTRACTED_DATA '
                f'shift.phasecenter=[{center[0]}deg,{center[1]}deg] avg.freqstep=8 avg.timestep={t_avg_factor}',
                log=p+'$nameMS_shiftavg.log', commandType='DP3')

    MSs_extract = lib_ms.AllMSs( glob.glob('mss-extract/shiftavg/'+p+'_*.MS-extract'), s )

    with w.if_todo('beamcorr_'+p):
        logger.info('Correcting beam...')
        MSs_extract.run('DP3 ' + parset_dir + '/DP3-beam.parset msin=$pathMS', log='$nameMS_beam.log', commandType='DP3')

    # apply init - closest DDE sol
    # TODO: this assumes phase000 and optionally, amplitude000
    with w.if_todo('apply_init_'+p):
        h5init = h5parm(dde_h5parm)
        solset_dde = h5init.getSolset('sol000')
        # get closest dir to target reg center
        dirs = np.array([solset_dde.getSou()[k] for k in solset_dde.getSoltab('phase000').dir])
        dir_dist = lib_util.distanceOnSphere(dirs[:,0], dirs[:,1],*np.deg2rad(center), rad=True)
        closest = solset_dde.getSoltab('phase000').dir[np.argmin(dir_dist)]
        logger.info('Init apply: correct closest DDE solutions ({})'.format(closest))
        logger.info('Correct init ph...')
        MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA '         
                'msout.datacolumn=CORRECTED_DATA cor.parmdb=' + dde_h5parm + ' cor.correction=phase000 cor.direction='+closest,
                log='$nameMS_init-correct.log', commandType='DP3')
        if 'amplitude000' in solset_dde.getSoltabNames():
            logger.info('Correct init amp...')
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                         cor.parmdb=' + dde_h5parm + ' cor.correction=amplitude000 cor.direction=' + closest,
                    log='$nameMS_init-correct.log', commandType='DP3')
        h5init.close()
    ### DONE

if no_selfcal: # finish here
    logger.info('No selfcal option is set in lilf.config. Done.')
    sys.exit(0)

MSs_extract = lib_ms.AllMSs(glob.glob('mss-extract/shiftavg/*.MS-extract'), s)

# initial imaging to get the model in the MODEL_DATA (could also be done using the Dico DDFacet model)
do_beam = len(close_pointings) > 1 # if > 1 pointing, correct beam every cycle, otherwise only at the end.
with w.if_todo('image_init'):
    logger.info('Initial imaging...')
    clean('init', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=do_beam, userReg=userReg, datacol='CORRECTED_DATA')

# Smoothing - ms:DATA -> ms:SMOOTHED_DATA
with w.if_todo('smooth'):
    logger.info('BL-based smoothing...')
    MSs_extract.run('BLsmooth.py -c 1 -n 8 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python', maxProcs=1)

# get initial noise and set iterators for timeint solutions
image = lib_img.Image(f'img/{targetname}-init-MFS-image.fits', userReg=userReg)
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

    h5ph = 'extract-files/cal-ph-c%02i.h5' % c
    solint_ph = next(iter_ph_solint)
    if doamp:
        h5amp1 = 'extract-files/cal-amp1-c%02i.h5' % c
        solint_amp = next(iter_amp_solint)
        h5amp2 = 'extract-files/cal-amp2-c%02i.h5' % c
        solint_amp2 = next(iter_amp2_solint)
        h5fj = 'extract-files/cal-fulljones-c%02i.h5' % c

    if phSolMode == 'phase':
        logger.info('Phase calibration...')
        with w.if_todo('cal-ph-c%02i' % c):
            MSs_extract.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
                     sol.mode=scalarphase sol.solint=' + str(solint_ph) + ' sol.smoothnessconstraint=5e6 sol.smoothnessreffrequency=54e6 \
                     sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]',
                     log='$nameMS_solGph-c%02i.log' % c, commandType='DP3')
            lib_util.run_losoto(s, 'ph', [ms + '/cal-ph.h5' for ms in MSs_extract.getListStr()],
                                [parset_dir + '/losoto-plot1.parset'],
                                plots_dir='extract-files/plots-%s' % c)
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
                                plots_dir='extract-files/plots-%s' % c)
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
                                    plots_dir='extract-files/plots-%s' % c)
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
                                    plots_dir='extract-files/plots-%s' % c)
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
                                commandType='python', maxProcs=1)
                ### DONE

            with w.if_todo('cal_fulljones_%02i' % c):
                logger.info('Solving full-Jones...')
                MSs_extract.run(f'DP3 {parset_dir}/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_CORRECTED_DATA '
                        f'sol.h5parm=$pathMS/cal-fulljones.h5 sol.mode=fulljones sol.smoothnessconstraint=5e6 sol.nchan=1 sol.solint={solint_amp2}',
                        log=f'$nameMS_solFulljones-c{c}.log', commandType="DP3")

                lib_util.run_losoto(s, 'fulljones', [ms + '/cal-fulljones.h5' for ms in MSs_extract.getListStr()], \
                                    [parset_dir + '/losoto-fulljones.parset'], plots_dir='extract-files/plots-%s' % c)
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
    extract_image = lib_img.Image(f'img/{targetname}-c%02i-MFS-image.fits' % c, userReg=userReg)
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
        h5ph = 'extract-files/cal-ph-c%02i.h5' % best_iter
        h5amp1 = 'extract-files/cal-amp1-c%02i.h5' % best_iter
        h5amp2 = 'extract-files/cal-amp2-c%02i.h5' % best_iter
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
        highres_image = lib_img.Image(f'img/{targetname}-sub-highres-MFS-image.fits')
        mask_highres = highres_image.imagename.replace('.fits', '_mask-highres.fits')
        highres_image.makeMask(threshpix=3, atrous_do=True, maskname=mask_highres, write_srl=True, write_ds9=True)

        fitsmask = f'{targetname}-sub-highres-MFS-image_mask-highres.fits'
        clean('compactmask', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), fits_mask=f'img/{fitsmask}',
              minuv=minuv_forsub, res='ultrahigh', datacol='CORRECTED_DATA')

    with w.if_todo('source_subtraction'):
        logger.info('Adding DIFFUSE_SUB column to datasets...')
        MSs_extract.addcol('DIFFUSE_SUB', 'DATA', usedysco=True, log='$nameMS_adddiffsubcol.log')
        MSs_extract.run('taql "update $pathMS set DIFFUSE_SUB=CORRECTED_DATA-MODEL_DATA"', log='$nameMS_hressubtract.log', commandType='general')

        logger.info('Final imaging with compact sources subtracted...')

        clean('sourcesubtracted', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), datacol='DIFFUSE_SUB',
              res='low', update_model=False)

if not os.path.exists('extract-images'):
    os.makedirs('extract-images')
os.system(f'mv img/{targetname}-final-MFS-image*.fits extract-images/')
os.system(f'mv img/{targetname}-highres-MFS-image*.fits extract-images/')
if os.path.exists(f'img/{targetname}-sourcesubtracted-MFS-image.fits'):
    os.system(f'mv img/{targetname}-sourcesubtracted-MFS-image.fits extract-images/')
if os.path.exists(f'pointinglist.txt'):
    os.system('rm pointinglist.txt')

w.alldone()

logfile = sorted(glob.glob('pipeline-extract_*.logger'))[-1]
with open(logfile, 'r') as f:
    last_line = f.readlines()[-1]
    if not "Done" in last_line:
        logger.error(f'Something went wrong in the extraction of {targetname} - check the logfile {targetname}/{logfile}.')
        raise lib_util.Exit
    else:
        logger.info(f'Target {targetname} has been extracted.')
        os.chdir('../')

        logger.info(f'Extracted datasets are in the mss-extract directory of {targetname}.')

        if sourcesub == True:
            logger.info(f'Nominal, high and source-subtracted low resolution images are in the extract-images/ directory of {targetname}.')
            logger.info('A new column "DIFFUSE_SUB" has been added to each extracted .MS file.')
            logger.info('It contains the compact-source-subtracted visibilities.')
        else:
            logger.info(f'Nominal and high resolution images are in the extract-images/ directory of {targetname}.')


