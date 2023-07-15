#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline for extraction of target region after LOFAR_dd.
# Provide a extractRegion in lilf.config. This pipeline will subtract
# sources outside of the region and perform subsequent self-calibration.
# Multiple pointings can be used for extraction.
# A userReg may be specified as clean mask.
# phSolMode can be used to solve either using phases or phaseandtec.

import sys, os, glob, re, argparse
import numpy as np
import lsmtool as lsm
from astropy.io import fits
from astropy.wcs import WCS
import astropy.wcs
from losoto.h5parm import h5parm
from pathlib import Path
import warnings
import pyrap.tables as pt
from astropy.cosmology import FlatLambdaCDM
from LiLF import lib_ms, lib_img, lib_util, lib_log


warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def get_ddf_parms_from_header(img):
    """
    Parse the HISTORY header of a DDFacet image and return a dict containing the options used to create the image.
    Will replace '-' by '_'.
    Parameters
    ----------

    img: str, filename of image

    Returns
    -------
    params_dict: dict,
    """
    params_dict = dict()
    hdr = fits.open(img)[0].header['HISTORY']
    for line in hdr:
        if 'MUFFIN' in line: continue
        if line.count('=') == 1 and line.count('-') > 0:
            _key, _value = line.replace(' ', '').split('=')
            _key = _key.replace('-', '_')
            params_dict[_key] = _value
        else:
            continue
    return params_dict


#######################################################
logger_obj = lib_log.Logger('pipeline-extract')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-extract.walker')

################################
##These two functions are to avoid excess printing from pyrap.tables.
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__
################################

def clean(p, MSs, res='normal', size=[1, 1], empty=False, userReg=None, apply_beam=False, datacol=None, minuv=30, numiter=100000, fitsmask=None):

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
        lib_util.run_wsclean(s, 'wscleanE-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
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

            lib_util.run_wsclean(s, 'wscleanA-' + str(p) + '.log', MSs.getStrWsclean(), name=imagename,
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
        imagenameM = 'img/extractM-' + str(p)
        if apply_beam:
            imsize[0] = imsize[1]
            arg_dict['use_idg'] = ''
            arg_dict['idg_mode'] = 'cpu'
            arg_dict['grid_with_beam'] = ''
            arg_dict['beam_aterm_update'] = 800
        else:
            arg_dict['baseline_averaging'] = ''
        if userReg:
            arg_dict['reuse_psf'] = imagename
            arg_dict['reuse_dirty'] = imagename
            arg_dict['fits_mask'] = mask + '.mask.fits'
        if fitsmask:
            arg_dict['fits_mask'] = fitsmask

        lib_util.run_wsclean(s, 'wscleanB-' + str(p) + '.log', MSs.getStrWsclean(), name=imagenameM,
                             do_predict=True,
                             size=imsize, scale=str(pixscale) + 'arcsec', weight=weight, niter=numiter,
                             no_update_model_required='', minuv_l=minuv, maxuv_l=maxuv_l, mgain=0.85, multiscale='',
                             parallel_deconvolution=512, auto_threshold=0.5, auto_mask=3.0, save_source_list='',
                             join_channels='', fit_spectral_pol=3, channels_out=ch_out, **arg_dict)  # , deconvolution_channels=3)

        os.system('cat '+logger_obj.log_dir+'/wscleanB-' + str(p) + '.log | grep "background noise"')


parser = argparse.ArgumentParser(description='Extraction of targets of interest from LBA survey observations.')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for observations.')

args = parser.parse_args()
pathdir = args.path

parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_extract'])))
parset_dir = parset.get('LOFAR_extract','parset_dir')
maxniter = parset.getint('LOFAR_extract','max_niter')
subtract_reg_file = parset.get('LOFAR_extract','subtract_region')  # default None - use only if you want to subtract individual sources which are in extractReg
phSolMode = parset.get('LOFAR_extract','ph_sol_mode')  # default: tecandphase
ampSolMode = parset.get('LOFAR_extract', 'amp_sol_mode') # default: diagonal
beam_cut = parset.getfloat('LOFAR_extract','beam_cut')  # default: 0.3
no_selfcal = parset.getboolean('LOFAR_extract','no_selfcal')  # Only extract, no selfcal?
userReg = parset.get('model','userReg')
ampcal = parset.get('LOFAR_extract','ampcal')

if ampcal.lower() not in ['false', 'true', 'auto']:
    logger.error('ampcal must be true, false or auto.')
    sys.exit()

if phSolMode not in ['tecandphase', 'phase']:
    logger.error('phSolMode {} not supported. Choose tecandphase, phase.')
    sys.exit()

if ampSolMode not in ['diagonal', 'fulljones']:
    logger.error('ampSolMode {} not supported. Choose diagonal, fulljones.')
    sys.exit()


ext_region_extent = 0.25 #deg. This is where we start to get pointings, then we can increase the radius depending on the flux density threshold.
data_temp = np.loadtxt('redshift_temp.txt', delimiter=' ', usecols=[0,1,2])
clname_temp = np.loadtxt('redshift_temp.txt', delimiter=' ', usecols=[3], dtype=np.str)

try:
    extreg_temp = np.loadtxt('redshift_temp.txt', delimiter=' ', usecols=[4], dtype=np.str)
    if str(extreg_temp) == 'None':
        extreg = 0
        maskreg = 0
    else:
        extreg=1
        try:
            mask_reg = np.loadtxt('redshift_temp.txt', delimiter=' ', usecols=[5], dtype=np.str)
            if str(mask_reg) == 'None':
                maskreg = 0
            else:
                maskreg = 1
        except:
            maskreg = 0
except:
    extreg = 0
    maskreg = 0

z = data_temp[0]
if z == -99: #Avoid source subtraction if no redshift info
    sourcesub = 1
    logger.info('Redshift information not found. Source subtraction will be skipped...')
else:
    sourcesub = 0

ra = data_temp[1]
dec = data_temp[2]
cluster = clname_temp

if maskreg == 1:
    if str(mask_reg).endswith('.reg'):
        userReg = str(mask_reg)
        fitsmask = None
    else:
        userReg = None
        fitsmask = str(mask_reg)
        openmask = fits.open(fitsmask)
        getheader = openmask[0].data
        wcs = WCS(openmask[0].header)
        naxis1 = openmask[0].header['NAXIS1']
        naxis2 = openmask[0].header['NAXIS2']
        pixtodeg1 =  openmask[0].header['CDELT1']
        pixtodeg2 = openmask[0].header['CDELT2']
        mask_width = round(naxis1*abs(pixtodeg1),2)
        mask_height = round(naxis2 * abs(pixtodeg2),2)
else:
    userReg = None
    fitsmask = None

if extreg==1:
    target_reg_file = str(extreg_temp)
    logger.info('Extraction region provided. No automatic region will be drawn...')
else:
    logger.info('Extraction region not set by the user. It will be created automatically...')
    target = lib_util.create_extregion(ra, dec, ext_region_extent)
    if os.path.exists('target.reg'):
        os.remove('target.reg')
    with open('target.reg', 'w') as f:
        f.write(target)
        target_reg_file = parset.get('LOFAR_extract','extractRegion')  # default 'target.reg'

target_reg = lib_util.Region_helper(target_reg_file)
center = target_reg.get_center() # center of the extract region

list_dirs = [_d for _d in Path(str(pathdir)).iterdir() if _d.is_dir()]
tocheck = []
for dir in list_dirs:
    if dir/'ddcal' in dir.iterdir() and dir/'mss-avg' in dir.iterdir():
        tocheck.append(dir)
close_pointings = []

if not os.path.exists('pointinglist.txt'):
    for pointing in tocheck:
        with fits.open(pointing/'ddcal/c01/images/wideDD-c01.MeanSmoothNorm.fits') as f:
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
    logger.error(f'Did not find any pointing covering coordinates {ra}, {dec} with primary beam response > {beam_cut} in {pathdir}.')
    logger.error(f'If this is somehow unexpected, check the path (-p) and coordinates and try again.')
    logger.error(f'If you wish to force the extraction, you can lower the beam sensitivity threshold (default = 0.3) in lilf.config.')
    sys.exit(1)

print('')
logger.info('The following pointings will be used:')
print('')
for name in close_pointings:
    if name != None:
        logger.info(f'Pointing {name};')
print('')

#Compute minuv for source subtraction. Assume to subtract everything below 400 kpc?
if sourcesub == 0:
    sourceLLS=0.25 #Mpc. TODO Should this be tuned for each cluster? Not sure how..
    oneradinmpc = cosmo.angular_diameter_distance(z) / (360. / (2. * np.pi))
    scalebarlengthdeg = sourceLLS / oneradinmpc.value
    minuv_forsub = 1./(scalebarlengthdeg*np.pi/180.)

with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('extract')
    lib_util.check_rm('img')
    lib_util.check_rm('mss-extract/shiftavg')
    os.makedirs('img')

    os.makedirs('mss-extract/shiftavg')
    for i, p in enumerate(close_pointings):
            os.makedirs('extract/init/'+p)
            os.system(f'cp {str(pathdir)}/{p}/ddcal/c01/images/wideDD-c01.app.restored.fits extract/init/{p}')  # copy ddcal images
            os.system(f'cp {str(pathdir)}/{p}/ddcal/c01/images/wideDD-c01.DicoModel extract/init/{p}')  # copy dico models
            os.system(f'cp {str(pathdir)}/{p}/ddcal/c01/solutions/interp.h5 extract/init/{p}')  # copy final dde sols
            lib_util.check_rm('mss-extract/'+p)
            if not os.path.exists('mss-extract/'+p):
                logger.info('Copying MS of '+p+'...')
                os.makedirs('mss-extract/' + p)
                os.system(f'cp -r {str(pathdir)}/{p}/mss-avg/* mss-extract/{p}')

if extreg != 1:
    for p in close_pointings:
        image_tocheck = 'extract/init/'+p+'/wideDD-c01.app.restored.fits'
        flux_check = lib_img.Image(image_tocheck)
        reg_flux = flux_check.calc_flux(image_tocheck, target_reg_file)
        flux_thresh = 5 #Jy. If flux is lower than this, the extent of the extraction region gets increased.
        param=1

        #LET IT CREATE THE EXT REGION
        while reg_flux < flux_thresh:
            #logger.info('Flux too low, increasing extraction region radius...')
            ext_region_extent += 0.084 #We add 5 arcmin every cycle
            if ext_region_extent < 0.75:
                param = 1
                target = lib_util.create_extregion(ra, dec, ext_region_extent)
                if os.path.exists('target.reg'):
                    os.remove('target.reg')
                with open('target.reg', 'w') as f:
                    f.write(target)
                target_reg_file = parset.get('LOFAR_extract', 'extractRegion')
                target_reg = lib_util.Region_helper(target_reg_file)
                center = target_reg.get_center()  # center of the extract region
                reg_flux = flux_check.calc_flux(image_tocheck, target_reg_file)
                #logger.info(f'Flux inside region of {round(ext_region_extent * 60)} arcmin radius: {round(reg_flux)}')

            else:
                param = 0
                ext_region_extent = 0.75
                target = lib_util.create_extregion(ra, dec, ext_region_extent)
                if os.path.exists('target.reg'):
                    os.remove('target.reg')
                with open('target.reg', 'w') as f:
                    f.write(target)
                target_reg_file = parset.get('LOFAR_extract', 'extractRegion')
                target_reg = lib_util.Region_helper(target_reg_file)
                center = target_reg.get_center()  # center of the extract region
                break

    if param == 0:
        logger.info('Low flux (<5 Jy) detected around the target in one or more pointings. Selfcalibration could be not optimal.')
        logger.info('A maximum radius of 45 arcmin was chosen for the extraction region.')
    else:
        logger.info(f'The extraction region will have a radius of {int(ext_region_extent * 60)} arcmin.')

for p in close_pointings:
    MSs = lib_ms.AllMSs( glob.glob('mss-extract/'+p+'/*MS'), s )
    ch_out = MSs.getChout(4e6)  # chout from dd
    fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
    phase_center = MSs.getListObj()[0].getPhaseCentre()
    # read image, h5parm, make mask
    wideDD_image = lib_img.Image('extract/init/'+p+'/wideDD-c01.app.restored.fits')
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


    with w.if_todo('predict_rest_'+p):

        # Add mock MODEL column to avoid DDFacet overflow
        # TODO: use MSs.addcol() to add MODEL_DATA w/o dysco?
        # MSs.addcol('MODEL_DATA', 'DATA', usedysco=False, log='$nameMS_addmodelcol.log')
        MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA', log='$nameMS_addmodelcol.log', commandType='python')

        # DDF predict+corrupt in MODEL_DATA of everything BUT the calibrator
        indico = wideDD_image.root + '.DicoModel'
        outdico = indico + '-' + target_reg_file.split('.')[0] # use prefix of target reg
        inmask = sorted(glob.glob(wideDD_image.root + '*_mask-ddcal.fits'))[-1]
        outmask = outdico + '.mask'
        lib_img.blank_image_reg(inmask, target_reg_file, outfile=outmask, inverse=False, blankval=0.)
        # if we have subtract reg, unmask that part again to predict+subtract it.
        if subtract_reg_file != '':
            os.system(f'cp ../{subtract_reg_file} .')
            logger.info(f"Re-adding sources in subtract-region {subtract_reg_file} to subtraction model.")
            lib_img.blank_image_reg(outmask, subtract_reg_file, inverse=False, blankval=1.)
        s.add('MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s' % (outmask, indico, outdico),
              log='MaskDicoModel.log', commandType='DDFacet', processors='max')
        s.run(check=True)

        # get DDF parameters used to create the image/model
        ddf_parms = get_ddf_parms_from_header(wideDD_image.imagename)
        h5init = h5parm(dde_h5parm)
        solset_dde = h5init.getSolset('sol000')
        # change for PREDICT
        ddf_parms['Data_MS'] = MSs.getStrDDF()
        ddf_parms['Data_ColName'] = 'CORRECTED_DATA'
        ddf_parms['Predict_ColName'] = 'MODEL_DATA'
        ddf_parms['Output_Mode'] = 'Predict'
        ddf_parms['Predict_InitDicoModel'] = outdico
        ddf_parms['Beam_Smooth'] = 1
        ddf_parms['Cache_Reset'] = 1

        if 'amplitude000' in solset_dde.getSoltabNames():
            ddf_parms['DDESolutions_DDSols'] = dde_h5parm + ':sol000/phase000+amplitude000'
        else:
            ddf_parms['DDESolutions_DDSols'] = dde_h5parm + ':sol000/phase000'
        if 'Misc_ParsetVersion' in ddf_parms.keys(): del ddf_parms['Misc_ParsetVersion']
        if 'Mask_External' in ddf_parms.keys(): del ddf_parms['Mask_External']

        logger.info('Predict corrupted rest-of-the-sky for '+p+'...')
        lib_util.run_DDF(s, 'ddfacet-pre.log', **ddf_parms)


    with w.if_todo('subtract_rest_'+p):

        # Remove corrupted data from CORRECTED_DATA
        logger.info('Add columns...')
        #MSs.run('addcol2ms.py -m $pathMS -c SUBTRACTED_DATA -i CORRECTED_DATA', log='$nameMS_addcol.log', commandType='python')
        MSs.addcol('SUBTRACTED_DATA', 'CORRECTED_DATA', log='$nameMS_addcol.log')
        logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
                log='$nameMS_subtract.log', commandType='general')

    # Phase shift in the target location
    with w.if_todo('phaseshift_'+p):
        t_avg_factor = int(round(32/MSs.getListObj()[0].getTimeInt()))
        logger.info('Phase shift and avg...')
        MSs.run(f'DP3 {parset_dir}/DP3-shiftavg.parset msin=$pathMS msout=mss-extract/shiftavg/{p}_$nameMS.MS-extract msin.datacolumn=SUBTRACTED_DATA '
                f'shift.phasecenter=[{center[0]}deg,{center[1]}deg\] avg.freqstep=8 avg.timestep={t_avg_factor}',
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
    clean('init', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=do_beam,
          userReg=userReg, datacol='CORRECTED_DATA')

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
                        f'cor.correction=fulljones cor.parmdb={h5fj} cor.soltab=\[amplitude000,phase000\]',
                        log=f'$nameMS_cor_gain-c{c:02}.log', commandType='DP3')

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
        logger.info('Best ieration: last cycle ({})'.format(best_iter))


with w.if_todo('imaging_final'):
    logger.info('Final imaging w/ beam correction...')
    clean('final', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=True, userReg=userReg, datacol='CORRECTED_DATA')# size 2 times radius  , apply_beam = c==maxniter

with w.if_todo('imaging_highres'):
     logger.info('Producing high resolution image...')
     clean('highres', MSs_extract, res='high', size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=True, userReg=userReg, datacol='CORRECTED_DATA')

if sourcesub == 0:
    logger.info('Do compact source subtraction + lowres imaging')
    with w.if_todo('find_compact_sources'):
        clean('sub-highres', MSs_extract, res='ultrahigh', size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()),  userReg=userReg, minuv=minuv_forsub, datacol='CORRECTED_DATA')

    with w.if_todo('produce_mask'):
        makemask='MakeMask.py'
        logger.info('Subtracting compact sources...')
        highimagename  = 'extractM-sub-highres-MFS-image.fits'
        os.system(f'MakeMask.py --RestoredIm img/{highimagename} --Th 3')
        fitsmask = highimagename + '.mask.fits'
        clean('compactmask', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), fitsmask='img/'+fitsmask,
              minuv=minuv_forsub, res='ultrahigh', datacol='CORRECTED_DATA')

    with w.if_todo('source_subtraction'):
        logger.info('Adding DIFFUSE_SUB column to datasets...')
        MSs_extract.addcol('DIFFUSE_SUB', 'DATA', usedysco=True, log='$nameMS_adddiffsubcol.log')
        MSs_extract.run('taql "update $pathMS set DIFFUSE_SUB=CORRECTED_DATA-MODEL_DATA"', log='$nameMS_hressubtract.log', commandType='general')

        logger.info('Final imaging with compact sources subtracted...')
        clean('sourcesubtracted', MSs_extract, size=(1.1 * target_reg.get_width(), 1.1 * target_reg.get_height()), apply_beam=True, datacol='DIFFUSE_SUB', res='low')

os.system('rm redshift_temp.txt')
logger.info('Done.')

