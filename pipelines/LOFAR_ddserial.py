#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, pickle, collections, fileinput
import numpy as np
from astropy.table import Table as astrotab
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import regions
import lsmtool
from losoto.h5parm import h5parm

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, lib_h5
logger_obj = lib_log.Logger('pipeline-ddserial')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-ddserial.walker')

# parse parset
parset = lib_util.getParset()
logger.info('Parset: '+str(dict(parset['LOFAR_ddserial'])))
parset_dir = parset.get('LOFAR_ddserial','parset_dir')
maxIter = parset.getint('LOFAR_ddserial','maxIter')
min_cal_flux60 = parset.getfloat('LOFAR_ddserial','minCalFlux60') # default: 0.8
solve_amp = parset.getboolean('LOFAR_ddserial','solve_amp')
manual_dd_cal = parset.get('LOFAR_ddserial','manual_dd_cal') # ds9 circle region file containing a manual dd-calibrator
develop = parset.getboolean('LOFAR_ddserial', 'develop') # for development, make more output/images
userReg = parset.get('model','userReg')

def clean(p, MSs, res='normal', size=[1,1], empty=False, imagereg='', masksigma=6.5):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    masksigma = sigma for masking
    """
    # set pixscale and imsize
    localpixscale = MSs.resolution

    if res == 'normal':
        localpixscale = float('%.1f'%(localpixscale/2.5))
    elif res == 'high':
        localpixscale = float('%.1f'%(localpixscale/3.5))
    elif res == 'low':
        pass # no change

    imsize = [int(size[0]*1.5/(localpixscale/3600.)), int(size[1]*1.5/(localpixscale/3600.))] # add 50%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 256: imsize[0] = 256
    if imsize[1] < 256: imsize[1] = 256

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(localpixscale))

    if res == 'normal':
        weight = 'briggs -0.5'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.7'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500
    else:
        logger.error('Wrong "res": %s.' % str(res))
        sys.exit()

    if empty:

        logger.info('Cleaning empty ('+str(p)+')...')
        imagename = 'img/empty-'+str(p)
        lib_util.run_wsclean(s, 'wscleanE-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, data_column='SUBTRACTED_DATA',
                size=imsize, scale=str(localpixscale)+'arcsec',
                weight=weight, niter=0, no_update_model_required='', minuv_l=30, mgain=0,
                baseline_averaging='')
    else:

        # clean 1
        logger.info('Cleaning ('+str(p)+')...')
        imagename = 'img/ddserial-'+str(p)
        lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename,
                size=imsize, scale=str(localpixscale)+'arcsec',
                weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85,
                baseline_averaging='', parallel_deconvolution=512, auto_threshold=5,
                join_channels='', fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3)
    
        # make mask
        if imagereg != '':
            s.add('breizorro.py -t %f -r %s -b 50 -o %s --merge %s' % (masksigma, imagename+'-MFS-image.fits', imagename+'-mask.fits', imagereg), 
                    log='makemask-'+str(p)+'.log', commandType='python')
            s.run(check=True)        
        else:
            s.add('breizorro.py -t %f -r %s -b 50 -o %s' % (masksigma, imagename+'-MFS-image.fits', imagename+'-mask.fits'), 
                    log='makemask-'+str(p)+'.log', commandType='python')
            s.run(check=True)        

        # clean 2
        # TODO: add deconvolution_channels when bug fixed
        # No beam correction as it is already applied when splitting the mss-dir
        logger.info('Cleaning w/ mask ('+str(p)+')...')
        imagenameM = 'img/ddserialM-'+str(p)
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagenameM, do_predict=True,
                size=imsize, save_source_list='', scale=str(localpixscale)+'arcsec', reuse_psf=imagename, reuse_dirty=imagename,
                weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85,
                multiscale='', multiscale_scale_bias=0.7, multiscale_scales='0,10,20,40,80', 
                baseline_averaging='',  auto_threshold=0.75, auto_mask=2.5, fits_mask=imagename+'-mask.fits',
                join_channels='', fit_spectral_pol=3, channels_out=ch_out)  #, deconvolution_channels=3) #local_rms

        os.system('cat '+logger_obj.log_dir+'/wscleanB-'+str(p)+'.log | grep "background noise"')

#############################################################
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('ddserial')
    os.makedirs('ddserial/init')
    os.system('cp ddparallel/skymodel/wideM-*model*.fits ddserial/init/')
    os.system('cp '+sorted(glob.glob("ddparallel/images/wideM-*image.fits"))[-1]+' ddserial/init/')
    if manual_dd_cal != '':
        os.system('cp '+sorted(glob.glob("ddparallel/images/wideM-*residual.fits"))[-1]+' ddserial/init/')
    lib_util.check_rm('img')
    os.makedirs('img')
    lib_util.check_rm('mss-avg')
### DONE

# use unaveraged MSs to be sure to get the same pixscale and imgsizepix of ddparallel
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s, check_consistency=True)
pixscale = MSs.getListObj()[0].getPixelScale() 
imgsizepix = int(1.85*max(MSs.getListObj()[0].getFWHM(freq='max', elliptical=True)) * 3600 / pixscale) # roughly to smallest null
if imgsizepix > 10000: imgsizepix = 10000 # keep SPARSE doable
if imgsizepix % 2 != 0: imgsizepix += 1  # prevent odd img sizes

# goes down to 8 seconds and multiple of 48 chans
# data should be at multiple of 48 channels already from timesplit, check is redundant and can be removed late 2025.
if not os.path.exists('mss-avg'):
    timeint = MSs.getListObj()[0].getTimeInt()
    avgtimeint = int(round(8/timeint))  # to 8 seconds
    nchan_init = MSs.getListObj()[0].getNchan()
    # chan: avg (x8) sol (x6) - we need a multiple of 8x6=48, the largest that is <nchan
    # survey after avg (x8): 60, final number of sol 10
    # pointed after avg (x8): 120, final number of sol 20
    nchan = nchan_init - nchan_init % 48
    os.makedirs('mss-avg')
    logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (timeint,timeint*avgtimeint,nchan_init,nchan))
    MSs.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msout=mss-avg/$nameMS.MS msin.datacolumn=CORRECTED_DATA msin.nchan='+str(nchan)+' \
            avg.timestep='+str(avgtimeint)+' avg.freqstep=1', log='$nameMS_initavg.log', commandType='DP3')

MSs = lib_ms.AllMSs(glob.glob('mss-avg/TC*[0-9].MS'), s, check_flags=True)
fwhm = max(MSs.getListObj()[0].getFWHM(freq='mid', elliptical=True))
workingReg = 'ddserial/workingRegion.reg' # sources outside of this region will be ignored (and not peeled)
MSs.getListObj()[0].makeBeamReg(workingReg, freq='min', to_pbval=0)
peelReg = 'ddserial/peelingRegion.reg' # sources outside of this region will be peeled
MSs.getListObj()[0].makeBeamReg(peelReg, freq='max', to_pbval=0.12) # this is slighly smaller than the null
freq_min = np.min(MSs.getFreqs())
freq_mid = np.mean(MSs.getFreqs())
phase_center = MSs.getListObj()[0].getPhaseCentre()
timeint = MSs.getListObj()[0].getTimeInt()
ch_out = MSs.getChout(4e6)  # for full band (48e6 MHz) is 12

facetregname = 'ddparallel/solutions/facets-c1.reg'
ddparallel_h5parm = 'ddparallel/solutions/cal-tec-merged-c1.h5'
interp_h5parm = ddparallel_h5parm # this variable points to the most recent dd-solutions and will be updated before each imaging step
correct_for = 'phase000'  # if needed add amplitudes000 before imaging

with w.if_todo('add_columns'):
    logger.info('Add columns...')
    # using mix of ms.addcol and addcol2ms because ms.addcol does not work with non-data columns
    MSs.addcol('CORRECTED_DATA', 'DATA', log='$nameMS_addcol.log')
    MSs.addcol('SUBTRACTED_DATA', 'DATA', log='$nameMS_addcol.log')
    MSs.run('addcol2ms.py -m $pathMS -c FLAG_BKP -i FLAG', log='$nameMS_addcol.log', commandType='python')

# print fractional flags
for MS in MSs.getListObj():
    logger.info(f'{MS.nameMS}: Fractional flags: {MS.fractionalFlag()*100:.1f}%.')

##############################################################
full_image = lib_img.Image('ddserial/init/wideM-1-MFS-image.fits', userReg=userReg)

for cmaj in range(maxIter):
    logger.info('Starting major cycle: %i' % cmaj)

    # cycle specific variables
    picklefile = 'ddserial/directions-c%02i.pickle' % cmaj
    mask_ddcal = 'ddserial/c%02i/skymodels/mask-ddcal-c%02i.fits' % (cmaj, cmaj)

    if not os.path.exists('ddserial/c%02i' % cmaj): os.makedirs('ddserial/c%02i' % cmaj)
    for subdir in ['plots','images','solutions','skymodels']:
        if not os.path.exists('ddserial/c%02i/%s' % (cmaj, subdir)): os.makedirs('ddserial/c%02i/%s' % (cmaj, subdir))

    if not os.path.exists(picklefile):
        directions = []

        if not os.path.exists(mask_ddcal.replace('fits', 'cat.fits')): # re-use if exists
            # making skymodel from image, used to find calibrators (mask-ddcal-c0x.cat.fits)
            full_image.makeMask(threshpix=4, atrous_do=False, maskname=mask_ddcal, write_srl=True, write_ds9=True)
        global_rms = full_image.getNoise()
        
        # locating DD-calibrators
        cal = astrotab.read(mask_ddcal.replace('fits','cat.fits'), format='fits')
        cal.remove_rows((cal['Total_flux'] < 0.01) | (cal['Isl_Total_flux'] < 0.1)) # remove very faint to speedup
        cal['Cluster_id'] = 'None           '
        cal['Flux_ratio'] = cal['Total_flux']/cal['Peak_flux']
        # grouping nearby sources
        grouper = lib_dd.Grouper(list(zip(cal['RA'],cal['DEC'])), cal['Total_flux'],
                                 look_distance=0.1, kernel_size=0.07, grouping_distance=0.03)
        grouper.run()
        grouper.grouping()
        grouper.plot()
        # assert sources in same island are in same group
        populated_isl = [isl for isl, n in collections.Counter(cal['Isl_id']).items() if n > 1]  # isls with more than 1 source
        ids_to_merge = [np.flatnonzero(cal['Isl_id'] == this_isl) for this_isl in populated_isl]  # list of lists of source ids that belong to populated cluster
        [grouper.merge_ids(ids) for ids in ids_to_merge]  # merge ids so to rejoin islands
        clusters = grouper.clusters
        os.system('mv grouping*png ddserial/c%02i/plots/' % cmaj)
        img_beam = full_image.getBeam()

        # reorder clusters based on flux
        fluxes = []
        for cluster_num, cluster_idxs in enumerate(clusters):
            fluxes.append(np.sum(cal['Total_flux'][cluster_idxs]))
        clusters = [x for _,x in sorted(zip(fluxes,clusters))][::-1]

        logger.info('Finding direction calibrators...')
        for cluster_num, cluster_idxs in enumerate(clusters):
            name = 'ddcal%04i' % cluster_num
            cal['Cluster_id'][cluster_idxs] = name  # just for debug
            fluxes = np.sum(cal['Total_flux'][cluster_idxs])
            spidx_coeffs = -0.8
            to_60_MHz = (60e6/freq_mid)**spidx_coeffs
            localrms = np.max(cal['Isl_rms'][cluster_idxs])


            # remove clusters with diffuse calibrators
            class ContinueI(Exception):
                pass
            continue_i = ContinueI()
            try:
                good_flux = 0
                # TODO the hard flux cuts here should scale with obs frequency
                for subcal in cal[cluster_idxs]:
                    # if compact - use all flux
                    if (subcal['Flux_ratio'] < 5):
                        good_flux += subcal['Total_flux']
                    # if not compact - two possibilities: 1. actually extended (low rms, discard). 2. iono smeared (high rms, calibrate)
                    elif (subcal['Peak_flux'] > 0.1*to_60_MHz) and (subcal['Isl_rms'] > 3 * global_rms):
                        good_flux += subcal['Total_flux']
                    else: # otherwise, only consider peak flux ( << total flux for true extended sources)
                        good_flux += subcal['Peak_flux'] 

                # if compact flux is present for less than 2 Jy then consider excluding it
                if (good_flux*to_60_MHz < 2) and (good_flux < 0.5*fluxes):
                    logger.debug("%s: found extended source; compact flux: %.0f%% (skip)" % (name,100*good_flux/fluxes))
                    cal['Cluster_id'][cluster_idxs] = '_'+name  # identify unused sources for debug
                    raise continue_i
            except ContinueI:
                continue

            d = lib_dd.Direction(name, soltypes=['ph-ddserial','ph1','amp1','amp2'])
            d.fluxes = fluxes
            d.spidx_coeffs = spidx_coeffs
            d.ref_freq = freq_mid
            d.localrms = localrms
            ra = np.mean(cal['RA'][cluster_idxs])
            dec = np.mean(cal['DEC'][cluster_idxs])
            d.set_position([ra, dec], phase_center)

            # skip faint directions
            if d.get_flux(60e6) < min_cal_flux60:
                logger.debug("%s: flux density @ 60 MHz: %.1f mJy (skip)" % (name, 1e3 * d.get_flux(60e6)))
                cal['Cluster_id'][cluster_idxs] = '_'+name  # identify unused sources for debug
            # skip sources that are too close to other dd-cals
            elif not lib_dd.distance_check( d, directions):
                logger.debug("%s: too close to another ddcal (skip)" % (name))
                cal['Cluster_id'][cluster_idxs] = '_'+name  # identify unused sources for debug
            # skip sources that are faint and in a low-rms region (no need to calibrate them)
            elif (d.get_flux(60e6) < 1.5) and (d.localrms < 2*global_rms) and (cmaj == 0):
                logger.debug("%s: faint source (%.1f mJy) in low-rms region (%.1f mJy) (skip)" % (name, 1e3*d.get_flux(60e6), 1e3*d.localrms))
                cal['Cluster_id'][cluster_idxs] = '_' + name  # identify unused sources for debug
            # skip if outside the mid-freq null (that should be empty)
            elif not d.is_in_region(workingReg, wcs=full_image.getWCS()):
                logger.debug("%s: outside the min-freq null region (skip)" % (name))
                cal['Cluster_id'][cluster_idxs] = '_'+name  # identify unused sources for debug
            else:
                logger.debug("%s: flux density @ 60 MHz: %.1f mJy (good)" % (name, 1e3 * d.get_flux(60e6)))
                #print('DEBUG:',name,fluxes,spidx_coeffs,gauss_area,freq_mid,size,img_beam,lsm.getColValues('MajorAxis')[idx])
                d.set_size(cal['RA'][cluster_idxs], cal['DEC'][cluster_idxs], cal['Maj'][cluster_idxs], img_beam[0]/3600)
                d.set_region(loc='ddserial/c%02i/skymodels' % cmaj)
                d.set_region_facets(facets_region_file=facetregname, loc='ddserial/c%02i/skymodels' % cmaj)
                if not d.is_in_region(peelReg, wcs=full_image.getWCS()): d.peel_off = True
                directions.append(d)

            # cap at 40 directions
            if len(directions) == 40:
                logger.warning('Too many directions found: stopping at 40.')
                break

        # take care of manually provided region to include in cal
        if manual_dd_cal != '':
            man_cal = regions.regions.Regions.read(manual_dd_cal)
            logger.warning(
                f'Using DS9 region {manual_dd_cal} for manual dd-calibrator (experimental).')
            assert len(man_cal) == 1
            # get flux in manual reg
            reg_flux = full_image.calc_flux(manual_dd_cal)
            logger.info(f'Flux in manual region: {reg_flux:.1f} Jy')
            full_residual_image = fits.open(full_image.imagename.replace('MFS-image','MFS-residual'))  # use residual so we can get data for rms
            hdr, data_res = lib_img.flatten(full_residual_image)
            wcs = WCS(hdr)
            for i, d in enumerate(directions):  # remove dd-cals that are already in region
                ra, dec = d.position
                sc = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
                if man_cal[0].contains(sc, wcs):
                    logger.info(f'{d.name} containted in manual region, remove auto-found direction.')
                    directions.pop(i)
            name = 'ddcal-' + manual_dd_cal.split('.')[0]
            d = lib_dd.Direction(name)
            d.fluxes = reg_flux
            d.spidx_coeffs = -0.8
            d.ref_freq = freq_mid
            # get rms in manual reg

            d.localrms = 15* np.std(man_cal[0].to_pixel(wcs).to_mask().cutout(data_res))
            ra, dec = man_cal[0].center.ra.to_value('deg'), man_cal[0].center.dec.to_value('deg')
            d.set_position([ra, dec], phase_center)
            d.set_size([ra], [dec], [man_cal[0].radius.to_value('deg')], img_beam[0] / 3600)
            d.set_region(loc='ddserial/c%02i/skymodels' % cmaj)
            d.set_region_facets(facets_region_file=facetregname, loc='ddserial/c%02i/skymodels' % cmaj)
            directions.insert(0, d)
            
        # create a concat region for debugging
        os.system('cat ddserial/c%02i/skymodels/ddcal[0-9][0-9][0-9][0-9].reg > ddserial/c%02i/skymodels/all-c%02i.reg' % (cmaj,cmaj,cmaj))
        # save catalogue for debugging
        cal.write('ddserial/c%02i/skymodels/cat-c%02i.fits' % (cmaj,cmaj), format='fits', overwrite=True)

        # order directions from the fluxiest one
        directions = [x for _, x in sorted(zip([d.get_flux(freq_mid) for d in directions],directions))][::-1]

        logger.info(f'Found {len(directions)} cals brighter than {min_cal_flux60} Jy (expected at 60 MHz):')
        logger.info(f'Global rms: {global_rms*1000:.2f} mJy')
        for d in directions:
            if not d.peel_off:
                logger.info('%s: flux: %.2f Jy (rms:%.2f mJy)' % (d.name, d.get_flux(freq_mid), d.localrms*1e3))
            else:
                logger.info('%s: flux: %.2f Jy (rms: %.2f mJy - peel off)' % (d.name, d.get_flux(freq_mid), d.localrms*1e3))

        pickle.dump( directions, open( picklefile, "wb" ) )
    else:
        directions = pickle.load( open( picklefile, "rb" ) )

    # do a full predict of the sky (if c>0, corrupted with solutions of prev cycle)
    with w.if_todo('c%02i-fullpredict' % cmaj):
        # wsclean predict - from ddparallel in cycle 0, otherwise from previous iteration
        logger.info('Predict full model...')
        s.add(f'wsclean -predict -padding 1.8 -name {full_image.root} -j {s.max_cpucores} -channels-out {ch_out} \
                -facet-regions {facetregname} -apply-facet-solutions {interp_h5parm} {correct_for} \
                -apply-facet-beam -use-differential-lofar-beam -facet-beam-update 120 -no-solution-directions-check \
                -reorder -parallel-reordering 4 {MSs.getStrWsclean()}',
              log='wscleanPRE-c' + str(cmaj) + '.log', commandType='wsclean')
        s.run(check=True)

        if cmaj == 0:
            logger.info('Corrupt MODEL_DATA with correct subfield solutions...')
            # LOFAR_ddparallel cycle 1 dd-solutions are on top of cycle 0 subfield solutions while data are corrected for cycle 1 subfield solutions.
            # Take this into account!
            # corrupt:
            MSs.run(f'DP3 {parset_dir}/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                      cor.parmdb=ddparallel/solutions/cal-tec-sf-merged-c0.h5 cor.correction=phase000 cor.invert=False', log='$nameMS_sf-correct.log', commandType='DP3')
            # correct:
            MSs.run(f'DP3 {parset_dir}/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                      cor.parmdb=ddparallel/solutions/cal-tec-sf-merged-c1.h5 cor.correction=phase000 cor.invert=True', log='$nameMS_sf-correct.log', commandType='DP3')

    ### DONE

    with w.if_todo('c%02i-fullsub' % cmaj):
        if cmaj == 0:
            # first cycle no corrections - ms:CORRECTED_DATA = DATA
            # cmaj>0 might have peeled sources, so we need to start from CORRECTED_DATA
            logger.info('Set CORRECTED_DATA = DATA...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = DATA"',
                    log='$nameMS_taql.log', commandType='general')
        # subtract - ms:SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA
        logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')

        ### TESTTESTTEST: empty image
        if develop and not os.path.exists('img/empty-init-c%02i-image.fits' % (cmaj)):
            clean('init-c%02i' % (cmaj), MSs, size=(fwhm*1.5, fwhm*1.5), res='normal', empty=True)
        imagenameEMPTY = 'img/wideDebug-empty-init'
        logger.info('Cleaning (empty with no sols image for debug)...')
        lib_util.run_wsclean(s, 'wscleanEMPTY-c' + str(cmaj) + '.log', MSs.getStrWsclean(), name=imagenameEMPTY,
                             data_column='SUBTRACTED_DATA',
                             size=imgsizepix, scale=str(pixscale) + 'arcsec', weight='briggs -0.5', niter=1,
                             gridder='wgridder',
                             parallel_gridding=32, minuv_l=30, mgain=0.85, parallel_deconvolution=1024,
                             join_channels='', fit_spectral_pol=3,
                             channels_out=str(ch_out), deconvolution_channels=3, pol='i',
                             no_update_model_required='', nmiter=12, auto_threshold=2.0, auto_mask=3.0,
                             local_rms='', local_rms_window=50, local_rms_strength=0.75,
                             concat_mss=True)
        os.system('mv %s-MFS-image.fits ddserial/c%02i/images' % (imagenameEMPTY, cmaj))
    ### DONE

    for dnum, d in enumerate(directions):

        logger.info('Working on direction %s/%s: %s (%f Jy - %f deg)' % (dnum+1, len(directions), d.name, d.get_flux(freq_mid), d.size))
        if d.size > 0.5: logger.warning('Patch size large: %f' % d.size)
        logstring = 'c%02i-%s' % (cmaj, d.name)

        # either faint sources that were not detected before or residuals of peeled sources - skip?
        if cmaj > 0 and d.peel_off:
            logger.info('This sources is far in the outkirts - skip.')
            continue

        with w.if_todo('%s-predict' % logstring):

            logger.info('Predict model...')
            # Predict - ms:MODEL_DATA
            s.add(f'wsclean -predict -padding 1.8 -name {full_image.root} -j {s.max_cpucores} -channels-out {ch_out} \
                -apply-facet-beam -use-differential-lofar-beam -facet-beam-update 120 \
                -facet-regions {d.get_region_facets()} -no-solution-directions-check -apply-facet-solutions {interp_h5parm} {correct_for} \
                -reorder -parallel-reordering 4 {MSs.getStrWsclean()}',
                log='wscleanPRE-'+logstring+'.log', commandType='wsclean')
            s.run(check=True)
            
            if cmaj == 0:
                logger.info('Corrupt MODEL_DATA with correct subfield solutions...')
                # LOFAR_ddparallel cycle 1 dd-solutions are on top of cycle 0 subfield solutions. Take this into account!
                MSs.run(f'DP3 {parset_dir}/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                      cor.parmdb=ddparallel/solutions/cal-tec-sf-merged-c0.h5 cor.correction=phase000 cor.invert=False',
                    log='$nameMS_sf-correct.log', commandType='DP3') # corrupt
                MSs.run(f'DP3 {parset_dir}/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA  \
                      cor.parmdb=ddparallel/solutions/cal-tec-sf-merged-c1.h5 cor.correction=phase000 cor.invert=True',
                    log='$nameMS_sf-correct.log', commandType='DP3') # correct

            # Add back the model previously subtracted for this dd-cal
            logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA + MODEL_DATA...')
            MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA + MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')
    
            ### TTESTTESTTEST: empty image but with the DD cal
            if develop and not os.path.exists('img/empty-butcal-%s-image.fits' % (logstring)):
                 clean('butcal-%s' % (logstring), MSs, size=(fwhm*1.5,fwhm*1.5), res='normal', empty=True)
        ### DONE

        # Determine parameters for correction of closest ddparallel facet
        closest = lib_h5.get_closest_dir(ddparallel_h5parm, d.position)
        with w.if_todo('%s-shift' % logstring):

            logger.info(f'Phase shift, apply phase of closest facet ({closest}) and avg...')
            lib_util.check_rm('mss-dir')
            os.makedirs('mss-dir')

            # Determine parameters for averaging - ms:SUBTRACTED_DATA -> ms-dir:DATA (->8/16/32 s and 1 chan every 2 SBs: tot of 60 or 120 chan)
            if d.get_flux(freq_mid) > 10: avgtimeint = int(round(8/timeint))
            elif d.get_flux(freq_mid) > 4: avgtimeint = int(round(16/timeint))
            elif d.get_flux(freq_mid) > 1: avgtimeint = int(round(32/timeint))
            else: avgtimeint = int(round(64/timeint))
            if d.size > 0.1/3600: # region larger than 0.1 deg -> average less to avoid smearing
                avgfreqint = int(round(MSs.getListObj()[0].getNchan() / MSs.getChout(size=2*0.195312e6))) # avg to 1 ch every 2 SBs
            else:
                avgfreqint = int(round(MSs.getListObj()[0].getNchan() / MSs.getChout(size=4*0.195312e6)))  # avg to 1 ch every 4 SBs
            if not (avgfreqint == 8 or avgfreqint == 16 or avgfreqint == 32):
                logger.warning('Strange averaging of channels (%i): %i -> %i' % (avgfreqint,MSs.getListObj()[0].getNchan(),int(MSs.getListObj()[0].getNchan()/avgfreqint)))
            MSs.run(f'DP3 {parset_dir}/DP3-shiftcorravg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=SUBTRACTED_DATA msout.datacolumn=DATA \
                    avg.timestep={avgtimeint} avg.freqstep={avgfreqint} cor.parmdb={ddparallel_h5parm} cor.correction=phase000 cor.direction={closest} \
                    shift.phasecenter=[{d.position[0]}deg,{d.position[1]}deg]', log='$nameMS_shiftcorravg-'+logstring+'.log', commandType='DP3')

            # save some info for debug
            d.avg_t = avgtimeint
            d.avg_f = avgfreqint
        ### DONE

        MSs_dir = lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s, check_flags=False )

        with w.if_todo('%s-flag' % logstring):
            logger.info('Flag on mindata...')
            MSs_dir.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')
        ### DONE

        # Correct for beam in that direction
        if not d.peel_off:
            with w.if_todo('%s-beamcorr' % logstring):
                logger.info('Correcting beam...')
                # Convince DP3 that DATA is corrected for the beam in the phase centre and correct for the new direction
                MSs_dir.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=DATA \
                        setbeam.direction=['+str(phase_center[0])+'deg,'+str(phase_center[1])+'deg] \
                        corrbeam.direction=['+str(d.position[0])+'deg,'+str(d.position[1])+'deg] corrbeam.invert=True',
                        log='$nameMS_beam-'+logstring+'.log', commandType='DP3')
            ### DONE

        # apply init - closest ddparallel sol - dd-phases
        with w.if_todo('%s-preimage' % logstring):
            logger.info('%s (pre): imaging...' % (d.name))
            clean('%s-pre' % logstring, MSs_dir, res='normal', size=[d.size,d.size], masksigma=5, imagereg=userReg)
        ### DONE

        # get initial noise and set iterators for timeint solutions
        image = lib_img.Image('img/ddserialM-%s-pre-MFS-image.fits' % logstring, userReg=userReg)
        d.rms_noise_init = image.getNoise(); rms_noise_pre = d.rms_noise_init
        d.mm_ratio_init = image.getMaxMinRatio(); mm_ratio_pre = d.mm_ratio_init
        d.set_model(image.root, typ='pre', apply_region=False)  # current best model

        doamp = False
        # usually there are 3600/32=112 or 3600/16=225 or 3600/8=450 timesteps and \
        # 60 (halfband)/120 (fullband) chans, try to use multiple numbers
        iter_ph_solint = lib_util.Sol_iterator([4, 2, 1])  # 32 or 16 or 8 * [4,2,1] s
        iter_amp1_solint = lib_util.Sol_iterator([30, 20, 10])  # 32 or 16 or 8 * [30,20,10] s
        iter_amp2_solint = lib_util.Sol_iterator([120, 60])
        iter_ph_soltype = 'diagonalphase' if (d.get_flux(freq_mid) > 5 and cmaj > 0) else 'scalarphase'
        datause = 'dual' if iter_ph_soltype == 'diagonalphase' else 'single'
        logger.info('RMS noise (init): %f' % (rms_noise_pre))
        logger.info('MM ratio (init): %f' % (mm_ratio_pre))

        d.fractional_flag_init = MSs_dir.meanFractionalFlag()
        logger.info(f'Mean fractional flag (init): {100*d.fractional_flag_init:.1f}%')

        for cdd in range(20):

            logger.info('c%02i - %s: Starting dd cycle: %02i' % (cmaj, d.name, cdd))
            logstringcal = logstring+'-cdd%02i' % cdd          

            ################################################################
            # Calibrate
            solint_ph = next(iter_ph_solint)
            dir_timeint = MSs_dir.getListObj()[0].getTimeInt()
            solint_ph_intermediate = int(min(8*solint_ph, 300/dir_timeint)) # no more than 5 minutes
            d.add_h5parm('ph-ddserial', 'ddserial/c%02i/solutions/cal-ph-ddserial%s.h5' % (cmaj,logstringcal) )
            d.add_h5parm('ph1', 'ddserial/c%02i/solutions/cal-ph1-%s.h5' % (cmaj,logstringcal) )
            if doamp:
                solint_amp1 = next(iter_amp1_solint)
                solch_amp1 = int(round(MSs_dir.getListObj()[0].getNchan() / ch_out / 4)) # TEST /4
                d.add_h5parm('amp1', 'ddserial/c%02i/solutions/cal-amp1-%s.h5' % (cmaj,logstringcal) )
                solint_amp2 = next(iter_amp2_solint)
                d.add_h5parm('amp2', 'ddserial/c%02i/solutions/cal-amp2-%s.h5' % (cmaj,logstringcal) )
            else:
                # not necessary but cleaner
                d.add_h5parm('amp1', None )
                d.add_h5parm('amp2', None )
   
            with w.if_todo('%s-calibrate' % logstringcal):
                
                # TODO: if to peel add sol.model_weighted_constraints=true

                # Smoothing - ms:DATA -> ms:SMOOTHED_DATA
                MSs_dir.run_Blsmooth(logstr=f'smooth-{logstringcal}')

                # Calibration - ms:SMOOTHED_DATA
                # Fast phase solutions
                logger.info('Distant RS phase calibration (solint: %i)...' % solint_ph)
                MSs_dir.run('DP3 '+parset_dir+'/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph-fast.h5 \
                            sol.mode='+iter_ph_soltype+' sol.datause='+datause+' sol.solint='+str(solint_ph)+' sol.smoothnessconstraint=6e6 sol.smoothnessreffrequency=54e6 \
                            sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS201LBA,CS301LBA,CS401LBA,CS501LBA,CS103LBA,CS302LBA]]',
                            log='$nameMS_solGphslow-'+logstringcal+'.log', commandType='DP3')
                # reset solutions for CS (anyways zero) and inner RS
                lib_util.run_losoto(s, 'ph-fast', [ms+'/cal-ph-fast.h5' for ms in MSs_dir.getListStr()],
                                    [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-resetph-close+mid.parset', parset_dir+'/losoto-plot-ph-fast.parset'], plots_dir='ddserial/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                # correct ph - ms:DATA -> ms:CORRECTED_DATA
                logger.info('Correct ph-fast...')
                MSs_dir.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                             cor.parmdb=cal-ph-fast.h5 cor.correction=phase000',
                            log='$nameMS_correct-'+logstringcal+'.log', commandType='DP3')

                # Smoothing - ms:CORRECTED_DATA -> ms:SMOOTHED_DATA
                MSs_dir.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-{logstringcal}')
                logger.info('Inner RS phase calibration (solint: %i)...' % (solint_ph_intermediate))
                MSs_dir.run('DP3 '+parset_dir+'/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph-slow.h5 \
                            sol.mode='+iter_ph_soltype+' sol.datause='+datause+' sol.solint='+str(solint_ph_intermediate)+' sol.smoothnessconstraint=10e6 sol.smoothnessreffrequency=54e6 \
                            sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS201LBA,CS301LBA,CS401LBA,CS501LBA,CS103LBA,CS302LBA]]',
                            log='$nameMS_solGphslow-'+logstringcal+'.log', commandType='DP3')
                # reset solutions for inner CS
                lib_util.run_losoto(s, 'ph-slow', [ms+'/cal-ph-slow.h5' for ms in MSs_dir.getListStr()],
                    [parset_dir+'/losoto-refph.parset', parset_dir+'/losoto-plot-ph-slow.parset'], plots_dir='ddserial/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                # correct ph - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                logger.info('Correct ph-slow...')
                MSs_dir.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                             cor.parmdb=cal-ph-slow.h5 cor.correction=phase000',
                            log='$nameMS_correct-'+logstringcal+'.log', commandType='DP3')

                # merge the individual h5parms.
                logger.info('Merge solutions (ph-fast + ph-slow -> ph-ddserial...')
                pol_param = '--no_pol' if iter_ph_soltype == 'scalarphase' else ''
                s.add('h5_merger.py --h5_out cal-ph-ddserial.h5 --h5_tables cal-ph-fast.h5 cal-ph-slow.h5 --h5_time_freq cal-ph-fast.h5 \
                      --no_antenna_crash %s' % (pol_param), log='h5_merger.log', commandType='python' )
                s.run(check=True)
                lib_util.run_losoto(s, 'ph-ddserial', 'cal-ph-ddserial.h5',
                                    [f'{parset_dir}/losoto-plot-ph-ddserial.parset'],
                                    plots_dir='ddserial/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                os.system('mv cal-ph-ddserial.h5 %s' % d.get_h5parm('ph-ddserial'))

                if doamp:
                   # Smoothing - ms:CORRECTED_DATA -> ms:SMOOTHED_DATA
                    MSs_dir.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-{logstringcal}')
                    logger.info('Gain amp calibration 1 (solint: %i, solch: %i)...' % (solint_amp1, solch_amp1))
                    # Calibration - ms:CORRECTED_DATA
                    MSs_dir.run('DP3 '+parset_dir+'/DP3-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-amp1.h5 \
                        sol.mode=diagonal sol.solint='+str(solint_amp1)+' sol.nchan='+str(solch_amp1)+' \
                        sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]]', \
                        log='$nameMS_solGamp1-'+logstringcal+'.log', commandType='DP3')

                    if d.peel_off:
                        losoto_parsets = [parset_dir+'/losoto-clip.parset', parset_dir+'/losoto-plot-amp1.parset']
                    else:
                        losoto_parsets = [parset_dir+'/losoto-norm.parset', parset_dir+'/losoto-plot-amp1.parset']
                    lib_util.run_losoto(s, 'amp1', [ms+'/cal-amp1.h5' for ms in MSs_dir.getListStr()], losoto_parsets,
                        plots_dir='ddserial/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                    os.system('mv cal-amp1.h5 %s' % d.get_h5parm('amp1'))

                    logger.info('Correct amp 1...')
                    # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                    MSs_dir.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                        cor.parmdb='+d.get_h5parm('amp1')+' cor.correction=amplitude000',
                        log='$nameMS_correct-'+logstringcal+'.log', commandType='DP3')

                   # Smoothing - ms:CORRECTED_DATA -> ms:SMOOTHED_DATA
                    MSs_dir.run_Blsmooth('CORRECTED_DATA', logstr=f'smooth-{logstringcal}')

                    logger.info('Gain amp calibration 2 (solint: %i)...' % solint_amp2)
                    # Calibration - ms:CORRECTED_DATA
                    MSs_dir.run('DP3 '+parset_dir+'/DP3-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp2.h5 \
                        sol.mode=diagonal sol.solint='+str(solint_amp2)+' sol.smoothnessconstraint=10e6 sol.model_weighted_constraints=true',
                        log='$nameMS_solGamp2-'+logstringcal+'.log', commandType='DP3')

                    if d.peel_off:
                        losoto_parsets = [parset_dir+'/losoto-clip2.parset', parset_dir+'/losoto-plot-amp2.parset']
                    else:
                        losoto_parsets = [parset_dir+'/losoto-norm.parset', parset_dir+'/losoto-plot-amp2.parset']

                    # ugly workaround to remove broken losoto plots in case of only 1 solution in time is present
                    if ((len(MSs_dir.getListObj()) == 1) and (solint_amp2 >= MSs_dir.getListObj()[0].getNtime())): losoto_parsets.pop()

                    lib_util.run_losoto(s, 'amp2', [ms+'/cal-amp2.h5' for ms in MSs_dir.getListStr()], losoto_parsets,
                        plots_dir='ddserial/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                    os.system('mv cal-amp2.h5 %s' % d.get_h5parm('amp2'))

                    logger.info('Correct amp 2...')
                    # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                    MSs_dir.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                        cor.parmdb='+d.get_h5parm('amp2')+' cor.correction=amplitude000',
                        log='$nameMS_correct-'+logstringcal+'.log', commandType='DP3')
            ### DONE

            # check flags
            d.fractional_flag.append(MSs_dir.meanFractionalFlag())
            logger.info(f'Mean fractional flag: {100*d.fractional_flag[-1]:.1f}% (initial: {100*d.fractional_flag_init:.1f}%)')
            if d.fractional_flag[-1] == 1.:
                logger.info('Breaking because the fractional flag increased too much...')
                break

            ###########################################################################
            # Imaging
            with w.if_todo('%s-image' % logstringcal):

                logger.info('%s (cdd: %02i): imaging...' % (d.name, cdd))
                clean('%s' % logstringcal, MSs_dir, res='normal', size=[d.size,d.size], imagereg=userReg)
            ### DONE

            image = lib_img.Image('img/ddserialM-%s-MFS-image.fits' % logstringcal, userReg=userReg)
            # something went wrong during last imaging, break
            if not os.path.exists(image.imagename):
                logger.warning('Breaking because the imaging diverged...')
                break

            # get noise, if larger than prev cycle: break
            rms_noise = image.getNoise()
            mm_ratio = image.getMaxMinRatio()
            d.add_rms_mm(rms_noise, mm_ratio) # track values for debug
            logger.info('RMS noise (cdd:%02i): %f' % (cdd,rms_noise))
            logger.info('MM ratio (cdd:%02i): %f' % (cdd,mm_ratio))

            if (cdd == 0 ) and ((rms_noise > 1.01 * rms_noise_pre and mm_ratio < 0.99 * mm_ratio_pre) or (rms_noise > 1.05 * rms_noise_pre)):
                logger.warning('Image quality decreased at cdd00! Possibly problematic ddcal.')
                # Predict - ms:MODEL_DATA (could also use wsclean but this seems easier)
                MSs_dir.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={d.get_model("pre")}-sources.txt',   log='$nameMS_pre-'+logstring+'.log', commandType='DP3')
                continue
            # if in cdd01 and cdd02 (before using amps) noise and mm ratio are worse or noise is a lot higher -> go back to previous ph_solint and restore current best model...
            # TODO technically if the cycle after we continue is converged, we would use the skipped cycle solutions since those are the "-2" ones... Fix that somehow.
            elif (cdd in [1,2,3]) and ((rms_noise > 1.01 * rms_noise_pre and mm_ratio < 0.99 * mm_ratio_pre) or (rms_noise > 1.05 * rms_noise_pre)):
                logger.warning(f'Image quality decreased after cdd{cdd}. Restore previous iteration ph_solint and model...')
                iter_ph_solint = lib_util.Sol_iterator([2*solint_ph]) # go from 1 to 2 or 2 to 4 and don't decrease further.
                # Predict - ms:MODEL_DATA (could also use wsclean but this seems easier)
                try:
                    MSs_dir.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={d.get_model("best")}-sources.txt',   log='$nameMS_pre-'+logstring+'.log', commandType='DP3')
                except:
                    # this happens if both ddcycle 0 and ddcycle 1 didn't improve
                    MSs_dir.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.sourcedb={d.get_model("pre")}-sources.txt',   log='$nameMS_pre-'+logstring+'.log', commandType='DP3')
                continue
            # if noise incresed and mm ratio decreased - or noise increased a lot!
            elif (cdd >= 4) and ((rms_noise > 0.99*rms_noise_pre and mm_ratio < 1.01*mm_ratio_pre) or rms_noise > 1.2*rms_noise_pre):
                # if (mm_ratio < 10 and cdd >= 2) or \
                # (mm_ratio < 20 and cdd >= 3) or \
                # (cdd >= 4):
                logger.debug('BREAK ddcal self cycle with noise: %f (noise_pre: %f) - mmratio: %f (mmratio_pre: %f)' % (rms_noise,rms_noise_pre,mm_ratio,mm_ratio_pre))
                break

            if (d.peel_off or cdd >= 3) and ((d.get_flux(freq_mid) > 1 and mm_ratio >= 30) or (d.get_flux(freq_mid) > 5)) and solve_amp:
                logger.debug('START AMP WITH MODE 1 - flux: %f - mmratio: %f - dist: %f' % (d.get_flux(freq_mid), mm_ratio, d.dist_from_centre))
                doamp = True
            # correct more amp in the outskirts
            elif (d.peel_off or cdd >= 3) and (d.dist_from_centre >= fwhm/4.) and ((d.get_flux(freq_mid) > 1 and mm_ratio >= 25) or d.get_flux(freq_mid) > 3) and solve_amp:
                logger.debug('START AMP WITH MODE 2 - flux: %f - mmratio: %f - dist: %f' % (d.get_flux(freq_mid), mm_ratio, d.dist_from_centre))
                doamp = True
            elif (d.peel_off or cdd >= 3) and (d.dist_from_centre >= fwhm/2.) and ((d.get_flux(freq_mid) > 1 and mm_ratio >= 20) or d.get_flux(freq_mid) > 2) and solve_amp:
                logger.debug('START AMP WITH MODE 3 - flux: %f - mmratio: %f - dist: %f' % (d.get_flux(freq_mid), mm_ratio, d.dist_from_centre))
                doamp = True

            d.set_model(image.root, typ='best', apply_region=False)  # current best model
            rms_noise_pre = rms_noise
            mm_ratio_pre = mm_ratio

        # End calibration cycle
        ##################################

        # if died the first cycle or diverged
        if cdd == 0 or ((rms_noise_pre >= d.rms_noise_init) and (mm_ratio_pre/2 < d.mm_ratio_init)):
            
            d.converged = False
            logger.warning('%s: something went wrong during the first self-cal cycle or noise did not decrease.' % (d.name))
            with w.if_todo('%s-subtract' % logstring):
                # Remove the ddcal to clean up the SUBTRACTED_DATA
                logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA')
                MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"',
                        log='$nameMS_taql.log', commandType='general')
            ### DONE

        # converged
        else:
            d.converged = True
            logger.info('%s: converged.' % d.name)
            # copy in the ddcal dir the best model
            model_skymodel = 'ddserial/c%02i/skymodels/%s-best-source.txt' % (cmaj, d.name)
            os.system('cp %s %s' % (d.get_model('best')+'-sources.txt', model_skymodel))

            # restrict to initial mask
            logger.info('Restrict model to initial region')
            os.system('cp '+d.model['best']+'-mask.fits'+' '+d.model['best']+'-mask-restricted.fits')
            lib_img.blank_image_reg(d.model['best']+'-mask-restricted.fits', d.get_region(), inverse=True, blankval=0.)
            lsm = lsmtool.load(model_skymodel)
            lsm.select('%s == True' % (d.model['best']+'-mask-restricted.fits'))
            lsm.write(model_skymodel, format='makesourcedb', clobber=True)

            # merge the final ddcal-h5parm with the nearest ddparallel facet
            with w.if_todo('%s-merge_h5' % logstring):
                # split out only the closest direction from the ddparallel h5parm
                logger.info(f'Split out dir {closest} of ddparallel solutions...')
                with h5parm(ddparallel_h5parm) as h5:
                    patches = h5.getSolset('sol000').getSou()
                p_idx = np.argwhere(np.array(list(patches.keys())) == closest) # get index of closest facet
                assert len(p_idx) == 1 # sanity check
                p_idx = p_idx[0]
                # name of split h5parm, need to remove square brackets for system commands...
                split_h5 = f'ddserial/c0{cmaj}/solutions/'+ddparallel_h5parm.split('/')[-1].replace('.h5','_'+closest.replace('[','').replace(']','')+'.h5')
                s.add(f'h5_merger.py --h5_out {split_h5} --h5_tables {ddparallel_h5parm} --filter_directions [{p_idx}] --no_pol --no_antenna_crash --no_weight_prop',
                    log='h5_merger.log', commandType='python')
                s.run(check=True)
                # change the direction of the closest facet to be exactly identical with the ddcal
                lib_h5.repoint_h5dir(split_h5, 'Dir00', d.position)
                logger.info(f'Merge final solutions for {d.name} with {closest}...')
                s.add(f'h5_merger.py --h5_out {d.get_h5parm("ph1", -2)} --h5_tables {d.get_h5parm("ph-ddserial", -2)} {split_h5} --h5_time_freq {split_h5} \
                       --no_antenna_crash', log='h5_merger.log', commandType='python')
                s.run(check=True)
                # plot the merged h5parm for debug
                lib_util.run_losoto(s, 'ph1', d.get_h5parm('ph1', -2),
                                    [f'{parset_dir}/losoto-plot-ph1.parset'],
                                    plots_dir='ddserial/c%02i/plots/plots-%s' % (cmaj,logstring))

            # remove the DD-cal from original dataset using new solutions
            with w.if_todo('%s-subtract' % logstring):

                # Predict - ms:MODEL_DATA
                logger.info('Add best model to MODEL_DATA...')
                MSs.run('DP3 '+parset_dir+'/DP3-predict.parset msin=$pathMS pre.sourcedb='+model_skymodel,
                    log='$nameMS_pre-'+logstring+'.log', commandType='DP3')

                # Store FLAGS - just for sources to peel as they might be visible only for a fraction of the band
                if d.peel_off:
                    MSs.run('taql "update $pathMS set FLAG_BKP = FLAG"',
                             log='$nameMS_taql.log', commandType='general')

                # Corrput now model - ms:MODEL_DATA -> MODEL_DATA
                # note that these are all scalar or diagonal, so they commute
                logger.info('Corrupt ph...')
                MSs.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                            cor.invert=False cor.parmdb='+d.get_h5parm('ph1',-2)+' cor.correction=phase000',
                            log='$nameMS_corrupt-'+logstring+'.log', commandType='DP3')
        
                if not d.get_h5parm('amp1',-2) is None:
                    logger.info('Corrupt amp...')
                    MSs.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                           cor.invert=False cor.parmdb='+d.get_h5parm('amp1',-2)+' cor.correction=amplitude000',
                           log='$nameMS_corrupt-'+logstring+'.log', commandType='DP3')
                if not d.get_h5parm('amp2',-2) is None:
                    MSs.run('DP3 '+parset_dir+'/DP3-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                            cor.invert=False cor.parmdb='+d.get_h5parm('amp2',-2)+' cor.correction=amplitude000',
                            log='$nameMS_corrupt-'+logstring+'.log', commandType='DP3')

                # sources to peel were not beam corrected when splitting
                if not d.peel_off:
                    # Corrupt for the beam
                    logger.info('Corrupting beam...')
                    # Convince DP3 that MODELDATA is corrected for the beam in the dd-cal direction, so I can corrupt.
                    # Then correct for the element to align with the phase centre situation
                    MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                           setbeam.direction=['+str(d.position[0])+'deg,'+str(d.position[1])+'deg] \
                           corrbeam.direction=['+str(d.position[0])+'deg,'+str(d.position[1])+'deg] corrbeam.invert=False', \
                           log='$nameMS_beam-'+logstring+'.log', commandType='DP3')
                    MSs.run('DP3 '+parset_dir+'/DP3-beam2.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                           corrbeam.direction=['+str(phase_center[0])+'deg,'+str(phase_center[1])+'deg] corrbeam.beammode=element corrbeam.invert=True', \
                           log='$nameMS_beam-'+logstring+'.log', commandType='DP3')

                if d.peel_off:
                    # Set MODEL_DATA = 0 where data are flagged, then unflag everything
                    MSs.run('taql "update $pathMS set MODEL_DATA[FLAG] = 0"',
                            log='$nameMS_taql.log', commandType='general')
                    # Restore of FLAGS
                    MSs.run('taql "update $pathMS set FLAG = FLAG_BKP"',
                            log='$nameMS_taql.log', commandType='general')

                    # if it's a source to peel, remove it from the data column used for imaging
                    logger.info('Source to peel: set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA')
                    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', \
                        log='$nameMS_taql.log', commandType='general')

                # Remove the ddcal again
                logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA')
                MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"',
                        log='$nameMS_taql.log', commandType='general')

            ### DONE

        ### TTESTTESTTEST: empty image
        if develop and not os.path.exists('img/empty-%02i-%s-image.fits' % (dnum, logstring)):
             clean('%02i-%s' % (dnum, logstring), MSs, size=(fwhm*1.5,fwhm*1.5), res='normal', empty=True)
        ###

    #####################################################
    # print a debug table
    logger.info("################################################")
    for d in directions:
        if d.peel_off and cmaj == 0:
            logger.info("### Direction (PEEL!): %s -- %.2f Jy" % (d.name, np.sum(d.fluxes)))
        elif d.peel_off and cmaj == 1:
            continue
        else:
            logger.info("### Direction: %s -- %.2f Jy" % (d.name, np.sum(d.fluxes)))
        logger.info("- Averaging: %i time - %i freq" % (d.avg_t, d.avg_f))
        logger.info("- Converged: %s" % str(d.converged))
        logger.info('init: Rms: %f, MMratio: %f' % (d.rms_noise_init,d.mm_ratio_init))
        for ic, (rms_noise, mm_ratio, fractional_flag) in enumerate(zip(d.rms_noise,d.mm_ratio,d.fractional_flag)):

            tables_to_print = '['
            for sol_type in ['ph1','amp1','amp2']: # fr
                if d.get_h5parm(sol_type, pos=ic) is not None:
                    tables_to_print += sol_type+','
            tables_to_print = tables_to_print[:-1] + ']'

            if ic == len(d.rms_noise)-2 and d.converged:
                logger.info('%02i: Rms: %f, MMratio: %f (flag:%.0f%%) - %s ***' % (ic,rms_noise,mm_ratio,100*fractional_flag,tables_to_print))
            else:
                logger.info('%02i: Rms: %f, MMratio: %f (flag:%.0f%%) - %s' % (ic,rms_noise,mm_ratio,100*fractional_flag,tables_to_print))

        # replace color in the region file to distinguish region that converged from those that didn't
        if d.converged:
            for line in fileinput.input('ddserial/c%02i/skymodels/all-c%02i.reg' % (cmaj,cmaj), inplace=True):
                if d.name in line:
                    if d.get_h5parm('amp1',-2) is not None:
                        print(line.replace('color=red','color=green'), end='')
                    else:
                        print(line.replace('color=red','color=yellow'), end='')
                else:
                    print(line, end='')

    logger.info("################################################")

    ######################################################
    # full imaging

    # be sure not to use flagged MS as ddf doesn't like them
    MSs = lib_ms.AllMSs(glob.glob('mss-avg/TC*[0-9].MS'), s, check_flags=True)
    
    imagename = 'img/wideDD-c%02i' % (cmaj)
    maskname = imagename+'_mask.fits'
    facetregname = 'ddserial/c%02i/solutions/facets-c%02i.reg' % (cmaj, cmaj)

    # combine the h5parms
    h5parms = {'ph':[], 'amp1':[], 'amp2':[]}
    for d in directions:
        # only those who converged
        if d.peel_off:
            continue
        if not d.converged:
            continue

        h5parms['ph'].append(d.get_h5parm('ph1',-2))
        if d.get_h5parm('amp1',-2) is not None:
            h5parms['amp1'].append(d.get_h5parm('amp1',-2))
        if d.get_h5parm('amp2',-2) is not None:
            h5parms['amp2'].append(d.get_h5parm('amp2',-2))

        log = '%s: Phase (%s)' % (d.name, d.get_h5parm('ph1',-2))
        log += ' Amp1 (%s)' % (d.get_h5parm('amp1',-2))
        log += ' Amp2 (%s)' % (d.get_h5parm('amp2',-2))
        logger.info(log)

    # If we hava amplitudes (should normally be the case), apply them during facet imaging
    if len(h5parms['amp1']) != 0: correct_for = 'phase000,amplitude000'
    else: correct_for = 'phase000'

    # update interp_h5parm to current cycle
    interp_h5parm = 'ddserial/c%02i/solutions/interp.h5' % cmaj
    with w.if_todo('c%02i-interpsol' % cmaj):
        logger.info("Imaging - preparing solutions:")

        for typ, h5parm_list in h5parms.items():
            # rename direction
            for h5parmFile in h5parm_list:
                dirname = h5parmFile.split('-')[3]
                lib_h5.repoint(h5parmFile, dirname)
                if typ == 'ph':
                    lib_h5.addpol(h5parmFile, 'phase000')
                    # this is a workaround for the different order of axes we get depending on scalar or diag phase solve, this will make h5parm_interpolator run smoothely.
                    lib_h5.reorder_axes(h5parmFile, ['time', 'ant', 'dir', 'freq', 'pol'], 'phase000')
                    # TODO here we have unit amplitude solutions from the h5parm_merger... delete them?
                    lib_h5.addpol(h5parmFile, 'amplitude000')
                    lib_h5.reorder_axes(h5parmFile, ['time', 'ant', 'dir', 'freq', 'pol'], 'amplitude000')
                    #lib_h5.addpol(h5parmFile, 'amplitude000')
                    s.add('losoto -v '+h5parmFile+' '+parset_dir+'/losoto-refph.parset ', log='h5parm_collector.log', commandType='python' )
                    s.run(check=True)
                    # reset high-res amplitudes in ph-solve
                    #s.add('losoto -v '+h5parmFile+' '+parset_dir+'/losoto-resetamp.parset ', log='h5parm_collector.log', commandType='python' )
                    #s.run(check=True)
                if typ == 'amp1' or typ == 'amp2':
                    lib_h5.reorder_axes(h5parmFile, ['time', 'ant', 'dir', 'freq', 'pol'], 'phase000')
                    lib_h5.reorder_axes(h5parmFile, ['time', 'ant', 'dir', 'freq', 'pol'], 'amplitude000')
                    s.add('losoto -v '+h5parmFile+' '+parset_dir+'/losoto-resetph.parset ', log='h5parm_collector.log', commandType='python' )
                    s.run(check=True)
    
        lib_util.check_rm(interp_h5parm)
        logger.info('Interpolating solutions...')
        s.add('H5parm_interpolator.py -o '+interp_h5parm+' '+' '.join(h5parms['ph']+h5parms['amp1']+h5parms['amp2']), log='h5parm_interpolator.log', commandType='python' )
        s.run(check=True)

    with w.if_todo('c%02i-imaging' % cmaj):

        logger.info('Preparing region file...')
        s.add('ds9_facet_generator.py --ms '+MSs.getListStr()[0]+' --h5 '+interp_h5parm+' --imsize '+str(int(1.1*imgsizepix))+' \
                --pixelscale '+str(pixscale)+' --writevoronoipoints --output '+facetregname,
                log='facet_generator.log', commandType='python')
        s.run(check=True)

        # non-circ beam at low dec
        if phase_center[1] < 23:
            logger.info(f'Low-declination observation ({phase_center[1]}deg). Use non-circular PSF')
            beam_kwargs = {}
        else:
            beam_kwargs = {'beam_size': 15}

        # masking
        #s.add('breizorro.py -t 3.0 -r %s -b 50 -o %s' % (full_image.imagename, maskname), 
        #        log='makemask-'+str(cmaj)+'.log', commandType='python' )
        #s.run(check=True)        

        # if defined, add userReg to the mask
        #if userReg != '': lib_img.blank_image_reg(maskname, userReg, blankval = 1.)

        # HE: What is optimal choice of subimage size and parallel gridding? Is cleaning to 3sigma enough?
        # TODO: remove fits mask and use only automasking (fits_mask=maskname)
        logger.info('Cleaning...')
        lib_util.run_wsclean(s, 'wsclean-c'+str(cmaj)+'.log', MSs.getStrWsclean(), name=imagename, data_column='CORRECTED_DATA',
                size=imgsizepix, scale=str(pixscale)+'arcsec', weight='briggs -0.5', niter=1000000, gridder='wgridder',
                parallel_gridding=len(h5parms['ph']), minuv_l=30, mgain=0.85, parallel_deconvolution=1024, join_channels='', fit_spectral_pol=3,
                channels_out=str(ch_out), deconvolution_channels=3,  multiscale='',  multiscale_scale_bias=0.65, pol='i',
                save_source_list='', no_update_model_required='',  nmiter=12, auto_threshold=2.0, auto_mask=3.0,
                apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='', facet_regions=facetregname,
                apply_facet_solutions=f'{interp_h5parm} {correct_for}', local_rms='', local_rms_window=50, local_rms_strength=0.75,
                concat_mss=True, **beam_kwargs)
 
        os.system(f'mv {imagename}*MFS-image*fits {imagename}*MFS-model*fits {imagename}*MFS-residual*fits \
                  {imagename}-0*image*fits {imagename}-0*model*fits ddserial/c{cmaj:02}/images')

    ### DONE

    full_image = lib_img.Image('ddserial/c%02i/images/%s-MFS-image.fits' % (cmaj, imagename.split('/')[-1]), userReg=userReg)
    full_image.nantozeroModel()
    # min_cal_flux60 *= 0.8  # go deeper

##############################################################################################################
### Calibration finished - additional images with scientific value

# Low res as this is relevant only for transient detection
with w.if_todo('output-timedep'):
    logger.info('Cleaning (time dep images)...')
    for tc, msfile in enumerate(MSs.getListStr()):
        imagenameT = 'img/wideDD-TC%02i-c%02i' % (tc, cmaj)
        lib_util.run_wsclean(s, 'wscleanTC'+str(tc)+'-c'+str(cmaj)+'.log', msfile, concat_mss=True, name=imagenameT, data_column='CORRECTED_DATA',
                size=int(imgsizepix/4), scale=str(pixscale*4)+'arcsec', taper_gaussian='60arcsec', weight='briggs 0', niter=1000000, gridder='wgridder',
                parallel_gridding=len(h5parms['ph']), minuv_l=30, mgain=0.85, parallel_deconvolution=512, join_channels='', fit_spectral_pol=3,
                channels_out=str(ch_out), deconvolution_channels=3,  multiscale='',  multiscale_scale_bias=0.65, pol='i',
                no_update_model_required='',  nmiter=12, auto_threshold=2.0, auto_mask=3.0,
                apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='', facet_regions=facetregname,
                apply_facet_solutions=f'{interp_h5parm} {correct_for}', local_rms='', local_rms_window=50, local_rms_strength=0.75,
                beam_size=15)

        os.system('mv %s-MFS-image*.fits %s-MFS-residual.fits ddserial/c%02i/images' % (imagenameT, imagenameT, cmaj))
### DONE

with w.if_todo('output-vstokes'):
    imagenameV = 'img/wideDD-v-c%02i' % (cmaj)
    logger.info('Cleaning (V-stokes)...')
    lib_util.run_wsclean(s, 'wscleanV-c'+str(cmaj)+'.log', MSs.getStrWsclean(), concat_mss=True, name=imagenameV, data_column='CORRECTED_DATA', size=int(imgsizepix/4), scale=str(pixscale*4)+'arcsec',
                taper_gaussian='60arcsec', weight='briggs 0', niter=1000000, gridder='wgridder', parallel_gridding=6, no_update_model_required='', minuv_l=30, mgain=0.85, parallel_deconvolution=512,
                auto_threshold=3.0, join_channels='', fit_spectral_pol=3, channels_out=6, deconvolution_channels=3,
                pol='v')

    os.system('mv %s-MFS-image*.fits %s-MFS-residual.fits ddserial/c%02i/images' % (imagenameV, imagenameV, cmaj))
### DONE

with w.if_todo('output-lres'):
    imagenameL = 'img/wideDD-lres-c%02i' % (cmaj)
    logger.info('Cleaning (low res)...')
    lib_util.run_wsclean(s, 'wscleanLR-c'+str(cmaj)+'.log', MSs.getStrWsclean(), concat_mss=True, name=imagenameL, data_column='CORRECTED_DATA',
                size=int(imgsizepix/4), scale=str(pixscale*4)+'arcsec', weight='briggs 0', taper_gaussian='60arcsec', niter=1000000, gridder='wgridder',
                parallel_gridding=len(h5parms['ph']), minuv_l=20, mgain=0.85, parallel_deconvolution=512, join_channels='', fit_spectral_pol=3,
                channels_out=str(ch_out), deconvolution_channels=3,  multiscale='',  multiscale_scale_bias=0.65, pol='i',
                no_update_model_required='',  nmiter=12, auto_threshold=2.0, auto_mask=3.0,
                apply_facet_beam='', facet_beam_update=120, use_differential_lofar_beam='', facet_regions=facetregname,
                apply_facet_solutions=f'{interp_h5parm} {correct_for}', local_rms='', local_rms_window=50, local_rms_strength=0.75, beam_size=60)

    os.system('mv %s-MFS-image*.fits %s-MFS-residual.fits ddserial/c%02i/images' % (imagenameL, imagenameL, cmaj))
### DONE

with w.if_todo('output-lressub'):

    # now make a low res and source subtracted map for masking extended sources
    logger.info('Predicting DD-corrupted...')
    s.add('wsclean -predict -padding 1.8 -name '+full_image.root+' -j '+str(s.max_cpucores)+' -channels-out '+str(ch_out)+' \
                    -apply-facet-beam -use-differential-lofar-beam -facet-beam-update 120 \
                    -facet-regions '+facetregname+' \
                    -apply-facet-solutions '+interp_h5parm+' amplitude000,phase000 \
                    -reorder -parallel-reordering 4 '+MSs.getStrWsclean(),
                    log='wscleanPRE4LR-c'+str(cmaj)+'.log', commandType='wsclean')
    s.run(check=True)

    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"',
        log='$nameMS_taql.log', commandType='general')

    imagenameLS = 'img/wideDD-lressub-c%02i' % (cmaj)
    logger.info('Cleaning (low res sub)...')
    lib_util.run_wsclean(s, 'wscleanLRS-c'+str(cmaj)+'.log', MSs.getStrWsclean(), concat_mss=True, name=imagenameLS, data_column='SUBTRACTED_DATA', size=int(imgsizepix/4), scale=str(pixscale*4)+'arcsec',
                weight='briggs 0', niter=1000000, gridder='wgridder', parallel_gridding=len(h5parms['ph']), no_update_model_required='', minuv_l=20, mgain=0.85, parallel_deconvolution=512,
                auto_threshold=3.0, join_channels='', fit_spectral_pol=3, channels_out=str(ch_out), deconvolution_channels=3,
                multiscale='', multiscale_scale_bias=0.65, pol='i', taper_gaussian='60arcsec',
                apply_facet_beam='', use_differential_lofar_beam='', facet_beam_update=120, facet_regions=facetregname, apply_facet_solutions=f'{interp_h5parm} {correct_for}')
 
    os.system('mv %s-MFS-image*.fits %s-MFS-residual.fits ddserial/c%02i/images' % (imagenameLS, imagenameLS, cmaj))
### DONE

with w.if_todo('output-debugempty'):
    imagenameEMPTY = 'img/wideDebug-empty-c%02i' % (cmaj)
    logger.info('Cleaning (empty with no sols image for debug)...')
    lib_util.run_wsclean(s, 'wscleanEMPTY-c'+str(cmaj)+'.log', MSs.getStrWsclean(), name=imagenameEMPTY, data_column='SUBTRACTED_DATA',
                size=imgsizepix, scale=str(pixscale)+'arcsec', weight='briggs -0.5', niter=1, gridder='wgridder',
                parallel_gridding=32, minuv_l=30, mgain=0.85, parallel_deconvolution=1024, join_channels='', fit_spectral_pol=3,
                channels_out=str(ch_out), deconvolution_channels=3, pol='i',
                no_update_model_required='',  nmiter=12, auto_threshold=2.0, auto_mask=3.0,
                local_rms='', local_rms_window=50, local_rms_strength=0.75,
                concat_mss=True)
    os.system('mv %s-MFS-image.fits ddserial/c%02i/images' % (imagenameEMPTY, cmaj))
### DONE

with w.if_todo('output_PB'):
    logger.info('Make primary beam...')
    s.add('makepb.py -o ddserial/primarybeam.fits -s 10 -p 120 %s' % MSs.getStrWsclean(), log='makepb.log', commandType='python')
    s.run(check=True)
### DONE

# remove unwanted columns
with w.if_todo('remove_col'):
    logger.info('Removing unwanted columns...')
    MSs.run('taql "update $pathMS set DATA = CORRECTED_DATA"', log='$nameMS_taql.log', commandType='general')
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CORRECTED_DATA,FLAG_BKP,MODEL_DATA"', log='$nameMS_delcol.log', commandType='python')
### DONE

w.alldone()
