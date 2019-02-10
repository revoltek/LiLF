#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-dd.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('dd','parset_dir')
maxniter = parset.getint('dd','maxniter')
userReg = parset.get('model','userReg')

####################################################
MSs_self = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs_self.getListObj()[0].getPhaseCentre()
MSs_self.getListObj()[0].makeBeamReg('self/beam.reg', to_null=True) # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null
beamReg = 'self/beam.reg'

##########################
logger.info('Cleaning...')
lib_util.check_rm('ddcal')
os.makedirs('ddcal/masks')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/skymodels')

def clean(p, MSs, size=2., apply_beam=False):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.getListObj()[0].getResolution()/2. # TODO: weighting will lower the resolutions a bit, therefore a /2 should be enough
    imsize = int(size/(pixscale/3600.))

    if imsize < 512:
        imsize = 512

    if imsize % 2 == 1: imsize += 1 # make even

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/ddcal-'+str(p)
    lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', baseline_averaging=5, minuv_l=30, mgain=0.85, \
            auto_threshold=20, join_channels='', fit_spectral_pol=2, channels_out=8)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 3)

    # clean 2
    # To be safe don't do -continue with beam
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/ddcalM-'+str(p)
    if apply_beam:
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            auto_threshold=0.1, fits_mask=im.maskname, join_channels='', fit_spectral_pol=2, channels_out=8, save_source_list='')
    else:
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', baseline_averaging=5, minuv_l=30, mgain=0.85, \
            auto_threshold=0.1, fits_mask=im.maskname, join_channels='', fit_spectral_pol=2, channels_out=8, save_source_list='')
    os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')


############################################################
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs_self.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
#logger.info('BL-based smoothing...')
#MSs.run('BLsmooth.py -f 1.0 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')

# setup initial model
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg, beamReg = beamReg)
mosaic_image.selectCC()
rms_noise_pre = np.inf

for c in xrange(maxniter):
    logger.info('Starting cycle: %i' % c)

    lib_util.check_rm('img')
    os.makedirs('img')
    os.makedirs('ddcal/images/c%02i/regions' % c)
    mask = 'ddcal/masks/facets%02i.fits' % c

    ### group into patches of similar flux
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group('tessellate', targetFlux='20Jy', root='Dir', applyBeam=False, method = 'wmean', pad_index=True)
    lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
    directions = lsm.getPatchPositions()
    logger.info("Created %i directions." % len(directions))

    # write file
    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
    skymodel_cl_plot = 'ddcal/skymodels/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # convert to blob
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_cl_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)

    ### create regions (using cluster directions)
    logger.info("Create regions.")
    lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, outdir_reg='ddcal/images/c%02i/regions/' % c, out_mask=mask, png='ddcal/skymodels/voronoi%02i.png' % c)
    lsm.group('facet', facet=mask, root='Dir')
    lsm.setPatchPositions(method='mid') # recalculate the patch centre as mid point for imaging
    directions = lsm.getPatchPositions()
    sizes = dict(zip(lsm.getPatchNames(), lsm.getPatchSizes(units='degree')))

    # write file
    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
    skymodel_voro_plot = 'ddcal/skymodels/skymodel%02i_voro.png' % c
    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')

    # convert to blob
    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_voro_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
    s.run(check=True)

    del lsm

#    ################################################################
#    # Calibration
#    logger.info('Calibrating...')
#    MSs.run('DPPP '+parset_dir+'/DPPP-solDD.parset msin=$pathMS ddecal.h5parm=$pathMS/cal-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
#            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')
#
#    # Plot solutions
#    lib_util.run_losoto(s, 'c'+str(c), [MS+'/cal-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot.parset'])
#    os.system('mv plots-c'+str(c)+'* ddcal/plots')

    ##############################################################
    # low S/N DIE corrections
    # TODO: add amp and FR sol + correction here after ft() a DDE-corrupted model

    ###########################################################
    # Empty the dataset
    logger.info('Subtraction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-sub.parset msin=$pathMS sub.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 sub.sourcedb='+skymodel_voro_skydb, \
                   log='$nameMS_sub-c'+str(c)+'.log', commandType='DPPP')#, maxThreads=1)

    for patch, phasecentre in directions.iteritems():

        # add back single path - ms:SUBTRACTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+patch+': add back...')
        MSs.run('DPPP '+parset_dir+'/DPPP-add.parset msin=$pathMS add.applycal.parmdb=$pathMS/cal-c'+str(c)+'.h5 add.sourcedb='+skymodel_voro_skydb+' add.directions=[['+patch+']]', \
                   log='$nameMS_add-c'+str(c)+'-p'+str(patch)+'.log', commandType='DPPP')

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+patch+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+patch+']', \
               log='$nameMS_cor-c'+str(c)+'-p'+str(patch)+'.log', commandType='DPPP')

        logger.info('Patch '+patch+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(phasecentre[0].degree)+'deg,'+str(phasecentre[1].degree)+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-p'+str(patch)+'.log', commandType='DPPP')
        
        logger.info('Patch '+patch+': imaging...')
        clean(patch, lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=sizes[patch], apply_beam = c==maxniter )

    ##############################################################
    # Mosaiching

    images = []
    for image, region in zip( sorted(glob.glob('img/ddcalM-Dir*MFS-image.fits')), sorted(glob.glob('ddcal/images/c%02i/regions/Dir*' % c)) ):
        images.append( lib_img.Image(image, facetReg = region, userReg = userReg) )

    logger.info('Mosaic: image...')
    image_files = ' '.join([image.imagename for image in images])
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    logger.info('Mosaic: residuals...')
    image_files = ' '.join([image.imagename.replace('image', 'residual') for image in images])
    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    # prepare new skymodel
    skymodels = []
    for image in images:
        image.selectCC() # restrict to facet
        skymodels.append(image.skymodel_cut)
    lsm = lsmtool.load(skymodels[0])
    for skymodel in skymodels[1:]:
        lsm2 = lsmtool.load(skymodel)
        lsm.concatenate(lsm2)
    lsm.group('single')
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)

    os.system('cp img/*M*MFS-image.fits img/mos-MFS-image.fits img/mos-MFS-residual.fits ddcal/images/c%02i' % c )
    mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise
