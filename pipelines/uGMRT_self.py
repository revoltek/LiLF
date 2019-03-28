#!/usr/bin/env python
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool
from astropy.coordinates import SkyCoord
import astropy.units as u

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-self.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

parset = lib_util.getParset()
parset_dir = parset.get('uGMRT_self','parset_dir')
sourcedb = parset.get('model','sourcedb') # relative to tgts dir "./tgts/xxx/model.txt
userReg = parset.get('model','userReg') # relative to tgts dir "./tgts/xxx/region.ref"

MSs = lib_ms.AllMSs( glob.glob('mss/*.MS'), s )

############################################################################
# Clear
logger.info('Cleaning...')
lib_util.check_rm('img')
os.makedirs('img')
lib_util.check_rm('self')
os.makedirs('self/images')
os.makedirs('self/solutions')
os.makedirs('self/plots')

is_wideband = len(MSs.getListObj()[0].getFreqs()) > 1000 # if > 1000 chans, it is wideband

# make beam
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg')
beamReg = 'self/beam.reg'

#################################################################
# Get online model
if sourcedb is None:
    if not os.path.exists('tgts.skydb'):
        fwhm = MSs.getListObj()[0].getFWHM()
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "http://172.104.228.177/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f"' % (radeg, decdeg, fwhm))
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        #Reduces the flux of clean component according to a primary beam function
        #NOTE: monochromatic approximation!
        center = SkyCoord(phasecentre[0]*u.deg, phasecentre[1]*u.deg)
        sources = SkyCoord( lsm.getColValues('RA')*u.deg, lsm.getColValues('Dec')*u.deg )
        d = center.separation(sources)
        # from http://www.aips.nrao.edu/cgi-bin/ZXHLP2.PL?PBCOR (converto to arcmin and multiply by freq in GHz)
        d = d.deg * 60 * np.mean(MSs.getListObj()[0].getFreqs())/1.e9
        I = lsm.getColValues('I')
        parm = [-3.397,47.192,-30.931,7.803] # 325 MHz GMRT
        I_corr = I * (1 + (parm[0]/10**3)*d**2 + (parm[1]/10**7)*d**4 + \
             (parm[2]/10**10)*d**6 + (parm[3]/10**13)*d**8)
        lsm.setColValues('I', I_corr)
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')

    sourcedb = 'tgts.skydb'

##################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for MS in MSs.getListStr():
    lib_util.check_rm(MS+'/'+sourcedb_basename)
    logger.debug('Copy: '+sourcedb+' -> '+MS)
    os.system('cp -r '+sourcedb+' '+MS)

logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c MODEL_DATA,SUBTRACTED_DATA,CORRECTED_DATA_DIE -i DATA', log="$nameMS_addcol.log", commandType="python")

logger.info('Add model to MODEL_DATA...')
MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=$pathMS/'+sourcedb_basename, log='$nameMS_pre.log', commandType='DPPP')

#####################################################################################################
# Self-cal cycle
for c in range(3):

    logger.info('Start selfcal cycle: '+str(c))

    # solve - concat*.MS:SMOOTHED_DATA
    if not is_wideband:
        logger.info('Solving G...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/gs.h5 \
                 sol.solint=1 sol.nchan=1 sol.mode=complexgain sol.smoothnessconstraint=1e6', \
                    log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
        lib_util.run_losoto(s, 'gs'+str(c), [MS+'/gs.h5' for MS in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-flag.parset', parset_dir+'/losoto-norm.parset'])
        os.system('mv plots-gs'+str(c)+'* self/solutions/')
        os.system('mv cal-gs'+str(c)+'*.h5 self/solutions/')
        # correct phases - MS:DATA -> MS:CORRECTED_DATA
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.correction=phase000', \
                    log='$nameMS_corGp-c'+str(c)+'.log', commandType='DPPP')
        if c >= 2:
            # correct amplitudes - MS:DATA -> MS:CORRECTED_DATA
            logger.info('Correcting Ga...')
            MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=self/solutions/cal-gs'+str(c)+'.h5 cor.correction=amplitude000', \
                    log='$nameMS_corGa-c'+str(c)+'.log', commandType='DPPP')
    else:
        logger.info('Solving TEC...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/tec.h5 \
                 sol.solint=1 sol.nchan=1 sol.mode=tec', \
                    log='$nameMS_solG-c'+str(c)+'.log', commandType='DPPP')
        lib_util.run_losoto(s, 'tec'+str(c), [MS+'/tec.h5' for MS in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-tec.parset'])
        os.system('mv plots-tec'+str(c)+'* self/plots/')
        os.system('mv cal-tec'+str(c)+'*.h5 self/solutions/')
        # correct phases - MS:DATA -> MS:CORRECTED_DATA
        logger.info('Correcting Gp...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=self/solutions/cal-tec'+str(c)+'.h5 cor.correction=tec000', \
                log='$nameMS_corGp-c'+str(c)+'.log', commandType='DPPP')

    # set image size at 1.5 * FWHM
    imgsizepix = int(1.5*MSs.getListObj()[0].getFWHM()/(2./3600))
    imgsizepix += imgsizepix % 2 # make even
    if c>=2: imgsizepix *= 2 # last cycle make a very large image to catch source in the sidelobes

    # clean mask clean
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsizepix, scale='2arcsec', \
            weight='briggs 0.', niter=10000, no_update_model_required='', mgain=0.9, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=20, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl=4, atrous_do=False)
    
    # baseline averaging possible as we cut longest baselines (also it is in time, where smearing is less problematic)
    # TODO: add -parallel-deconvolution=256 when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imgsizepix, scale='2arcsec', \
            weight='briggs 0.', niter=100000, no_update_model_required='', mgain=0.9, \
            #multiscale='', multiscale_scales='0,5,10,20,40', \
            baseline_averaging=5, auto_threshold=1, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=2, channels_out=8)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')

    if c != 2:

        im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
        im.makeMask(threshisl=5, atrous_do=False)
        im.selectCC()

        # predict - ms: MODEL_DATA
        logger.info('Predict model...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA pre.sourcedb='+im.skydb, \
                log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

# Copy images
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in range(3) ]
os.system('mv img/wideM-2-sources.txt self/images')

# final, large self-cal image (used to get the skymodel)
image_field = lib_img.Image('self/images/wideM-2-MFS-image.fits', userReg=userReg)
image_field.makeMask(threshisl=5, atrous_do=True)
image_field.selectCC()
# small self-cal image (used only to define the size of the final mosaic)
image_small = lib_img.Image('self/images/wideM-1-MFS-image.fits', userReg=userReg)
image_small.makeMask(threshisl=5, atrous_do=True)

# Move DIE-corrected data into CORRECTED_DATA_DIE
logger.info('Set CORRECTED_DATA_DIE = CORRECTED_DATA...')
MSs.run('taql "update $pathMS set CORRECTED_DATA_DIE = CORRECTED_DATA"', log='$nameMS_taql2.log', commandType='general')

# TESTTESTTEST
imgsizepix = int(1.5*MSs.getListObj()[0].getFWHM()/(2./3600))
imgsizepix += imgsizepix % 2 # make even
lib_util.run_wsclean(s, 'wscleanTEST-c'+str(c)+'.log', MSs.getStrWsclean(), name='img/testinit', size=imgsizepix, scale='2arcsec', \
        data_column='CORRECTED_DATA_DIE', \
        weight='briggs 0.', niter=100000, no_update_model_required='', mgain=0.9, \
        baseline_averaging=5, auto_threshold=1, \
        join_channels='', fit_spectral_pol=2, channels_out=8)

logger.info('== Starting DDE cal ==')
logger.info('Cleaning...')
lib_util.check_rm('ddcal')
os.makedirs('ddcal/masks')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/skymodels')

#####################################################################################################
# DDE-cal cycle
directions = []
rms_noise_pre = np.inf
for c in range(3):

    logger.info('Start DDE-cal cycle: '+str(c))

    # Prepare facets
    if not os.path.exists('ddcal/masks/regions-c%02i' % c): os.makedirs('ddcal/masks/regions-c%02i' % c)
    if not os.path.exists('ddcal/images/c%02i' % c): os.makedirs('ddcal/images/c%02i' % c)
    mask_voro = 'ddcal/masks/facets%02i.fits' % c

    ### group into patches corresponding to the mask islands
    lsm = lsmtool.load(image_field.skymodel_cut)
    lsm.group(image_field.maskname, root='Isl')

    ### select bright sources
    # TODO: aggregate nearby sources
    lsm.select('I >= 0.2 Jy', aggregate='sum')
    lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
    for name, flux in zip(lsm.getPatchNames(), lsm.getColValues('I', aggregate='sum')):
        direction = lib_dd.Direction(name)
        position = [ lsm.getPatchPositions()[name][0].deg, lsm.getPatchPositions()[name][1].deg ]
        direction.set_position( position, cal=True )
        direction.set_flux(flux, cal=True)
        directions.append(direction)

    logger.info("Created %i bright sources" % len(directions))
    tot_flux = np.sum([d.flux_cal for d in directions])
    logger.info("Total flux of bright sources %i Jy" % tot_flux)

    # write file (skymodel cluster: has only cals, no facets)
    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
    skymodel_cl_plot = 'ddcal/skymodels/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # convert to blob
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_cl_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)

    ### select the rest of the sources to be subtracted
    lsm = lsmtool.load(image_field.skymodel_cut)
    lsm.group(image_field.maskname, root='Isl')
    lsm.select('I < 0.2 Jy', aggregate='sum')
    lsm.ungroup()
    rest_field = lsm.getColValues('I')
    rest_field = np.sum(rest_field)
    logger.info("Total flux in rest field %i Jy" % rest_field)

    # write file (skymodel rest: has only facets, no cals)
    skymodel_rest = 'ddcal/skymodels/skymodel%02i_rest.txt' % c
    lsm.write(skymodel_rest, format='makesourcedb', clobber=True)
    skymodel_rest_plot = 'ddcal/skymodels/skymodel%02i_rest.png' % c
    lsm.plot(fileName=skymodel_rest_plot, labelBy='patch')

    # convert to blob
    skymodel_rest_skydb = skymodel_rest.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_rest_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_rest, skymodel_rest_skydb), log='makesourcedb_rest.log', commandType='general')
    s.run(check=True)

    ### create regions (using cluster directions)
    logger.info("Create regions.")
    lsm = lsmtool.load(image_field.skymodel_cut)
    # use the cycle=1 image that is as large as the beam
    lib_dd.make_voronoi_reg(directions, image_small.maskname, \
        outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/masks/voronoi%02i.png' % c)
    # TODO: check if group ignore sources outside mask_voro
    lsm.group('facet', facet=mask_voro, root='Isl_patch')
    [ d.add_mask_voro(mask_voro) for d in directions ]

    # write file (skymodel voro: has all cals+facets)
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

    logger.debug("Island info:")
    for i, d in enumerate(directions):
        logger.info("%s: Flux=%f (coord: %s - size: %s deg)" % ( d.name, d.flux_cal, str(d.position_cal), str(d.size) ) )

    ###############################################################
    # Calibration

    #Predict - ms:MODEL_DATA
    logger.info('Add rest_field to MODEL_DATA...')
    MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb,log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')

    # Empty dataset from faint sources (TODO: better corrupt with DDE solutions when available before subtract)
    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA_DIE - MODEL_DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA_DIE - MODEL_DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')

    # TESTTESTTEST
    imgsizepix = int(1.5*MSs.getListObj()[0].getFWHM()/(2./3600))
    imgsizepix += imgsizepix % 2 # make even
    lib_util.run_wsclean(s, 'wscleanTEST-c'+str(c)+'.log', MSs.getStrWsclean(), name='img/testempty', size=imgsizepix, scale='2arcsec', \
            data_column='SUBTRACTED_DATA', \
            weight='briggs 0.', niter=100000, no_update_model_required='', mgain=0.9, \
            baseline_averaging=5, auto_threshold=1, \
            join_channels='', fit_spectral_pol=2, channels_out=8)

    # Smoothing - ms:SUBTRACTED_DATA -> ms:SMOOTHED_DATA
    # TODO: check if it makes sense
    #logger.info('BL-based smoothing...')
    #MSs.run('BLsmooth.py -f 1.0 -r -i SUBTRACTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solGdd.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA sol.h5parm=$pathMS/cal-dd-c'+str(c)+'.h5 sol.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    lib_util.run_losoto(s, 'dd-c'+str(c), [MS+'/cal-dd-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot-amp-dd.parset', parset_dir+'/losoto-plot-ph-dd.parset'])
    os.system('mv plots-dd-c'+str(c)+'* ddcal/plots')

    ###########################################################
    # Empty the dataset
    logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA_DIE...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA_DIE"', log='$nameMS_taql2-c'+str(c)+'.log', commandType='general')

    logger.info('Subtraction...')

    for i, d in enumerate(directions):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d.name+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name,log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        # TODO: corrupt also for amplitudes?
        logger.info('Patch '+d.name+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-dd-c'+str(c)+'.h5 cor.direction=['+d.name+'] cor.correction=phase000 cor.invert=False \
                 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA', \
                log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        # subtract - ms:SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA
        logger.info('Patch '+d.name+': subtract...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql3-c'+str(c)+'-'+d.name+'.log', commandType='general')

    ###########################################################
    # Facet imaging
    for i, d in enumerate(directions):
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d.name+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                   log='$nameMS_pre2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        # TODO: corrupt also for amplitudes?
        logger.info('Patch '+d.name+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-dd-c'+str(c)+'.h5 cor.direction=['+d.name+'] cor.correction=phase000 cor.invert=False \
                 msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA', \
                 log='$nameMS_corrupt2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': add...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql4-c'+str(c)+'-'+d.name+'.log', commandType='general')

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+d.name+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.parmdb=$pathMS/cal-dd-c'+str(c)+'.h5 cor.direction=['+d.name+'] cor.correction=phase000 \
               msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA', \
               log='$nameMS_cor-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(d.position_facet[0])+'deg,'+str(d.position_facet[1])+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': imaging...')

        # set pixscale and imsize
        MSs_shift = lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s )
    
        imsize = [0,0]
        imsize[0] = int(d.size[0]*1.05/(2/3600.)) # add 5%
        imsize[1] = int(d.size[1]*1.05/(2/3600.)) # add 5%
        imsize[0] += imsize[0] % 2
        imsize[1] += imsize[1] % 2
        if imsize[0] < 64: imsize[0] == 64
        if imsize[1] < 64: imsize[1] == 64
    
        logger.debug('Image size: '+str(imsize))
    
        # clean 1
        logger.info('Cleaning ('+d.name+')...')
        imagename = 'img/ddcal-'+d.name
        lib_util.run_wsclean(s, 'wscleanA-'+d.name+'.log', MSs_shift.getStrWsclean(), name=imagename, size=imsize, scale='2arcsec', \
                weight='briggs 0.', niter=10000, no_update_model_required='', minuv_l=30, mgain=0.85, \
                baseline_averaging=5, parallel_deconvolution=256, \
                auto_threshold=20, join_channels='', fit_spectral_pol=2, channels_out=8)
    
        # make mask
        im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
        im.makeMask(threshisl = 3)
    
        # clean 2
        # TODO: add -parallel-deconvolution when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
        logger.info('Cleaning w/ mask ('+d.name+')...')
        imagename = 'img/ddcalM-'+d.name
        lib_util.run_wsclean(s, 'wscleanB-'+d.name+'.log', MSs_shift.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale='2arcsec', \
                weight='briggs 0.', niter=100000, no_update_model_required='', minuv_l=30, mgain=0.85, \
                baseline_averaging=5, auto_threshold=0.1, fits_mask=im.maskname, \
                join_channels='', fit_spectral_pol=2, channels_out=8)
    
        os.system('cat logs/wscleanA-'+d.name+'.log logs/wscleanB-'+d.name+'.log | grep "background noise"')

    ##############################################################
    # Mosaiching
    isl_nums = [d.isl_num for d in directions]
    directions = [d for _, d in sorted(zip(isl_nums,directions))]

    directions_formosaic = []
    for d in directions:
        d.image = lib_img.Image('img/ddcalM-%s-MFS-image.fits' % d.name, userReg = userReg)
        d.image_res = lib_img.Image('img/ddcalM-%s-MFS-residual.fits' % d.name, userReg = userReg)
        # restrict skymodel to facet
        d.image.selectCC()
        lsm = lsmtool.load(d.image.skymodel_cut)
        if d.cal_has_facet:
            lsm.group('facet', facet=mask_voro, root='Isl_patch' )
            lsm.select('Patch = Isl_patch_%i' % d.isl_num )
            directions_formosaic.append(d)
        lsm.write(d.image.skymodel_cut, format='makesourcedb', clobber=True)

    logger.info('Mosaic: image...')
    image_files = ' '.join([d.image_res.imagename for d in directions_formosaic])
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    logger.info('Mosaic: residuals...')
    image_files = ' '.join([d.image_res.imagename for d in directions_formosaic])
    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    # prepare new skymodel
    lsm = lsmtool.load(directions[0].image.skymodel_cut)
    for image in [d.image for d in directions[1:]]:
        lsm2 = lsmtool.load(image.skymodel_cut)
        lsm.concatenate(lsm2)
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)

    os.system('cp img/*M*MFS-image.fits img/mos-MFS-image.fits img/mos-MFS-residual.fits ddcal/images/c%02i' % c )
    image_field = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)
    image_field.makeMask(threshisl=3, atrous_do=True)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise

logger.info("Done.")
