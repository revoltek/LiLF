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
parset_dir = parset.get('LOFAR_dd-incr','parset_dir')
maxniter = parset.getint('LOFAR_dd','maxniter')
calFlux = parset.getfloat('LOFAR_dd','calFlux')
userReg = parset.get('model','userReg')

####################################################
MSs_self = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs_self.getListObj()[0].getPhaseCentre()
fwhm = MSs_self.getListObj()[0].getFWHM(freq='mid')

##########################
logger.info('Cleaning...')
lib_util.check_rm('ddcal')
os.makedirs('ddcal/masks')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/solutions')
os.makedirs('ddcal/skymodels')

def clean(p, MSs, size, res='normal', apply_beam=False):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.getListObj()[0].getResolution() 
    if res == 'normal':
        pixscale = float('%.1f'%(pixscale/2.))
    elif res == 'high':
        pixscale = float('%.1f'%(pixscale/4.))
    elif res == 'low':
        pass # no change

    imsize = [0,0]
    imsize[0] = int(size[0]*1.05/(pixscale/3600.)) # add 5%
    imsize[1] = int(size[1]*1.05/(pixscale/3600.)) # add 5%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 64: imsize[0] == 64
    if imsize[1] < 64: imsize[1] == 64

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    if res == 'normal':
        weight = 'briggs 0'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.7'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/ddcal-'+str(p)
    lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.8, \
            baseline_averaging=5, parallel_deconvolution=256, auto_threshold=3, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 4)

    # clean 2
    # TODO: add -parallel-deconvolution when source lists can be saved (https://sourceforge.net/p/wsclean/tickets/141/)
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/ddcalM-'+str(p)
    if apply_beam:

        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, save_source_list='', size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.8, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            multiscale='', \
            auto_threshold=0.5, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

        logger.info('Cleaning V ('+str(p)+')...')
        imagename = 'img/ddcalV-'+str(p)
        lib_util.run_wsclean(s, 'wscleanV-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imgsize, scale=str(pixscale)+'srcsec', pol='V', \
            weight='briggs 0.', niter=1000, no_update_model_required='', minuv_l=30, maxuv_l=5000, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400)

    else:

        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
            weight=weight, niter=50000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.8, \
            multiscale='', \
            auto_threshold=0.5, fits_mask=im.maskname, \
            baseline_averaging=5, join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')


############################################################
# TODO: use SUBTRACTED_DATA (no pre-correction) or CORRECTED_DATA (DIE iono correction)?
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs_self.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg)
mosaic_image.selectCC()

# TEST:
#mosaic_image = lib_img.Image('ddcal/images/c00/mos-MFS-image.fits', userReg = userReg)
rms_noise_pre = np.inf

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)
    if c>=1: directions_old = directions
    directions = []

    lib_util.check_rm('img')
    os.makedirs('img')
    if not os.path.exists('ddcal/masks/regions-c%02i' % c): os.makedirs('ddcal/masks/regions-c%02i' % c)
    if not os.path.exists('ddcal/images/c%02i' % c): os.makedirs('ddcal/images/c%02i' % c)
    mask_voro = 'ddcal/masks/facets%02i.fits' % c
    if c>=1: mask_voro_old = 'ddcal/masks/facets%02i.fits' % (c-1)

    ### TTESTTESTTEST: DIE image
    #if c == 0:
    #    clean('init', MSs, size=(fwhm,fwhm), res='normal')
    ###

    ### group into patches corresponding to the mask islands
    mask_cl = mosaic_image.imagename.replace('image.fits', 'mask-cl.fits')
    # this mask is with no user region, done isolate only bight compact sources
    if not os.path.exists(mask_cl): 
        mosaic_image.makeMask(threshisl=7, atrous_do=False, remove_extended_cutoff=0.001, maskname=mask_cl)
    
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group(mask_cl, root='Isl')
    # this removes all sources not in the mask-cl
    lsm.select('Patch = Isl.*', useRegEx=True)
    # this regroup sources
    x = lsm.getColValues('RA',aggregate='wmean')
    y = lsm.getColValues('Dec',aggregate='wmean')
    flux = lsm.getColValues('I',aggregate='sum')
    grouper = lib_dd.Grouper(zip(x,y),flux)
    grouper.run()
    clusters = grouper.grouping()
    grouper.plot()
    os.system('mv grouping*png ddcal/skymodels/')
    patchNames = lsm.getPatchNames()

    logger.info('Merging nearby sources...')
    for cluster in clusters:
        patches = patchNames[cluster]
        #print ('merging:', cluster, patches)
        if len(patches) > 1:
            lsm.merge(patches.tolist())

    lsm.select('I >= %f Jy' % calFlux, aggregate='sum')

    # keep track of CC names used for calibrators so not to subtract them afterwards
    cal_names = lsm.getColValues('Name')

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
    
    # write file
    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
    skymodel_cl_plot = 'ddcal/masks/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # convert to blob
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_cl_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)
    del lsm
    
    ### select the rest of the sources to be subtracted
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    names = lsm.getColValues('Name')
    lsm.remove( np.array([ i for i, name in enumerate(names) if name in cal_names ]) )
    lsm.ungroup()
    logger.info("Total flux in rest field %i Jy" % np.sum(lsm.getColValues('I')) )
    
    # when possible regroup in patches using old DD-calibrators
    if c>=1: lsm.group('facet', facet=mask_voro_old, root='Isl_patch')

    # write file
    skymodel_rest = 'ddcal/skymodels/skymodel%02i_rest.txt' % c
    lsm.write(skymodel_rest, format='makesourcedb', clobber=True)
    skymodel_rest_plot = 'ddcal/masks/skymodel%02i_rest.png' % c
    lsm.plot(fileName=skymodel_rest_plot, labelBy='patch')
       
    # convert to blob
    skymodel_rest_skydb = skymodel_rest.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_rest_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_rest, skymodel_rest_skydb), log='makesourcedb_rest.log', commandType='general')
    s.run(check=True)
    
    ### create regions (using cluster directions)
    logger.info("Create regions.")
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, \
            outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/masks/voronoi%02i.png' % c)
    lsm.group('facet', facet=mask_voro, root='Isl_patch')
    [ d.add_mask_voro(mask_voro) for d in directions ]

    # write file
    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
    skymodel_voro_plot = 'ddcal/masks/skymodel%02i_voro.png' % c
    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')

    # convert to blob
    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_voro_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
    s.run(check=True)
    del lsm

    logger.debug("Islands' info:")
    for i, d in enumerate(directions):
        logger.info("%s: Flux=%f (coord: %s - size: %s deg)" % ( d.name, d.flux_cal, str(d.position_cal), str(d.size) ) )

    ################################################################
    # Calibrate TEC
    logger.info('Subtraction rest_field...')

    # Fist cycle remove other sources with no DD corrections, when DD correction is available use it
    if c == 0:
        # Predict - ms:MODEL_DATA
        logger.info('Add rest_field to MODEL_DATA...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb,log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')
    
        # Empty dataset from faint sources
        logger.info('Set SUBTRACTED_DATA = DATA - MODEL_DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    
    else:

        # Copy DATA -> SUBTRACTED_DATA
        logger.info('Set SUBTRACTED_DATA = DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        for i, d in enumerate(directions_old):
            # predict - ms:MODEL_DATA
            logger.info('Patch '+d.name+': predict...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb+' pre.sources='+d.name, \
                    log='$nameMS_pre0-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
            # corrupt (note: use previous cal table) - ms:MODEL_DATA -> ms:MODEL_DATA
            logger.info('Patch '+d.name+': corrupt...')
            MSs.run('DPPP '+parset_dir+'/DPPP-corrupt2.parset msin=$pathMS \
                    corC.parmdb=ddcal/solutions/cal-core-c'+str(c-1)+'.h5   corC.direction=['+d.name+'] \
                    corR.parmdb=ddcal/solutions/cal-remote-c'+str(c-1)+'.h5 corR.direction=['+d.name+']', \
                    log='$nameMS_corrupt0-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
            
            logger.info('Patch '+d.name+': subtract...')
            MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

    ### TESTTESTTEST: empty image with cals
    #MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    #clean('onlycals-c'+str(c), MSs, size=(fwhm,fwhm), res='normal')
    ###

    # Smoothing - ms:SUBTRACTED_DATA -> ms:SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -f 1.0 -r -i SUBTRACTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')    

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Core calibration...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS msin.baseline="[CR]*&&;!RS210LBA;!RS310LBA;!RS509LBA;!RS508LBA;!RS409LBA;!RS208LBA;!RS307LBA" \
            sol.antennaconstraint=[[CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA]] \
            sol.solint=5 sol.nchan=4 sol.h5parm=$pathMS/cal-core-c'+str(c)+'.h5 sol.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solTECcore-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    for MS in MSs.getListObj():
        lib_util.run_losoto(s, 'core-c'+str(c)+'-'+MS.nameMS, MS.pathMS+'/cal-core-c'+str(c)+'.h5', \
                [parset_dir+'/losoto-jump.parset', parset_dir+'/losoto-resetremote.parset', parset_dir+'/losoto-plot-tec.parset'])
    os.system('mv plots-core-c'+str(c)+'* ddcal/plots')
    s.add('H5parm_collector.py -V -s sol000 -o ddcal/solutions/cal-core-c'+str(c)+'.h5 '+' '.join(glob.glob('cal-core-c'+str(c)+'-*.h5')),\
                            log='losoto-collector-c'+str(c)+'.log', commandType="python", processors='max')
    s.run(check = True)
    lib_util.check_rm('cal-core-c'+str(c)+'-*.h5')

    # Calibration - ms:SMOOTHED_DATA
    logger.info('Remote calibration...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solTEC.parset msin=$pathMS \
            sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS406LBA,RS407LBA,RS503LBA]] \
            sol.applycal.parmdb=ddcal/solutions/cal-core-c'+str(c)+'.h5 sol.applycal.correction=tec000 \
            sol.solint=1 sol.nchan=4 sol.h5parm=$pathMS/cal-remote-c'+str(c)+'.h5 sol.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solTECremote-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    for MS in MSs.getListObj():
        lib_util.run_losoto(s, 'remote-c'+str(c)+'-'+MS.nameMS, MS.pathMS+'/cal-remote-c'+str(c)+'.h5', \
                [parset_dir+'/losoto-jump.parset', parset_dir+'/losoto-plot-tec.parset'])
    os.system('mv plots-remote-c'+str(c)+'* ddcal/plots')
    s.add('H5parm_collector.py -V -s sol000 -o ddcal/solutions/cal-remote-c'+str(c)+'.h5 '+' '.join(glob.glob('cal-remote-c'+str(c)+'-*.h5')),\
                            log='losoto-collector-c'+str(c)+'.log', commandType="python", processors='max')
    s.run(check = True)
    lib_util.check_rm('cal-remote-c'+str(c)+'-*.h5')


    sys.exit()
    
    ##############################################################
    # low S/N DIE corrections
    if c>=0:
        logger.info('DIE calibration...')
        # predict and corrupt each facet
        logger.info('Reset MODEL_DATA...')
        MSs.run('taql "update $pathMS set MODEL_DATA = 0"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
       
        for i, d in enumerate(directions):
            # predict - ms:MODEL_DATA
            logger.info('Patch '+d.name+': predict...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS msout.datacolumn=MODEL_DATA_DIR pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
    
            # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
            logger.info('Patch '+d.name+': corrupt...')
            MSs.run('DPPP '+parset_dir+'/DPPP-corrupt2.parset msin=$pathMS msin.datacolumn=MODEL_DATA_DIR msout.datacolumn=MODEL_DATA_DIR \
                    corC.parmdb=ddcal/solutions/cal-core-c'+str(c)+'.h5   corC.direction=['+d.name+'] \
                    corR.parmdb=ddcal/solutions/cal-remote-c'+str(c)+'.h5 corR.direction=['+d.name+']', \
                log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
            logger.info('Patch '+d.name+': subtract...')
            MSs.run('taql "update $pathMS set MODEL_DATA = MODEL_DATA + MODEL_DATA_DIR"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

        # Smoothing - ms:DATA -> ms:SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -f 1.0 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')    

        # Convert to circular - SMOOTHED_DATA -> SMOOTHED_DATA
        logger.info('Converting to circular...')
        MSs.run('mslin2circ.py -i $pathMS:SMOOTHED_DATA -o $pathMS:SMOOTHED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)

        # FR Calibration - ms:SMOOTHED_DATA
        logger.info('Solving DIE FR...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.h5parm=$pathMS/calG1-c'+str(c)+'.h5', \
                log='$nameMS_solG1-c'+str(c)+'.log', commandType='DPPP')
    
        # Plot solutions
        lib_util.run_losoto(s, 'G1-c'+str(c), [MS+'/calG1-c'+str(c)+'.h5' for MS in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-fr.parset'])
        os.system('mv plots-G1-c'+str(c)+'* ddcal/plots')

        # Correct DIE FR - ms:DATA -> CORRECTED_DATA
        logger.info('DIE FR correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor1.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-G1-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
               log='$nameMS_corFR-c'+str(c)+'.log', commandType='DPPP')

        # Smoothing - ms:CORRECTED_DATA -> ms:SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -f 1.0 -r -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')    
    
        # DIE Calibration - ms:SMOOTHED_DATA
        logger.info('Solving DIE AMP...')
        MSs.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS sol.h5parm=$pathMS/calG2-c'+str(c)+'.h5', \
                log='$nameMS_solG2-c'+str(c)+'.log', commandType='DPPP')
    
        # Plot solutions
        lib_util.run_losoto(s, 'G2-c'+str(c), [MS+'/calG2-c'+str(c)+'.h5' for MS in MSs.getListStr()], \
                [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-amp.parset'])
        os.system('mv plots-G2-c'+str(c)+'* ddcal/plots')

        # Correct DIE AMP - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('DIE AMP correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor1.parset msin=$pathMS cor.parmdb=cal-G2-c'+str(c)+'.h5 cor.correction=amplitudeSmooth', \
               log='$nameMS_corAMP-c'+str(c)+'.log', commandType='DPPP')

        # Copy CORRECTED_DATA -> SUBTRACTED_DATA
        logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

        ### TESTTESTTEST: init image with DIE correction
        #clean('die-c'+str(c), MSs, size=(fwhm,fwhm), res='normal')
        ###

    else:
        # Copy DATA -> SUBTRACTED_DATA
        logger.info('Set SUBTRACTED_DATA = DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')

    ###########################################################
    # Empty the dataset

    logger.info('Subtraction...')
    for i, d in enumerate(directions):
        
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d.name+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
    
        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+d.name+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt2.parset msin=$pathMS \
                    corC.parmdb=ddcal/solutions/cal-core-c'+str(c)+'.h5   corC.direction=['+d.name+'] \
                    corR.parmdb=ddcal/solutions/cal-remote-c'+str(c)+'.h5 corR.direction=['+d.name+']', \
                log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
        logger.info('Patch '+d.name+': subtract...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

    ### TESTTESTTEST: empty image
    MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    clean('empty-c'+str(c), MSs, size=(fwhm,fwhm), res='normal')
    ###

    ###########################################################
    # Add back 
    logger.info('Facet imaging...')
    for i, d in enumerate(directions):
        #TODO: see if we can phase shift and average before predict-corrupt=correct
        # predict - ms:MODEL_DATA
        logger.info('Patch '+d.name+': predict...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                   log='$nameMS_pre2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
        logger.info('Patch '+d.name+': corrupt...')
        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt2.parset msin=$pathMS \
                 corC.parmdb=ddcal/solutions/cal-core-c'+str(c)+'.h5   corC.direction=['+d.name+'] \
                 corR.parmdb=ddcal/solutions/cal-remote-c'+str(c)+'.h5 corR.direction=['+d.name+']', \
                 log='$nameMS_corrupt2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': add...')
        MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+d.name+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor2.parset msin=$pathMS \
                corC.parmdb=ddcal/solutions/cal-core-c'+str(c)+'.h5   corC.direction=['+d.name+'] \
                corR.parmdb=ddcal/solutions/cal-remote-c'+str(c)+'.h5 corR.direction=['+d.name+']', \
                log='$nameMS_cor-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

        logger.info('Patch '+d.name+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(d.position_facet[0])+'deg,'+str(d.position_facet[1])+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
        
        logger.info('Patch '+d.name+': imaging...')
        clean(d.name, lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size, apply_beam = c==maxniter )

        # TEST: if one wants to make a low-res patch
        if c>=2:
            logger.info('Patch '+d.name+': imaging high-res...')
            clean(d.name+'-high', lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size, res='high')
            logger.info('Patch '+d.name+': predict high-res...')
            # predict - ms:MODEL_DATA
            s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % \
                    ('img/ddcalM-'+d.name+'-high-sources.txt', 'img/ddcalM-'+d.name+'-high-sources.skydb'), log='makesourcedb_'+d.name+'.log', commandType='general' )
            s.run(check=True)
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb=img/ddcalM-'+d.name+'-high-sources.skydb pre.sources='+d.name, \
                    log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
            logger.info('Patch '+d.name+': subtract high-res...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')
            logger.info('Patch '+d.name+': imaging low-res...')
            clean(d.name+'-low', lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size, res='low', apply_beam = c==maxniter )

    ##############################################################
    # Mosaiching

    # reorder in increasing isl_num order
    isl_nums = [d.isl_num for d in directions]
    directions = [d for _, d in sorted(zip(isl_nums,directions))]

    for d in directions:
        d.image = lib_img.Image('img/ddcalM-%s-MFS-image.fits' % d.name, userReg = userReg)
        d.image_res = lib_img.Image('img/ddcalM-%s-MFS-residual.fits' % d.name, userReg = userReg)
        d.image_low = lib_img.Image('img/ddcalM-%s-low-MFS-image.fits' % d.name, userReg = userReg)
        d.image_high = lib_img.Image('img/ddcalM-%s-high-MFS-image.fits' % d.name, userReg = userReg)

        # restrict skymodel to facet
        d.image.makeMask(threshisl=5)
        d.image.selectCC()
        lsm = lsmtool.load(d.image.skymodel_cut)
        lsm.group('facet', facet=mask_voro, root='Isl_patch' )
        lsm.select('Patch = Isl_patch_%i' % d.isl_num )
        lsm.write(d.image.skymodel_cut, format='makesourcedb', clobber=True)

    logger.info('Mosaic: image...')
    image_files = ' '.join([d.image.imagename for d in directions])
    mosaic_image_file = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_image_file, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    logger.info('Mosaic: residual image...')
    image_files = ' '.join([d.image_res.imagename for d in directions])
    mosaic_residual_file = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual_file, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    if c>=2:
        logger.info('Mosaic: low-res image...')
        image_files = ' '.join([d.image_low.imagename for d in directions])
        mosaic_image_low_file = 'img/mos-low-MFS-image.fits'
        s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_image_low_file, log='mosaic-img-low-c'+str(c)+'.log', commandType='python')
        s.run(check=True)
    
        logger.info('Mosaic: high-res image...')
        image_files = ' '.join([d.image_high.imagename for d in directions])
        mosaic_image_high_file = 'img/mos-high-MFS-image.fits'
        s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_image_high_file, log='mosaic-img-high-c'+str(c)+'.log', commandType='python')
        s.run(check=True)

    # prepare new skymodel
    lsm = lsmtool.load(directions[0].image.skymodel_cut)
    lsm.ungroup()
    for image in [d.image for d in directions[1:]]:
        lsm2 = lsmtool.load(image.skymodel_cut)
        lsm2.ungroup()
        lsm.concatenate(lsm2, keep='all')
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)
    del lsm

    os.system('cp img/*M*MFS-image.fits img/mos*.fits ddcal/images/c%02i' % c )
    mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)
    mosaic_image.makeMask(threshisl=3, atrous_do=True) # used in the faceting function

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = mosaic_image.getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > rms_noise_pre: break
    rms_noise_pre = rms_noise
