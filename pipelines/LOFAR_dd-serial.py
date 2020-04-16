#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re, pickle
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
logger_obj = lib_log.Logger('pipeline-dd-serial.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-dd-serial.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_dd-serial','parset_dir')
userReg = parset.get('model','userReg')
aterm_imaging = False

def clean(p, MSs, res='normal', size=[1,1]):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs.getListObj()[0].getResolution() 
    if res == 'normal':
        pixscale = float('%.1f'%(pixscale/2.5))
    elif res == 'high':
        pixscale = float('%.1f'%(pixscale/3.5))
    elif res == 'low':
        pass # no change

    imsize = [int(size[0]*1.5/(pixscale/3600.)), int(size[1]*1.5/(pixscale/3600.))] # add 50%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 256: imsize[0] = 256
    if imsize[1] < 256: imsize[1] = 256

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    if res == 'normal':
        weight = 'briggs -0.3'
        maxuv_l = None
    elif res == 'high':
        weight = 'briggs -0.6'
        maxuv_l = None
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500
    else:
        logger.error('Wrong "res": %s.' % str(res))
        sys.exit()

    # clean 1
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/ddcal-'+str(p)
    lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            baseline_averaging=5, parallel_deconvolution=512, auto_threshold=5, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 5)

    # clean 2
    logger.info('Cleaning w/ mask ('+str(p)+')...')
    imagename = 'img/ddcalM-'+str(p)
    lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, do_predict=True, \
            size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            multiscale='', multiscale_scale_bias=0.75, multiscale_scales='0,10,20,40,80', 
            baseline_averaging=5, parallel_deconvolution=512, local_rms='', auto_threshold=0.75, auto_mask=1.5, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

    os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')


#############################################################
if w.todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('ddcal')
    os.makedirs('ddcal/masks')
    os.makedirs('ddcal/plots')
    os.makedirs('ddcal/images')
    os.makedirs('ddcal/solutions')
    os.makedirs('ddcal/skymodels')

    w.done('cleaning')
### DONE

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
logger.info('Add columns...')

MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
MSs.getListObj()[0].makeBeamReg('ddcal/beam.reg', freq='mid')
beamReg = 'ddcal/beam.reg'
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg)
if not os.path.exists(mosaic_image.skymodel_cut): mosaic_image.selectCC()

for C in range(2):
    logger.info('Starting major cycle: %i' % C)
    
    if w.todo('delimg-c%02i' % C):
        lib_util.check_rm('img')
        os.makedirs('img')
        w.done('delimg-c%02i' % C)
    ### DONE

    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % C
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    #skymodel_rest = 'ddcal/skymodels/skymodel%02i_rest.txt' % C
    #skymodel_rest_skydb = skymodel_rest.replace('.txt','.skydb')
    #skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % C
    #skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')

    picklefile = 'ddcal/directions-c%02i.pickle' % C

    if not os.path.exists(picklefile):
        directions = []

        if not os.path.exists('ddcal/masks/regions-c%02i' % C): os.makedirs('ddcal/masks/regions-c%02i' % C)
        if not os.path.exists('ddcal/images/c%02i' % C): os.makedirs('ddcal/images/c%02i' % C)
        #mask_voro = 'ddcal/masks/facets%02i.fits' % C
    
        ### TTESTTESTTEST: DIE image
        #if c == 0:
        #    clean('init', MSs, size=(fwhm*1.5,fwhm*1.5), res='normal')
        ###
    
        ### group into patches corresponding to the mask islands
        mask_cl = mosaic_image.imagename.replace('image.fits', 'mask-cl.fits')
        # this mask is with no user region, done to isolate only bight compact sources
        if not os.path.exists(mask_cl): 
            mosaic_image.beamReg = 'ddcal/beam.reg'
            mosaic_image.makeMask(threshisl=7, atrous_do=False, remove_extended_cutoff=0.001, maskname=mask_cl, only_beam=True)
        
        lsm = lsmtool.load(mosaic_image.skymodel_cut)
        #lsm.group(mask_cl, root='Isl')
        lsm.group('tessellate', targetFlux=0.1, root='Isl') # test to keep all sources
        # this removes all sources not in the mask-cl
        #lsm.select('Patch = Isl.*', useRegEx=True)
        # this regroup sources
        x = lsm.getColValues('RA',aggregate='wmean')
        y = lsm.getColValues('Dec',aggregate='wmean')
        flux = lsm.getColValues('I',aggregate='sum')
        grouper = lib_dd.Grouper(list(zip(x,y)), flux, look_distance=0.2, kernel_size=0.1, grouping_distance=0.03)
        grouper.run()
        clusters = grouper.grouping()
        grouper.plot()
        os.system('mv grouping*png ddcal/plots/')
        patchNames = lsm.getPatchNames()
    
        logger.info('Merging nearby sources...')
        for cluster in clusters:
            patches = patchNames[cluster]
            print('Merging:', patches)
            if len(patches) > 1:
                lsm.merge(patches.tolist())
    
        # keep track of CC names used for calibrators so not to subtract them afterwards
        cal_names = lsm.getColValues('Name')
    
        lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
        positions = lsm.getPatchPositions()
        for name, flux, size in \
                zip( lsm.getPatchNames(), lsm.getColValues('I', aggregate='sum'), lsm.getPatchSizes(units='deg') ):
            direction = lib_dd.Direction(name)
            position = [positions[name][0].deg, positions[name][1].deg ]
            direction.set_position( position, cal=True )
            direction.set_flux(flux, cal=True)
            direction.set_size([size,size], cal=True)
            directions.append(direction)
        directions = [x for _,x in sorted(zip([d.flux_cal for d in directions],directions))][::-1] # reorder with flux

        # write file
        lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
        skymodel_cl_plot = 'ddcal/masks/skymodel%02i_cluster.png' % C
        lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')
        lsm.setColValues('name', [x.split('_')[-1] for x in lsm.getColValues('patch')]) # just for the region - this makes this lsm useless
        lsm.write('ddcal/masks/regions-c%02i/cluster.reg' % C, format='ds9', clobber=True)
        del lsm
    
        # convert to blob
        lib_util.check_rm(skymodel_cl_skydb)
        s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
        s.run(check=True)
        
        #### select the rest of the sources to be subtracted
        #lsm = lsmtool.load(mosaic_image.skymodel_cut)
        #names = lsm.getColValues('Name')
        #lsm.remove( np.array([ i for i, name in enumerate(names) if name in cal_names ]) )
        #lsm.ungroup()
        #logger.info("Total flux in rest field %i Jy" % np.sum(lsm.getColValues('I')) )
        #
        ## write file
        #lsm.write(skymodel_rest, format='makesourcedb', clobber=True)
        #skymodel_rest_plot = 'ddcal/masks/skymodel%02i_rest.png' % c
        #lsm.plot(fileName=skymodel_rest_plot, labelBy='patch')
        #   
        ## convert to blob
        #lib_util.check_rm(skymodel_rest_skydb)
        #s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_rest, skymodel_rest_skydb), log='makesourcedb_rest.log', commandType='general')
        #s.run(check=True)
        #
        #### create regions (using cluster directions)
        #logger.info("Create regions.")
        #lsm = lsmtool.load(mosaic_image.skymodel_cut)
        #lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, \
        #        outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/masks/voronoi%02i.png' % c)
        #lsm.group('facet', facet=mask_voro, root='Isl_patch')
        #[ d.add_mask_voro(mask_voro) for d in directions ]
    
        ## write file
        #lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
        #skymodel_voro_plot = 'ddcal/masks/skymodel%02i_voro.png' % c
        #lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')
        #lsm.setColValues('name', [x.split('_')[-1] for x in lsm.getColValues('patch')]) # just for the region - this makes this lsm useless
        #lsm.write('ddcal/masks/regions-c%02i/voro.reg' % c, format='ds9', clobber=True)
        #del lsm
    
        ## convert to blob
        #lib_util.check_rm(skymodel_voro_skydb)
        #s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
        #s.run(check=True)

        pickle.dump( directions, open( picklefile, "wb" ) )

    else:
        directions = pickle.load( open( picklefile, "rb" ) )

    for d in directions:
        logger.info('Working on direction: %s (%f Jy - %f deg)' % (d.name, d.flux_cal, d.size_cal[0]))

        if w.todo('%s-subtract' % d.name):
            logger.info('%s: Subtraction rest_field...' % d.name)

            patches_sub = [xd.name for xd in directions if xd.name != d.name]
    
            # Predict - ms:MODEL_DATA
            logger.info('Add rest_field to MODEL_DATA...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_cl_skydb+' pre.sources=['+','.join(patches_sub)+']', \
                    log='$nameMS_pre-'+d.name+'.log', commandType='DPPP')
            
            # Empty dataset from faint sources
            logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
            MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"', \
                    log='$nameMS_taql-'+d.name+'.log', commandType='general')
            
            w.done('%s-subtract' % d.name)
        ### DONE

        if w.todo('%s-shift' % d.name):
            logger.info('%s: phase shift and avg...' % d.name)
            
            lib_util.check_rm('mss-dir')
            os.makedirs('mss-dir')
            # Shift - ms:SUBTRACTED_DATA -> ms:DATA
            MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=SUBTRACTED_DATA \
                    shift.phasecenter=['+str(d.position_cal[0])+'deg,'+str(d.position_cal[1])+'deg\]', \
                    log='$nameMS_shift-'+d.name+'.log', commandType='DPPP')
 
            w.done('%s-shift' % d.name)
        ### DONE

        MSs_dir = lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s )

        if w.todo('%s-smooth' % d.name):
            logger.info('%s: BL-based smoothing...' % d.name)

            # Smoothing - ms:DATA -> ms:SMOOTHED_DATA
            MSs_dir.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', \
                    log='$nameMS_smooth-'+d.name+'.log', commandType='python')    
 
            w.done('%s-smooth' % d.name)
        ### DONE

        if w.todo('%s-predict' % d.name):
            logger.info('Add ddcal model to MODEL_DATA...')

            # Predict - ms:MODEL_DATA
            MSs_dir.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_cl_skydb+' pre.sources=['+d.name+']', \
                    log='$nameMS_pre2-'+d.name+'.log', commandType='DPPP')

            w.done('%s-predict' % d.name)
        ### DONE

        if w.todo('%s-preimage' % d.name):

            logger.info('%s: pre-imaging...' % d.name)
            clean('%s-pre' % d.name, MSs_dir, res='normal', size=d.size_cal)

            w.done('%s-preimage' % d.name)
        ### DONE
        
        # get initial noise
        image = lib_img.Image('img/ddcalM-%s-pre-MFS-image.fits' % d.name)
        rms_noise_pre = image.getNoise()
        logger.info('RMS noise: %f' % rms_noise_pre)

        for c in range(10):

            logger.info('%s: Starting cycle: %02i' % (d.name, c))

            ################################################################
            # Calibrate
   
            if w.todo('calibrate-%s-c%02i' % (d.name, c)):
                logger.info('%s: Calibrate...' % d.name)

                # Calibration - ms:SMOOTHED_DATA
                logger.info('Gain calibration...')
                try: solint = [20,10,5,2][c]
                except: solint = 1
                MSs_dir.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS \
                    sol.h5parm=$pathMS/cal-g.h5 sol.solint='+str(solint), \
                    log='$nameMS_solG-'+d.name+'-c'+str(c)+'.log', commandType='DPPP')
    
                # Plot solutions
                lib_util.run_losoto(s, 'g', [ms+'/cal-g.h5' for ms in MSs_dir.getListStr()], \
                    [parset_dir+'/losoto-amp.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-plot-ph.parset'])
                os.system('mv plots-g ddcal/plots/plots-g-%s-c%i' % (d.name, c))
                os.system('mv cal-g.h5 ddcal/solutions/cal-g-%s-c%i.h5' % (d.name, c))

                # correct G - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                logger.info('Correct ph...')
                MSs_dir.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS \
                             cor.parmdb=ddcal/solutions/cal-g-'+d.name+'-c'+str(c)+'.h5 cor.correction=phase000 cor.direction=['+d.name+']', \
                             log='$nameMS_correct-'+d.name+'-c'+str(c)+'.log', commandType='DPPP')
                if c>3:
                    logger.info('Correct amp...')
                    MSs_dir.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS \
                        cor.parmdb=ddcal/solutions/cal-g-'+d.name+'c'+str(c)+'.h5 cor.correction=amplitude000 cor.direction=['+d.name+']', \
                        log='$nameMS_correct-'+d.name+'-c'+str(c)+'.log', commandType='DPPP') 

                w.done('calibrate-%s-c%02i' % (d.name, c))
            ### DONE

            if w.todo('image-%s-c%02i' % (d.name, c)):

                logger.info('%s: imaging...' % d.name)
                clean('%s-c%02i' % (d.name, c), MSs_dir, res='normal', size=d.size_cal)

                w.done('image-%s-c%02i' % (d.name, c))
            ### DONE
        
            # get noise, if larger than 95% of prev cycle: break
            image = lib_img.Image('img/ddcalM-%s-c%02i-MFS-image.fits' % (d.name, c))
            rms_noise = image.getNoise()
            logger.info('RMS noise: %f' % rms_noise)
            if rms_noise > rms_noise_pre: break
            rms_noise_pre = rms_noise
