#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

# TODO: remove regions and move to masks

import sys, os, glob, re, pickle
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd, lib_h5
logger_obj = lib_log.Logger('pipeline-dd-serial.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-dd-serial.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_dd-serial','parset_dir')
userReg = parset.get('model','userReg')
maxIter = parset.getint('LOFAR_dd-serial','maxIter')
min_cal_flux60 = parset.getfloat('LOFAR_dd-serial','minCalFlux60')
removeExtendedCutoff = parset.getfloat('LOFAR_dd-serial','removeExtendedCutoff')

def clean(p, MSs, res='normal', size=[1,1], empty=False, imagereg=None):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    try:
        pixscale = MSs.getListObj()[0].getResolution() 
    except:
        logger.warning('Fail evaluate pixscale, probably the MS is fully flagged.')
        return

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

    if empty:

        logger.info('Cleaning empty ('+str(p)+')...')
        imagename = 'img/empty-'+str(p)
        lib_util.run_wsclean(s, 'wscleanE-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, data_column='SUBTRACTED_DATA', \
                size=imsize, scale=str(pixscale)+'arcsec', \
                weight=weight, niter=0, no_update_model_required='', minuv_l=30, mgain=0, \
                baseline_averaging='')
 
    else:

        # clean 1
        logger.info('Cleaning ('+str(p)+')...')
        imagename = 'img/ddcal-'+str(p)
        lib_util.run_wsclean(s, 'wscleanA-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, \
                size=imsize, scale=str(pixscale)+'arcsec', \
                weight=weight, niter=10000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
                baseline_averaging='', parallel_deconvolution=512, auto_threshold=5, \
                join_channels='', fit_spectral_pol=3, channels_out=ch_out, deconvolution_channels=3)
    
        # make mask
        im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
        try:
            im.makeMask(threshisl = 7, rmsbox=(70,5))
        except:
            logger.warning('Fail to create mask for %s.' % imagename+'-MFS-image.fits')
            return

        # restrict to inside the dd-region
        if imagereg is not None:
            lib_img.blank_image_reg(im.maskname, imagereg, inverse=True, blankval=0.,)
    
        # clean 2
        # TODO: add deconvolution_channels when bug fixed
        logger.info('Cleaning w/ mask ('+str(p)+')...')
        imagename = 'img/ddcalM-'+str(p)
        lib_util.run_wsclean(s, 'wscleanB-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, do_predict=True, \
                size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
                weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
                multiscale='', multiscale_scale_bias=0.65, multiscale_scales='0,10,20,40,80', 
                baseline_averaging='', parallel_deconvolution=512, local_rms='', auto_threshold=0.75, auto_mask=1.5, fits_mask=im.maskname, \
                join_channels='', fit_spectral_pol=3, channels_out=ch_out) #, deconvolution_channels=3)

        os.system('cat logs/wscleanA-'+str(p)+'.log logs/wscleanB-'+str(p)+'.log | grep "background noise"')


#############################################################
if w.todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('ddcal')
    os.makedirs('ddcal/init')
    lib_util.check_rm('img')
    os.makedirs('img')

    w.done('cleaning')
### DONE

# goes down to 8 seconds/4 chans
if not os.path.exists('mss-avg'):
    MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )
    timeint = MSs.getListObj()[0].getTimeInt()
    avgtimeint = int(round(10/timeint))
    nchan_init = MSs.getListObj()[0].getNchan()
    # avg (x8) sol (x6) - we need a multiple of 8x6=48, the largest that is <nchan
    # survey after avg (x8): 60, final number of sol 10
    # pointed after avg (x8): 120, final number of sol 20
    nchan = nchan_init - nchan_init%48
    os.makedirs('mss-avg')
    logger.info('Averaging in time (%is -> 10s - %ich -> %ich)' % (timeint,nchan_init,nchan))
    MSs.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-avg/$nameMS.MS msin.datacolumn=CORRECTED_DATA msin.nchan='+str(nchan)+' \
            avg.timestep='+str(avgtimeint)+' avg.freqstep=1', \
            log='$nameMS_initavg.log', commandType='DPPP')

MSs = lib_ms.AllMSs( glob.glob('mss-avg/TC*[0-9].MS'), s, check_flags=False )

fwhm = MSs.getListObj()[0].getFWHM(freq='mid')
detectability_dist = MSs.getListObj()[0].getFWHM(freq='max')*1.8/2. # 1.8 to go to close to the null
freq_min = np.min(MSs.getFreqs())
freq_mid = np.mean(MSs.getFreqs())
phase_center = MSs.getListObj()[0].getPhaseCentre()
timeint = MSs.getListObj()[0].getTimeInt()
ch_out = MSs.getChout(4e6) # for full band (48e6 MHz) is 12
ch_out_idg = 12 # better 24, but slow
#MSs.getListObj()[0].makeBeamReg('ddcal/beam.reg', freq='mid')
#beamReg = 'ddcal/beam.reg'

logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')
MSs.run('addcol2ms.py -m $pathMS -c FLAG_BKP -i FLAG', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
os.system('cp self/images/wideM-1* ddcal/init/')
full_image = lib_img.Image('ddcal/init/wideM-1-MFS-image.fits', userReg = userReg)

for cmaj in range(maxIter):
    logger.info('Starting major cycle: %i' % cmaj)
    
    # cycle specific variables
    picklefile = 'ddcal/directions-c%02i.pickle' % cmaj
    interp_h5parm = 'ddcal/c%02i/solutions/interp.h5' % cmaj
    aterm_config_file = 'ddcal/c%02i/aterm/aterm.config' % cmaj
    mask_cl = full_image.imagename.replace('.fits', '_mask-cl.fits')
    mask_ext = full_image.imagename.replace('.fits', '_mask-ext.fits')
    
    if not os.path.exists('ddcal/c%02i' % cmaj): os.makedirs('ddcal/c%02i' % cmaj)
    for subdir in ['plots','images','solutions','skymodels']:
        if not os.path.exists('ddcal/c%02i/%s' % (cmaj, subdir)): os.makedirs('ddcal/c%02i/%s' % (cmaj, subdir))

    if not os.path.exists(picklefile):
        directions = []

        ### group into patches corresponding to the mask islands
        if not os.path.exists(mask_cl): 
            full_image.makeMask(threshisl=7, atrous_do=False, remove_extended_cutoff=removeExtendedCutoff, only_beam=False, maskname=mask_cl, write_srl=True)
        if not os.path.exists(mask_ext) and cmaj == 0: 
            full_image.makeMask(threshisl=4, atrous_do=True, remove_extended_cutoff=0, only_beam=False, maskname=mask_ext)

        # the txt skymodel is used only to find directions
        if cmaj > 0:
            full_image.skymodel_cut = mask_cl.replace('fits','skymodel')
        elif not os.path.exists(full_image.skymodel_cut): 
            full_image.selectCC(checkBeam=False, maskname=mask_cl)

        # cleanup model
        if cmaj == 0:
            logger.info('Cleanup model...')
            for model_file in glob.glob(full_image.root+'*model.fits'):
                lib_img.blank_image_fits(model_file, mask_ext, model_file, inverse=True, blankval=0.)
        
        # locating DD-calibrators
        lsm = lsmtool.load(full_image.skymodel_cut)
        lsm.group(mask_cl, root='Isl')
        # This regroup nearby sources
        x = lsm.getColValues('RA',aggregate='wmean')
        y = lsm.getColValues('Dec',aggregate='wmean')
        flux = lsm.getColValues('I',aggregate='sum')
        grouper = lib_dd.Grouper(list(zip(x,y)), flux, look_distance=0.07, kernel_size=0.05, grouping_distance=0.03)
        grouper.run()
        clusters = grouper.grouping()
        grouper.plot()
        os.system('mv grouping*png ddcal/c%02i/plots/' % cmaj)
        patchNames = lsm.getPatchNames()
    
        logger.info('Merging nearby sources...')
        for cluster in clusters:
            patches = patchNames[cluster]
            if len(patches) > 1:
                lsm.merge(patches.tolist())
   
        lsm.setPatchPositions(method='mid')
        img_beam = full_image.getBeam()
        for name, size, ra, dec in \
                zip( lsm.getPatchNames(), lsm.getPatchSizes(units='deg'), \
                     lsm.getPatchPositions(asArray=True)[0], lsm.getPatchPositions(asArray=True)[1] ):

            # keep track of the spidx of sources
            idx = lsm.getRowIndex(name)
            fluxes = lsm.getColValues('I')[idx]
            spidx_coeffs = lsm.getColValues('SpectralIndex')[idx]
            ref_freq = lsm.getColValues('ReferenceFrequency')[idx]
            gauss_area = (lsm.getColValues('MajorAxis')[idx]*lsm.getColValues('MinorAxis')[idx])/(img_beam[0]*img_beam[1]) # in beams
            for i in range(len(idx)):
                if gauss_area[i] > 1:
                    fluxes[i] /= gauss_area[i] # reduce the fluxes for gaussians to the peak value

            d = lib_dd.Direction(name)
            d.set_flux(fluxes, spidx_coeffs, ref_freq )
            # skip faint directions
            if d.get_flux(freq_min) < min_cal_flux60*(freq_min/60e6)**(-0.8) or d.get_flux(freq_mid) < min_cal_flux60*(freq_mid/60e6)**(-0.8): continue
            if size < 2*img_beam[0]/3600:
                size = 4*img_beam[0]/3600

            #print('DEBUG:',name,fluxes,spidx_coeffs,gauss_area,ref_freq,size,img_beam,lsm.getColValues('MajorAxis')[idx])

            d.set_position( [ra, dec], distance_peeloff=detectability_dist, phase_center=phase_center )
            d.set_size(size*1.2) # size increased by 20%
            d.set_region(loc='ddcal/c%02i/skymodels' % cmaj)
            model_root = 'ddcal/c%02i/skymodels/%s-init' % (cmaj, name)
            for model_file in glob.glob(full_image.root+'*[0-9]-model.fits'):
                os.system('cp %s %s' % (model_file, model_file.replace(full_image.root, model_root)))
            d.set_model(model_root, typ='init', apply_region=True)

            directions.append(d)

        # order directions from the fluxiest
        directions = [x for _,x in sorted(zip([d.get_flux(freq_min) for d in directions],directions))][::-1] # reorder with flux

        for d in directions:
            if not d.peel_off:
                logger.info( '%s: min: %.2f Jy; mid: %.2f Jy' % (d.name, d.get_flux(freq_min), d.get_flux(freq_mid)) )
            else:
                logger.info( '%s: min: %.2f Jy; mid: %.2f Jy (peel off)' % (d.name, d.get_flux(freq_min), d.get_flux(freq_mid)) )

        # write file
        skymodel_cl = 'ddcal/c%02i/skymodels/cluster.skymodel' % cmaj
        lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
        lsm.setColValues('name', [x.split('_')[-1] for x in lsm.getColValues('patch')]) # just for the region - this makes this lsm useless
        lsm.write('ddcal/c%02i/skymodels/cluster.reg' % cmaj, format='ds9', clobber=True)
        del lsm
   
        pickle.dump( directions, open( picklefile, "wb" ) )
    else:
        directions = pickle.load( open( picklefile, "rb" ) )

    if cmaj == 0:
        if w.todo('c%02i-fullpredict' % cmaj):
            # wsclean predict
            logger.info('Predict full model...')
            if cmaj == 0:
                s.add('wsclean -predict -name '+full_image.root+' -j '+str(s.max_processors)+' -channels-out '+str(ch_out)+' \
                        -reorder -parallel-reordering 4 '+MSs.getStrWsclean(), \
                        log='wscleanPRE-c'+str(cmaj)+'.log', commandType='wsclean', processors='max')
                s.run(check=True)
    
            w.done('c%02i-fullpredict' % cmaj)
        ### DONE
    
        if w.todo('c%02i-fullsub' % cmaj):
            # subtract - ms:SUBTRACTED_DATA = DATA - MODEL_DATA
            logger.info('Set SUBTRACTED_DATA = DATA - MODEL_DATA...')
            MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"', \
                        log='$nameMS_taql.log', commandType='general')
            # reset - ms:CORRECTED_DATA = DATA
            logger.info('Set CORRECTED_DATA = DATA...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = DATA"', \
                        log='$nameMS_taql.log', commandType='general')
    
            w.done('c%02i-fullsub' % cmaj)
        ### DONE

        ### TESTTESTTEST: empty image
        #if not os.path.exists('img/empty-init-c'+str(cmaj)+'-image.fits'):
        #    clean('init-c'+str(cmaj), MSs, size=(fwhm*1.5,fwhm*1.5), res='normal', empty=True)
        ###

    for dnum, d in enumerate(directions):

        logger.info('c%02i - Working on direction: %s (%f Jy - %f deg)' % (cmaj, d.name, d.get_flux(freq_mid), d.size))
        if d.size > 0.5: logger.warning('Patch size large: %f' % d.size)
        logstring = 'c%02i-%s' % (cmaj,d.name)

        if w.todo('%s-predict' % logstring):

            if cmaj == 0:
                # Predict - ms:MODEL_DATA
                logger.info('Predict model...')
                s.add('wsclean -predict -name '+d.get_model('init')+' -j '+str(s.max_processors)+' -channels-out '+str(ch_out)+' '+MSs.getStrWsclean(), \
                        log='wscleanPRE-'+logstring+'.log', commandType='wsclean', processors='max')
                s.run(check=True)
    
                # Add back the model previously subtracted for this dd-cal
                logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA + MODEL_DATA...')
                MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', \
                        log='$nameMS_taql.log', commandType='general')
    
            else:

                # these dd-cal are not in the data anymore
                if d.peel_off:
                    logger.info('This sources has been peeled, skip.')
                    continue

                # DDF predict+corrupt in MODEL_DATA of everything BUT the calibrator
                indico = full_image.root+'.DicoModel'
                outdico = indico+'-'+d.name
                outmask = outdico+'.mask'
                lib_img.blank_image_reg(mask_ext, d.get_region(), outfile=outmask, inverse=False, blankval = 0.)
                s.add('MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s' % (outmask, indico, outdico), \
                       log='MaskDicoModel-'+logstring+'.log', commandType='python', processors='max')
                s.run(check=True)

                ddf_parms = {
                    'Data_MS':MSs.getStrDDF(),
                    'Data_ColName':'CORRECTED_DATA',
                    'Data_Sort':1,
                    'Output_Mode':'Predict',
                    'Predict_InitDicoModel':outdico,
                    'Predict_ColName':'MODEL_DATA',
                    'Deconv_Mode':'HMP',
                    #'Deconv_CycleFactor':0,
                    #'Deconv_MaxMinorIter':1000000,
                    #'Deconv_RMSFactor':3.0,
                    #'Deconv_PeakFactor':0.005,
                    #'Deconv_FluxThreshold':0.0,
                    #'Weight_Robust':-0.5,
                    'Image_NPix':8750,
                    'CF_wmax':50000,
                    'CF_Nw':100,
                    'Beam_CenterNorm':1,
                    'Beam_Smooth':1,
                    'Beam_Model':'LOFAR',
                    'Beam_LOFARBeamMode':'A',
                    'Beam_NBand':1,
                    'Beam_DtBeamMin':5,
                    'Output_Also':'onNeds',
                    'Image_Cell':3.,
                    'Freq_NDegridBand':ch_out,
                    'Freq_NBand':ch_out,
                    #'Mask_Auto':1,
                    #'Mask_SigTh':5.0,
                    'GAClean_MinSizeInit':10,
                    'Facets_DiamMax':1.5,
                    'Facets_DiamMin':0.1,
                    'Weight_ColName':'WEIGHT_SPECTRUM',
                    #'Output_Name':imagename,
                    'Comp_BDAMode':1,
                    'DDESolutions_DDModeGrid':'AP',
                    'DDESolutions_DDModeDeGrid':'AP',
                    'RIME_ForwardMode':'BDA-degrid',
                    #'Output_RestoringBeam':15.,
                    'DDESolutions_DDSols':interp_h5parm.replace('c%02i' % cmaj, 'c%02i' % (cmaj-1))+':sol000/phase000+amplitude000'
                    }
                logger.info('Predict corrupted rest-of-the-sky...')
                lib_util.run_DDF(s, 'ddfacet-pre-'+logstring+'.log', **ddf_parms, Cache_Reset=1)

                #                run("DDF.py --Output-Name=image_full_ampphase_di_m.NS_SUB --Data-ChunkHours=" + str(args['chunkhours']) + " --Data-MS=" + args['mslist'] + " --Deconv-PeakFactor 0.001000 --Data-ColName " + args['column'] + " --Parallel-NCPU="+str(ncpu) + " --Facets-CatNodes=" + clustercat + " --Beam-CenterNorm=1 --Deconv-Mode SSD --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust " + str(robust) +" --Image-NPix=" + str(imagenpix) + " --CF-wmax 50000 --CF-Nw 100 --Output-Also onNeds --Image-Cell "+ str(imagecell) + " --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=3.000000 --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=. --Log-Memory 1 --Cache-Weight=reset --Output-Mode=Predict --Output-RestoringBeam 6.000000 --Freq-NBand=2 --RIME-DecorrMode=FT --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0 --Mask-Auto=1 --Mask-SigTh=5.00 --Mask-External=" + outmask + " --DDESolutions-GlobalNorm=None --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --DDESolutions-DDSols=[" + ddsolstr + "] --Predict-InitDicoModel=" + outdico + " --Selection-UVRangeKm=" + uvsel + " --GAClean-MinSizeInit=10 --Cache-Reset 1 --Beam-Smooth=1 --Predict-ColName='PREDICT_SUB' --DDESolutions-SolsDir=SOLSDIR")

                # Remove corrupted data from CORRECTED_DATA
                logger.info('Set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA...')
                MSs.run('taql "update $pathMS set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"', \
                        log='$nameMS_taql.log', commandType='general')

                ### TTESTTESTTEST: empty image
                if not os.path.exists('img/almostempty-%02i-%s-image.fits' % (dnum, logstring)):
                    clean('%02i-%s' % (dnum, logstring), MSs, size=(fwhm*1.5,fwhm*1.5), res='normal', empty=True)
    
            w.done('%s-predict' % logstring)
 
        ### DONE

        if w.todo('%s-shift' % logstring):
            logger.info('Phase shift and avg...')

            lib_util.check_rm('mss-dir')
            os.makedirs('mss-dir')

            # Shift - ms:SUBTRACTED_DATA -> ms:DATA
            if d.get_flux(freq_mid) > 4: avgtimeint = int(round(15/timeint))
            else: avgtimeint = int(round(30/timeint))
            MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=SUBTRACTED_DATA msout.datacolumn=DATA \
                    avg.timestep='+str(avgtimeint)+' avg.freqstep=8 shift.phasecenter=\['+str(d.position[0])+'deg,'+str(d.position[1])+'deg\]', \
                    log='$nameMS_shift-'+logstring+'.log', commandType='DPPP')

            w.done('%s-shift' % logstring)
        ### DONE

        MSs_dir = lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s, check_flags=False )

        if w.todo('%s-flag' % logstring):
            logger.info('Flag on mindata...')
            MSs_dir.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

            w.done('%s-flag' % logstring)
        ### DONE

        # Correct for beam in that direction
        if not d.peel_off:
            if w.todo('%s-beamcorr' % logstring):
                logger.info('Correcting beam...')
                # Convince DPPP that DATA is corrected for the beam in the phase centre
                MSs_dir.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=DATA \
                        setbeam.direction=\['+str(phase_center[0])+'deg,'+str(phase_center[1])+'deg\] \
                        corrbeam.direction=\['+str(d.position[0])+'deg,'+str(d.position[1])+'deg\] corrbeam.invert=True', \
                        log='$nameMS_beam-'+logstring+'.log', commandType='DPPP')
    
                w.done('%s-beamcorr' % logstring)
            ### DONE

        if w.todo('%s-preimage' % logstring):
            logger.info('Pre-imaging...')
            clean('%s-pre' % logstring, MSs_dir, res='normal', size=[d.size,d.size], imagereg=d.get_region())

            w.done('%s-preimage' % logstring)
        ### DONE
        
        # get initial noise and set iterators for timeint solutions
        image = lib_img.Image('img/ddcalM-%s-pre-MFS-image.fits' % logstring)
        rms_noise_pre = image.getNoise(); rms_noise_init = rms_noise_pre
        mm_ratio_pre = image.getMaxMinRatio(); mm_ratio_init = mm_ratio_pre
        doamp = False
        # usually there are 3600/30=120 or 3600/15=240 timesteps, try to use multiple numbers
        iter_ph_solint = lib_util.Sol_iterator([4,1])
        iter_amp_solint = lib_util.Sol_iterator([120,60,30])
        iter_amp2_solint = lib_util.Sol_iterator([120,60])
        logger.info('RMS noise (init): %f' % (rms_noise_pre))
        logger.info('MM ratio (init): %f' % (mm_ratio_pre))

        for cdd in range(10):

            logger.info('c%02i - %s: Starting dd cycle: %02i' % (cmaj, d.name, cdd))
            logstringcal = logstring+'-cdd%02i' % cdd

            ################################################################
            # Calibrate
            solint_ph = next(iter_ph_solint)
            d.add_h5parm('ph', 'ddcal/c%02i/solutions/cal-ph-%s.h5' % (cmaj,logstringcal) )
            if doamp:
                solint_amp = next(iter_amp_solint)
                d.add_h5parm('amp1', 'ddcal/c%02i/solutions/cal-amp1-%s.h5' % (cmaj,logstringcal) )
                solint_amp2 = next(iter_amp2_solint)
                d.add_h5parm('amp2', 'ddcal/c%02i/solutions/cal-amp2-%s.h5' % (cmaj,logstringcal) )
   
            if w.todo('%s-calibrate' % logstringcal):

                if cdd == 0:
                    logger.info('BL-based smoothing...')
                    # Smoothing - ms:DATA -> ms:SMOOTHED_DATA
                    MSs_dir.run('BLsmooth.py -r -i DATA -o SMOOTHED_DATA $pathMS', \
                        log='$nameMS_smooth-'+logstringcal+'.log', commandType='python')    
 
                # Calibration - ms:SMOOTHED_DATA
                # possible to put nchan=6 if less channels are needed in the h5parm (e.g. for IDG)
                logger.info('Gain phase calibration...')
                MSs_dir.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/cal-ph.h5 \
                    sol.mode=scalarcomplexgain sol.solint='+str(solint_ph)+' sol.nchan=1 sol.smoothnessconstraint=5e6 \
                    sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]]', \
                    log='$nameMS_solGph-'+logstringcal+'.log', commandType='DPPP')
                lib_util.run_losoto(s, 'ph', [ms+'/cal-ph.h5' for ms in MSs_dir.getListStr()], \
                    [parset_dir+'/losoto-plot1.parset'], plots_dir='ddcal/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                os.system('mv cal-ph.h5 %s' % d.get_h5parm('ph'))

                # correct ph - ms:DATA -> ms:CORRECTED_DATA
                logger.info('Correct ph...')
                MSs_dir.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=CORRECTED_DATA \
                             cor.parmdb='+d.get_h5parm('ph')+' cor.correction=phase000', \
                             log='$nameMS_correct-'+logstringcal+'.log', commandType='DPPP')
                #MSs_dir.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                #             cor.parmdb=ddcal/solutions/cal-ph-'+logstringcal+'.h5 cor.correction=tec000', \
                #             log='$nameMS_correct-'+logstringcal+'.log', commandType='DPPP')

                if doamp:
                    
                    logger.info('Gain amp calibration 1...')
                    # Calibration - ms:CORRECTED_DATA
                    # possible to put nchan=6 if less channels are needed in the h5parm (e.g. for IDG)
                    MSs_dir.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp1.h5 \
                        sol.mode=diagonal sol.solint='+str(solint_amp)+' sol.nchan=1 sol.uvmmin=200 sol.smoothnessconstraint=4e6 \
                        sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]]', \
                        log='$nameMS_solGamp1-'+logstringcal+'.log', commandType='DPPP')
                    lib_util.run_losoto(s, 'amp1', [ms+'/cal-amp1.h5' for ms in MSs_dir.getListStr()], \
                        [parset_dir+'/losoto-clip.parset', parset_dir+'/losoto-norm.parset', parset_dir+'/losoto-plot2.parset'], \
                        plots_dir='ddcal/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                    os.system('mv cal-amp1.h5 %s' % d.get_h5parm('amp1'))

                    logger.info('Correct amp 1...')
                    # correct amp - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                    MSs_dir.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                        cor.parmdb='+d.get_h5parm('amp1')+' cor.correction=amplitude000', \
                        log='$nameMS_correct-'+logstringcal+'.log', commandType='DPPP') 

                    logger.info('Gain amp calibration 2...')
                    # Calibration - ms:SMOOTHED_DATA
                    MSs_dir.run('DPPP '+parset_dir+'/DPPP-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/cal-amp2.h5 \
                        sol.mode=diagonal sol.solint='+str(solint_amp2)+' sol.nchan=6 sol.uvmmin=200 sol.smoothnessconstraint=10e6', \
                        log='$nameMS_solGamp2-'+logstringcal+'.log', commandType='DPPP')
                    lib_util.run_losoto(s, 'amp2', [ms+'/cal-amp2.h5' for ms in MSs_dir.getListStr()], \
                        [parset_dir+'/losoto-clip2.parset', parset_dir+'/losoto-norm.parset', parset_dir+'/losoto-plot3.parset'], \
                        plots_dir='ddcal/c%02i/plots/plots-%s' % (cmaj,logstringcal))
                    os.system('mv cal-amp2.h5 %s' % d.get_h5parm('amp2'))

                    logger.info('Correct amp 2...')
                    # correct amp2 - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
                    MSs_dir.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA \
                        cor.parmdb='+d.get_h5parm('amp2')+' cor.correction=amplitude000', \
                        log='$nameMS_correct-'+logstringcal+'.log', commandType='DPPP') 

                w.done('%s-calibrate' % logstringcal)
            ### DONE

            ###########################################################################
            # Imaging
            if w.todo('%s-image' % logstringcal):

                logger.info('%s (cdd: %02i): imaging...' % (d.name, cdd))
                clean('%s' % logstringcal, MSs_dir, res='normal', size=[d.size,d.size], imagereg=d.get_region())

                w.done('%s-image' % logstringcal)
            ### DONE

            image = lib_img.Image('img/ddcalM-%s-MFS-image.fits' % logstringcal, userReg=userReg)
            # something went wrong during last imaging, break
            if not os.path.exists(image.imagename):
                break
            d.set_model(image.root, typ='best', apply_region=False) # currently best model
            # get noise, if larger than prev cycle: break
            rms_noise = image.getNoise()
            mm_ratio = image.getMaxMinRatio()
            logger.info('RMS noise (cdd:%02i): %f' % (cdd,rms_noise))
            logger.info('MM ratio (cdd:%02i): %f' % (cdd,mm_ratio))
            if rms_noise > 0.99*rms_noise_pre and mm_ratio < 1.01*mm_ratio_pre:
                if   mm_ratio < 10 and cdd >= 2: break
                elif mm_ratio < 20 and cdd >= 3: break
                elif mm_ratio < 30 and cdd >= 4: break
                elif cdd >= 5: break

            if cdd >= 4 and mm_ratio >= 30:
                doamp = True

            rms_noise_pre = rms_noise
            mm_ratio_pre = mm_ratio

        # End calibration cycle
        ##################################

        # if divergency or died the first cycle, don't subtract
        if cdd == 0:
            d.converged = False
            logger.warning('%s: something went wring during the first self-cal cycle in this direction.' % (d.name))
            d.clean()
            continue
        elif rms_noise_pre > rms_noise_init:
            d.converged = False
            logger.warning('%s: noise did not decresed (%f -> %f), do not subtract source.' % (d.name, rms_noise_init, rms_noise_pre))
            d.clean()
            continue
        else:
            d.converged = True
            # copy in the ddcal dir the best model
            model_skymodel = 'ddcal/c%02i/skymodels/%s-best-source.txt' % (cmaj, d.name)
            model_skydb = 'ddcal/c%02i/skymodels/%s-best-source.skydb' % (cmaj, d.name)
            os.system('cp %s %s' % (d.get_model('best')+'-sources.txt', model_skymodel))
            s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (model_skymodel, model_skydb), log='makesourcedb_cl.log', commandType='general' )
            s.run()
        
        # remove the DD-cal from original dataset using new solutions
        if w.todo('%s-subtract' % logstring):

            # Predict - ms:MODEL_DATA
            logger.info('Add best model to MODEL_DATA...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+model_skydb, \
                    log='$nameMS_pre-'+logstring+'.log', commandType='DPPP')

            # Store FLAGS
            MSs.run('taql "update $pathMS set FLAG_BKP = FLAG"', \
                    log='$nameMS_taql.log', commandType='general')

            # Corrput now model - ms:MODEL_DATA -> MODEL_DATA
            logger.info('Corrupt ph...')
            MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.invert=False cor.parmdb='+d.get_h5parm('ph',-2)+' cor.correction=phase000', \
                        log='$nameMS_corruping'+logstring+'.log', commandType='DPPP')
            #MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
            #            cor.invert=False cor.parmdb='+d.get_h5parm('ph',-2)+' cor.correction=tec000', \
            #            log='$nameMS_corrupt-'+logstring+'.log', commandType='DPPP')

            if not d.get_h5parm('amp1',-2) is None:
                logger.info('Corrupt amp...')
                MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                       cor.invert=False cor.parmdb='+d.get_h5parm('amp1',-2)+' cor.correction=amplitude000', \
                       log='$nameMS_corrupt-'+logstring+'.log', commandType='DPPP') 
            if not d.get_h5parm('amp2',-2) is None:
                MSs.run('DPPP '+parset_dir+'/DPPP-correct.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                       cor.invert=False cor.parmdb='+d.get_h5parm('amp2',-2)+' cor.correction=amplitude000', \
                       log='$nameMS_corrupt-'+logstring+'.log', commandType='DPPP') 

            if not d.peel_off:
                # Corrupt for the beam
                logger.info('Corrupting beam...')
                # Convince DPPP that MODELDATA is corrected for the beam in the dd-cal direction, so I can corrupt
                MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                       setbeam.direction=\['+str(d.position[0])+'deg,'+str(d.position[1])+'deg\] \
                       corrbeam.direction=\['+str(d.position[0])+'deg,'+str(d.position[1])+'deg\] corrbeam.invert=False', \
                       log='$nameMS_beam-'+logstring+'.log', commandType='DPPP')
                #[MS.delBeamInfo(col='MODEL_DATA') for MS in MSs.getListObj()]
                MSs.run('DPPP '+parset_dir+'/DPPP-beam2.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA steps=corrbeam \
                       corrbeam.direction=\['+str(phase_center[0])+'deg,'+str(phase_center[1])+'deg\] corrbeam.beammode=element corrbeam.invert=True', \
                       log='$nameMS_beam-'+logstring+'.log', commandType='DPPP')

            # Set MODEL_DATA = 0 where data are flagged, then unflag everything
            MSs.run('taql "update $pathMS set MODEL_DATA[FLAG] = 0"', \
                    log='$nameMS_taql.log', commandType='general')

            # Restore of FLAGS
            MSs.run('taql "update $pathMS set FLAG = FLAG_BKP"', \
                    log='$nameMS_taql.log', commandType='general')

            # Remove the ddcal again
            logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA')
            MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', \
                    log='$nameMS_taql.log', commandType='general')

            # if it's a source to peel, remove it from the data column used for imaging
            if d.peel_off:
                logger.info('Source to peel: set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA')
                MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', \
                        log='$nameMS_taql.log', commandType='general')

            w.done('%s-subtract' % logstring)
        ### DONE

        ### TTESTTESTTEST: empty image
        if not os.path.exists('img/empty-%02i-%s-image.fits' % (dnum, logstring)):
            clean('%02i-%s' % (dnum, logstring), MSs, size=(fwhm*1.5,fwhm*1.5), res='normal', empty=True)
        ###

    ######################################################
    # full imaging
    
    imagename = 'img/wideDD-c%02i' % (cmaj)

    if w.todo('c%02i-imaging' % cmaj):

        # combine the h5parms
        h5parms = {'ph':[], 'amp1':[], 'amp2':[]}
        for d in directions:
            # only those who converged
            if d.peel_off:
                continue
            if not d.converged:
                continue
    
            h5parms['ph'].append(d.get_h5parm('ph',-2))
            if d.get_h5parm('amp1',-2) is not None:
                h5parms['amp1'].append(d.get_h5parm('amp1',-2))
            if d.get_h5parm('amp2',-2) is not None:
                h5parms['amp2'].append(d.get_h5parm('amp2',-2))
            
            log = '%s: Phase (%s)' % (d.name, d.get_h5parm('ph',-2))
            log += ' Amp1 (%s)' % (d.get_h5parm('amp1',-2))
            log += ' Amp2 (%s)' % (d.get_h5parm('amp2',-2))
            logger.info(log)
    
    
        for typ, h5parm_list in h5parms.items():
            # rename direction
            for h5parmFile in h5parm_list:
                dirname = h5parmFile.split('-')[3]
                lib_h5.repoint(h5parmFile, dirname)
    
                if typ == 'ph':
                    lib_h5.addpol(h5parmFile, 'phase000')
                    lib_h5.addpol(h5parmFile, 'amplitude000')
                    # reset high-res amplitudes in ph-solve
                    s.add('losoto -v '+h5parmFile+' '+parset_dir+'/losoto-resetamp.parset ', log='h5parm_collector.log', commandType='python' )
                    s.run()
                    s.add('losoto -v '+h5parmFile+' '+parset_dir+'/losoto-refph.parset ', log='h5parm_collector.log', commandType='python' )
                    s.run()
    
                if typ == 'amp1' or typ == 'amp2':
                    s.add('losoto -v '+h5parmFile+' '+parset_dir+'/losoto-resetph.parset ', log='h5parm_collector.log', commandType='python' )
                    s.run()
    
        lib_util.check_rm(interp_h5parm)
        logger.info('Interpolating solutions...')
        s.add('H5parm_interpolator.py -o '+interp_h5parm+' '+' '.join(h5parms['ph']+h5parms['amp1']+h5parms['amp2']), log='h5parm_interpolator.log', commandType='python' )
        s.run()
    
        idg = False
        if idg:
    
            for i, MS in enumerate(MSs.getStrObj()):
                s.add('/home/fdg/scripts/LiLF/scripts/make_gain_screen.py \
                        -m '+MS.path+' -p '+interp_h5parm+' -o ddcal/c%02i/aterm/TC%02i' % (cmaj, i), \
                        log='h5parm_interpolator.log', commandType='python')
                s.run()
     
            # create aterm config file
            with open(aterm_config_file, 'w') as file:  # Use file to refer to the file object
                file.write('aterms = [diagonal, beam]')
                file.write('diagonal.images = ['+' '.join(sorted(glob.glob('ddcal/c%02i/aterm/TC*fits' % cmaj)))+']')
                file.write('diagonal.window = tukey\n diagonal.update_interval  = 48.066724')
                file.write('beam.differential = true\n beam.update_interval = 120\n beam.usechannelfreq = true')
        
            # run the imager
            lib_util.run_wsclean(s, 'wsclean-c'+str(cmaj)+'.log', MSs.getStrWsclean(), name=imagename, size='6000 6000', save_source_list='', scale='5arcsec', \
                        weight='briggs -0.3', niter=2000, no_update_model_required='', minuv_l=30, mgain=0.85, \
                        multiscale='', multiscale_scale_bias=0.65, multiscale_scales='0,10,20,40,80',
                        parallel_deconvolution=512, local_rms='', auto_threshold=0.5, auto_mask=1.5, \
                        join_channels='', fit_spectral_pol=3, channels_out=ch_out_idg, deconvolution_channels=3, \
                        temp_dir='./', pol='I', use_idg='', aterm_config=aterm_config_file, aterm_kernel_size=45, nmiter=4 )

            sys.exit('Not implemente further on...')
    
        else:
    
            ddf_parms = {
                    'Data_MS':MSs.getStrDDF(),
                    'Data_ColName':'CORRECTED_DATA',
                    'Data_Sort':1,
                    'Output_Mode':'Clean',
                    'Deconv_CycleFactor':0,
                    'Deconv_MaxMinorIter':1000000,
                    'Deconv_Mode':'HMP',
                    'Deconv_RMSFactor':3.0,
                    'Deconv_PeakFactor':0.005,
                    'Deconv_FluxThreshold':0.0,
                    'Weight_Robust':-0.5,
                    'Image_NPix':8750,
                    'CF_wmax':50000,
                    'CF_Nw':100,
                    'Beam_CenterNorm':1,
                    'Beam_Smooth':1,
                    'Beam_Model':'LOFAR',
                    'Beam_LOFARBeamMode':'A',
                    'Beam_NBand':1,
                    'Beam_DtBeamMin':5,
                    'Output_Also':'onNeds',
                    'Image_Cell':3.,
                    #'Facets_NFacets':25,
                    'Freq_NDegridBand':ch_out,
                    'Freq_NBand':ch_out,
                    'Mask_Auto':1,
                    'Mask_SigTh':5.0,
                    'GAClean_MinSizeInit':10,
                    'Facets_DiamMax':1.5,
                    'Facets_DiamMin':0.1,
                    'Weight_ColName':'WEIGHT_SPECTRUM',
                    'Output_Name':imagename,
                    'Comp_BDAMode':1,
                    'DDESolutions_DDModeGrid':'AP',
                    'DDESolutions_DDModeDeGrid':'AP',
                    'RIME_ForwardMode':'BDA-degrid',
                    'Output_RestoringBeam':15.,
                    'DDESolutions_DDSols':interp_h5parm+':sol000/phase000+amplitude000'
                    }
            
            logger.info('Cleaning 1...')
            lib_util.run_DDF(s, 'ddfacet-c'+str(cmaj)+'.log', **ddf_parms,
                    Deconv_MaxMajorIter=1,
                    Cache_Reset=1
                    )
    
            # make mask
            im = lib_img.Image(imagename+'.app.restored.fits', userReg=userReg)
            im.makeMask(threshisl = 3, rmsbox=(70,5))
    
            logger.info('Cleaning 2...')
            lib_util.run_DDF(s, 'ddfacet-c'+str(cmaj)+'.log', **ddf_parms,
                    Cache_Reset=0,
                    Cache_Dirty='forcedirty',
                    Cache_PSF='force',
                    Cache_SmoothBeam='force',
                    Deconv_MaxMajorIter=3,
                    Predict_InitDicoModel=imagename+'.DicoModel',
                    Mask_External=im.maskname
                    )
 
            os.system('mv %s* ddcal/c%02i/images' % (imagename, cmaj))
        w.done('c%02i-imaging' % cmaj)
    ### DONE

    full_image = lib_img.Image('ddcal/c%02i/images/%s.app.restored.fits' % (cmaj,imagename.split('/')[-1]), userReg = userReg)
    min_cal_flux60 *= 0.8 # go a bit deeper

if w.todo('upload'):
    
    logger.info('Save final images...')
    targetname = os.getcwd().split('/')[-1]
    logger.info('Copy: ddcal/c0*/images/img/wideDD-c*... -> lofar.herts.ac.uk:/beegfs/lofar/lba/products/%s' % targetname)
    os.system('ssh lofar.herts.ac.uk "rm -rf /beegfs/lofar/lba/products/%s"' % targetname)
    os.system('ssh lofar.herts.ac.uk "mkdir /beegfs/lofar/lba/products/%s"' % targetname)
    os.system('scp -q ddcal/c0*/images/img/wideDD-c*.app.restored.fits lofar.herts.ac.uk:/beegfs/lofar/lba/products/%s' % targetname)
    os.system('scp -q ddcal/c0*/images/img/wideDD-c*.int.restored.fits lofar.herts.ac.uk:/beegfs/lofar/lba/products/%s' % targetname)
    
    w.done('upload')
### DONE

logger.info("Done.")
