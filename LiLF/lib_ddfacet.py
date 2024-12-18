import os, sys, inspect
from functools import wraps
import numpy as np
from astropy.io import fits

from LiLF.lib_log import logger
from LiLF import lib_img


def killms_data(s, logfile, MSs, imagename, outsols, clusterfile=None, colname='DATA', niterkf=6, dicomodel=None,
                uvrange=None, wtuv=None, robust=None, dt=None, cache_dir='ddfcal/cache', sols_dir='ddfcal/solutions',
                SolverType="KAFCA", PolMode="Scalar", MergeSmooth=False, NChanSols=1, 
                DISettings=None, EvolutionSolFile=None, CovQ=0.1, InterpToMSListFreqs=None, 
                SkipSmooth=False, PreApplySols=None, SigmaFilterOutliers=None):

    if not os.path.isdir(cache_dir):
        raise RuntimeError('Missing cache dir.')
    if not os.path.isdir(sols_dir):
        raise RuntimeError('Missing sols dir.')

    # run killms individually on each MS -- allows restart if it failed in the middle
    for f in MSs.getListStr():

        MSName  = os.path.abspath(f).split("/")[-1]
        solname = os.path.abspath(sols_dir)+"/"+MSName+'/killMS.'+outsols+'.sols.npz'
        checkname=solname

        #checkname=f+'/killMS.'+outsols+'.sols.npz'
        if os.path.isfile(checkname):

            logger.warning('Solutions file '+checkname+' already exists, not running killMS step.')
            
        else:
            runcommand = "kMS.py --MSName %s --SolverType %s --PolMode %s --BaseImageName %s --dt %f --NIterKF %i --CovQ %f --LambdaKF 0.5 --NCPU %i --OutSolsName %s --InCol %s --DoBar 1 --SolsDir %s" \
                    % (f, SolverType, PolMode, imagename, dt, niterkf, CovQ, s.max_processors, outsols, colname, sols_dir)

            # set weights
            if robust is None:
                runcommand += ' --Weighting Natural'
            else:
                runcommand += ' --Weighting Briggs --Robust=%f' % robust

            # set uv-range
            if uvrange is not None:
                if wtuv is not None:
                    runcommand+=' --WTUV=%f --WeightUVMinMax=%f,%f' % (wtuv, uvrange[0], uvrange[1])
                else:
                    runcommand+=' --UVMinMax=%f,%f' % (uvrange[0], uvrange[1])

            if PreApplySols:
                runcommand+=' --PreApplySols=[%s]'%PreApplySols

            if DISettings is None:
                runcommand+=' --NChanSols %i' % NChanSols
                runcommand+=' --BeamMode LOFAR --LOFARBeamMode=A --DDFCacheDir=%s'%cache_dir
                if clusterfile is not None:
                    runcommand+=' --NodesFile '+clusterfile
                if dicomodel is not None:
                    runcommand+=' --DicoModel '+dicomodel
                if EvolutionSolFile is not None:
                    runcommand+=' --EvolutionSolFile '+EvolutionSolFile
                    
            else:
                runcommand+=" --SolverType %s --PolMode %s --SkyModelCol %s --OutCol %s --ApplyToDir 0"%DISettings
                _,_,ModelColName,_=DISettings
                _,dt,_,n_df=give_dt_dnu(f,
                                        DataCol=colname,
                                        ModelCol=ModelColName,
                                        T=10.)
                runcommand+=" --dt %f --NChanSols %i"%(dt+1e-4,n_df)
                
                
            s.add(runcommand, log=logfile, commandType='singularity')
            s.run(check=True)

            # Clip anyway - on IMAGING_WEIGHT by default
            if DISettings is not None:
                ClipCol=DISettings[-1]
            else:
                ClipCol=colname

            runcommand="ClipCal.py --MSName %s --ColName %s"%(f,ClipCol)
            s.add(runcommand, log=logfile, commandType='singularity')
            s.run(check=True)

    if MergeSmooth:
        print('ADD SMOOTHSOL')
        #outsols=smooth_solutions(mslist,outsols,dryrun=o['dryrun'],InterpToMSListFreqs=InterpToMSListFreqs,
        #                         SkipSmooth=SkipSmooth,SigmaFilterOutliers=SigmaFilterOutliers)
        
    return outsols


def ddf_image(s, logfile, MSs, imagename, cleanmask=None, cleanmode='HMP', ddsols=None, applysols=None, threshold=None, majorcycles=3, use_dicomodel=False, robust=0, beamsize=None, beamsize_minor=None, beamsize_pa=None, reuse_psf=False, reuse_dirty=False, verbose=False, saveimages=None, imsize=None, cellsize=None, uvrange=None, colname='CORRECTED_DATA', peakfactor=0.1, dicomodel_base=None, do_decorr=None, normalization=None, dirty_from_resid=False, clusterfile=None, HMPsize=10, automask=True, automask_threshold=10.0, smooth=False, noweights=False, cubemode=False, apply_weights=True, use_weightspectrum=False, rms_factor=3.0, predict_column=None, conditional_clearcache=False, PredictSettings=None, RMSFactorInitHMP=1., MaxMinorIterInitHMP=10000, OuterSpaceTh=None, AllowNegativeInitHMP=False, phasecenter=None, polcubemode=False, channels=None, startchan=None, endchan=None, stokes=None):

    # prepare MSsfile
    mslist = ','.join(MSs.getListStr())
    #with open(mslist, 'w') as the_file:
    #    for MS in MSs.getListObj():
    #        the_file.write(MS.pathMS+'\n')

    # saveimages lists _additional_ images to save
    if saveimages is None:
        saveimages=''
    saveimages+='onNeds'
    
    #if options is None:
    #    options=o # attempt to get global if it exists

    if do_decorr is None:
        do_decorr=True
    if beamsize is None:
        beamsize=30 # arcsec
    if imsize is None:
        imsize=8750
    if cellsize is None:
        cellsize=3
        
    cache_dir = 'ddfcal/cache'
    if not os.path.isdir(cache_dir):
        raise RuntimeError('Missing cache dir.')

    if majorcycles>0:
        fname=imagename+'.app.restored.fits'
    else:
        fname=imagename+'.dirty.fits'

    if PredictSettings is not None and PredictSettings[0]=="Predict":
        fname="_has_predicted_OK.%s.info"%imagename

    ncpu = s.maxThreads

    runcommand = "DDF.py --Output-Name=%s --Data-MS=%s --Deconv-PeakFactor %f --Data-ColName %s --Parallel-NCPU=%i --Beam-CenterNorm=1 --Deconv-CycleFactor=0 --Deconv-MaxMinorIter=1000000 --Deconv-MaxMajorIter=%s --Deconv-Mode %s --Beam-Model=LOFAR --Beam-LOFARBeamMode=A --Weight-Robust %f --Image-NPix=%i --CF-wmax 50000 --CF-Nw 100 --Output-Also %s --Image-Cell %f --Facets-NFacets=11 --SSDClean-NEnlargeData 0 --Freq-NDegridBand 1 --Beam-NBand 1 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Deconv-RMSFactor=%f --SSDClean-ConvFFTSwitch 10000 --Data-Sort 1 --Cache-Dir=%s --Log-Memory 1"%(imagename,mslist,peakfactor,colname,ncpu,majorcycles,cleanmode,robust,imsize,saveimages,float(cellsize),rms_factor,cache_dir)

    runcommand += " --GAClean-RMSFactorInitHMP %f"%RMSFactorInitHMP
    runcommand += " --GAClean-MaxMinorIterInitHMP %f"%MaxMinorIterInitHMP
    if AllowNegativeInitHMP:
        runcommand += " --GAClean-AllowNegativeInitHMP True"
    if OuterSpaceTh is not None:
        runcommand += " --HMP-OuterSpaceTh %f"%OuterSpaceTh
        
    runcommand+=' --DDESolutions-SolsDir=ddfcal/solutions'
    runcommand+=' --Cache-Weight=reset'

    
    if PredictSettings is None:
        runcommand += " --Output-Mode=Clean"
    else:
        if len(PredictSettings) == 2:
            runcommand += " --Output-Mode=%s --Predict-ColName %s"%PredictSettings
        elif len(PredictSettings) == 3:
            runcommand += " --Output-Mode=%s --Predict-ColName %s --Predict-MaskSquare [0,%i]"%PredictSettings
        else:
            raise RuntimeError('PredictSettings has the wrong dimensions %s '%PredictSettings)

    if beamsize_minor is not None:
        runcommand += ' --Output-RestoringBeam %f,%f,%f'%(beamsize,beamsize_minor,beamsize_pa)
    elif beamsize is not None:
        runcommand += ' --Output-RestoringBeam %f'%(beamsize)
    
    if apply_weights:
        runcommand+=' --Weight-ColName="IMAGING_WEIGHT"'
    else:
        if not use_weightspectrum:
            runcommand+=' --Weight-ColName="None"'
        else:
            runcommand+=' --Weight-ColName="WEIGHT_SPECTRUM"'

#    if cubemode:
#        # number of channels equals number of distinct freqs in data
#        freqs=[]
#        mss=[l.rstrip() for l in open(mslist).readlines()]
#        for ms in mss:
#            t = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
#            freq=t[0]['REF_FREQUENCY']
#            if freq not in freqs:
#                freqs.append(freq)
#        channels=len(freqs)
#        runcommand+=' --Output-Cubes I --Freq-NBand=%i' % channels
#
#    if polcubemode:
#        runcommand+=' --Output-Cubes=dD --RIME-PolMode=QU --Output-Mode=Dirty  --Freq-NBand=%i --Selection-ChanStart=%s --Selection-ChanEnd=%s' % (channels,startchan,endchan)

    if not cubemode and not polcubemode:
        runcommand+=' --Freq-NBand=2'
    if stokes:
        runcommand +=' --RIME-PolMode=%s --Output-Mode=Dirty'%stokes

    if do_decorr:
        runcommand += ' --RIME-DecorrMode=FT'

    if cleanmode == 'SSD':
        runcommand += ' --SSDClean-SSDSolvePars [S,Alpha] --SSDClean-BICFactor 0'
    if clusterfile is not None:
        runcommand += ' --Facets-CatNodes=%s' % clusterfile
    if automask:
        runcommand += ' --Mask-Auto=1 --Mask-SigTh=%.2f' % automask_threshold
    if cleanmask is not None:
        runcommand += ' --Mask-External=%s'%cleanmask
    if applysols is not None:
        if normalization is not None:
            if normalization[:3]=='Abs':
                normalization='Mean'+normalization # backward compat. hack
            runcommand += ' --DDESolutions-GlobalNorm='+normalization
        runcommand += ' --DDESolutions-DDModeGrid=%s --DDESolutions-DDModeDeGrid=%s --DDESolutions-DDSols=%s'%(applysols,applysols,ddsols)
    if use_dicomodel:
        if dicomodel_base is not None:
            runcommand += ' --Predict-InitDicoModel=%s.DicoModel' % dicomodel_base
        else:
            raise RuntimeError('use_dicomodel is set but no dicomodel supplied')
    if dicomodel_base is None and use_dicomodel:
        raise RuntimeError('that s wrong')
        
    if threshold is not None:
        runcommand += ' --Deconv-FluxThreshold=%f'%threshold
    if uvrange is not None:
        runcommand += ' --Selection-UVRangeKm=[%f,%f]' % (uvrange[0],uvrange[1])
    if dirty_from_resid and reuse_dirty:
        raise RuntimeError('Cannot combine reuse_dirty and dirty_from_resid')
    if dirty_from_resid:
        # possible that crashes could destroy the cache, so need to check
        if os.path.exists(cache_dir+'/'+mslist+'.ddfcache/LastResidual'):
            runcommand += ' --Cache-Dirty forceresidual'
    if reuse_dirty:
        if os.path.exists(cache_dir+'/'+mslist+'.ddfcache/Dirty'):
            runcommand += ' --Cache-Dirty forcedirty'
    if reuse_psf:
        if os.path.exists(cache_dir+'/'+mslist+'.ddfcache/PSF'):
            runcommand += ' --Cache-PSF force'

    if HMPsize is not None:
        runcommand += ' --GAClean-MinSizeInit=%i' % HMPsize

    #if options['nobar']:
    #    runcommand += ' --Log-Boring=1'

    if smooth:
        runcommand += ' --Beam-Smooth=1'

    if predict_column is not None:
        runcommand += ' --Predict-ColName=%s' % predict_column
        
    if phasecenter is not None:
        runcommand += " --Image-PhaseCenterRADEC=[%s,%s]"%(phasecenter[0],phasecenter[1])

    if conditional_clearcache:
        clearcache(mslist,options)

    s.add(runcommand, log=logfile, commandType='singularity')
    s.run(check=True)

    return imagename


def smooth_solutions(mslist, ddsols, catcher=None, dryrun=False, SkipSmooth=False, SigmaFilterOutliers=None):
    filenames=[l.strip() for l in open(mslist,'r').readlines()]
    full_sollist = []
    start_times = []
    SolsDir=o["SolsDir"]
    if SolsDir is None or SolsDir=="":
        for fname in filenames:
            solname =fname+'/killMS.'+ddsols+'.sols.npz'
            t0,t1 = get_solutions_timerange(solname)
            start_times.append(t0)
            full_sollist.append(solname)
            f.write('%s\n'%(solname))
    else:
        for fname in filenames:
            MSName=os.path.abspath(fname).split("/")[-1]
            solname =os.path.abspath(SolsDir)+"/"+MSName+'/killMS.'+ddsols+'.sols.npz'
            t0,t1 = get_solutions_timerange(solname)
            start_times.append(t0)
            full_sollist.append(solname)

    Ustart_times = np.unique(start_times)

    for start_time in Ustart_times:
        with open('solslist_%s.txt'%start_time,'w') as f:
            for i in range(0,len(full_sollist)):
                if start_times[i] == start_time:
                    solname = full_sollist[i]
                    f.write('%s\n'%(solname))
        
        checkname='%s_%s_merged.npz'%(ddsols,start_time)
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running MergeSols step')
        else:
            ss='MergeSols.py --SolsFilesIn=solslist_%s.txt --SolFileOut=%s_%s_merged.npz'%(start_time,ddsols,start_time)
            if SigmaFilterOutliers:
                ss+=" --SigmaFilterOutliers %f"%SigmaFilterOutliers
            run(ss,dryrun=dryrun)
            
        checkname='%s_%s_smoothed.npz'%(ddsols,start_time)
        if o['restart'] and os.path.isfile(checkname):
            warn('Solutions file '+checkname+' already exists, not running SmoothSols step')
        elif SkipSmooth:
            warn('Skipping smoothing Solutions file')
        else:
            run('SmoothSols.py --SolsFileIn=%s_%s_merged.npz --SolsFileOut=%s_%s_smoothed.npz --InterpMode=%s'%(ddsols,start_time,ddsols,start_time,o['smoothingtype']),dryrun=dryrun)

        smoothoutname='%s_%s_smoothed.npz'%(ddsols,start_time)

        #if InterpToMSListFreqs:
        #    interp_outname="%s_%s_interp.npz"%(smoothoutname,start_time)
        #    checkname=interp_outname
        #    if o['restart'] and os.path.isfile(checkname):
        #        warn('Solutions file '+checkname+' already exists, not running MergeSols step')
        #    else:
        #        command="InterpSols.py --SolsFileIn %s --SolsFileOut %s --MSOutFreq %s"%(smoothoutname,interp_outname,InterpToMSListFreqs)
        #        run(command,dryrun=dryrun)
        
        for i in range(0,len(full_sollist)):
            if start_times[i] == start_time:
                if not SkipSmooth:
                    symsolname = full_sollist[i].replace(ddsols,ddsols+'_smoothed')
                else:
                    symsolname = full_sollist[i].replace(ddsols,ddsols+'_merged')                 
                # always overwrite the symlink to allow the dataset to move -- costs nothing
                if os.path.islink(symsolname):
                    logger.warning('Symlink ' + symsolname + ' already exists, recreating')
                    os.unlink(symsolname)

                if not SkipSmooth:
                    os.symlink(os.path.abspath('%s_%s_smoothed.npz'%(ddsols,start_time)),symsolname)
                else:
                    os.symlink(os.path.abspath('%s_%s_merged.npz'%(ddsols,start_time)),symsolname)
                    
                    
        if SkipSmooth:
            outname = ddsols + '_merged'
        else:
            outname = ddsols + '_smoothed'

    return outname
