import glob, os
from configparser import ConfigParser

from LiLF.lib_log import logger

def getParset(parsetFile=''):
    """
    Get parset file and return dict of values
    """

    def add_default(section, option, val):
        if not config.has_option(section, option): config.set(section, option, val)

    if parsetFile == '':
        matched_conf_files = glob.glob('[Ll][Ii][Ll][Ff].conf*') + glob.glob('../[Ll][Ii][Ll][Ff].conf*')
        if len(matched_conf_files) > 1:
            raise LookupError(f'Found more than one configuration file: {matched_conf_files}')
        elif len(matched_conf_files) == 1:
            parsetFile = matched_conf_files[0]
            logger.info(f'Found config file: {parsetFile}')
        else:
            parsetFile = os.path.dirname(__file__) + '/../LiLF.conf'
            logger.info(f'No config file found - using default: {parsetFile}')

    config = ConfigParser(defaults=None)
    config.read(parsetFile)

    # add pipeline sections and defaul parset dir:
    for pipeline in glob.glob(os.path.dirname(__file__) + '/../parsets/*'):
        pipeline = os.path.basename(pipeline)
        if not config.has_section(pipeline): config.add_section(pipeline)
        if not config.has_option(pipeline, 'parset_dir'):
            config.set(pipeline, 'parset_dir', os.path.dirname(__file__) + '/../parsets/' + pipeline)
    # add other sections
    if not config.has_section('flag'): config.add_section('flag')
    if not config.has_section('model'): config.add_section('model')
    if not config.has_section('PiLL'): config.add_section('PiLL')
    if not config.has_section("scheduler"): config.add_section("scheduler")

    
    print(config.sections())
    ### General ###
    # scheduler
    add_default('scheduler', 'use_shm', 'False')  # use /dev/shm for temporary files, if available
    add_default('scheduler', 'backend', 'local') # 'local' or 'slurm'
    add_default('scheduler', 'max_cpucores', '') # default: use all on node
    add_default('scheduler', 'slurm_max_workers', '1') # (only for slurm): max parallel workers
    add_default('scheduler', 'slurm_mem_per_cpu', '16') # (only for slurm): max parallel workers
    # flag
    add_default('flag', 'stations', '')
    # model
    add_default('model', 'sourcedb', '')
    add_default('model', 'fits_model', '')
    add_default('model', 'apparent', 'False')
    add_default('model', 'userReg', '')

    ### LOFAR ###

    # PiLL
    add_default('PiLL', 'working_dir', os.getcwd())
    add_default('PiLL', 'redo_cal', 'False')  # re-do the calibrator although it is in the archive
    add_default('PiLL', 'download_file', '')  # html.txt file to use instead of staging
    add_default('PiLL', 'project', '')
    add_default('PiLL', 'target', '')
    add_default('PiLL', 'obsid', '')  # unique ID
    add_default('PiLL', 'minmaxhrs', '0,9999')  # min and max hours for an obs to be selected
    add_default('PiLL', 'logfile', '')  # logfile for PiLL
    # preprocess
    add_default('LOFAR_preprocess', 'fix_table', 'True')  # fix bug in some old observations
    add_default('LOFAR_preprocess', 'renameavg', 'True')
    add_default('LOFAR_preprocess', 'keep_IS', 'True')
    add_default('LOFAR_preprocess', 'backup_full_res', 'False')
    add_default('LOFAR_preprocess', 'demix_sources',
                '')  # Demix  sources in these patches (e.g. [VirA,TauA], default: No demix
    add_default('LOFAR_preprocess', 'demix_skymodel', '')  # Use non-default demix skymodel.
    add_default('LOFAR_preprocess', 'demix_field_skymodel',
                'LOTSS-DR3')  # Provide a custom target skymodel instead of online gsm model. Set to '' to ignore target.
    add_default('LOFAR_preprocess', 'run_aoflagger',
                'False')  # run aoflagger on individual sub-bands, only in cases where this was not one by the observatory!
    add_default('LOFAR_preprocess', 'tar', 'True')  # Tar MS files at the end
    add_default('LOFAR_preprocess', 'data_dir', '')
    # cal
    add_default('LOFAR_cal', 'data_dir', 'data-bkp/')
    add_default('LOFAR_cal', 'skymodel', '')  # by default use calib-simple.skydb for LBA and calib-hba.skydb for HBA
    add_default('LOFAR_cal', 'imaging', 'False')
    add_default('LOFAR_cal', 'fillmissingedges', 'True')
    add_default('LOFAR_cal', 'less_aggressive_flag',
                'False')  # change flagging so that we can handle data with alternating SBs only or many flagged points
    add_default('LOFAR_cal', 'develop', 'False')  # if true prevents the deletion of files
    add_default('LOFAR_cal', 'use_GNSS', 'False')  # Use GNSS satellite data for pre correcting TEC and FR
    # timesplit
    add_default('LOFAR_timesplit', 'data_dir', 'data-bkp/')
    add_default('LOFAR_timesplit', 'cal_dir', '')  # by default the repository is tested, otherwise ../obsid_3[c|C]*
    add_default('LOFAR_timesplit', 'ngroups', '1')
    add_default('LOFAR_timesplit', 'initc', '0')
    add_default('LOFAR_timesplit', 'fillmissingedges', 'True')
    add_default('LOFAR_timesplit', 'apply_fr',
                'False')  # Also transfer rotationmeasure sols? (E.g. for nearby calibrator and target)
    add_default('LOFAR_timesplit', 'no_aoflagger',
                'False')  # TEST: Skip aoflagger (e.g. for observations of A-Team sources)
    add_default('LOFAR_timesplit', 'ateam_clip',
                '')  # [CygA, CasA] or [CasA] or [CygA] or '' - the code clips the specified ateams (if not demixed from the observatory)
    add_default('LOFAR_timesplit', 'use_GNSS', 'False')  # Use GNSS satellite data for pre correcting TEC and FR
    # ddparallel
    add_default('LOFAR_ddparallel', 'maxIter', '2')
    add_default('LOFAR_ddparallel', 'subfield',
                '')  # possible to provide a ds9 box region customized sub-field. DEfault='' -> Automated detection using subfield_min_flux.
    add_default('LOFAR_ddparallel', 'subfield_min_flux', '20')  # min flux within calibration subfield
    add_default('LOFAR_ddparallel', 'ph_sol_mode', 'phase')  # phase or tecandphase
    add_default('LOFAR_ddparallel', 'remove3c', 'True')
    add_default('LOFAR_ddparallel', 'fulljones', 'False')
    add_default('LOFAR_ddparallel', 'min_facets', '')
    add_default('LOFAR_ddparallel', 'max_facets', '')
    add_default('LOFAR_ddparallel', 'develop', 'False')  # if true make more debug images (slower)
    add_default('LOFAR_ddparallel', 'data_dir', '')
    # ddserial
    add_default('LOFAR_ddserial', 'maxIter', '1')
    add_default('LOFAR_ddserial', 'minCalFlux60', '0.8')
    add_default('LOFAR_ddserial', 'solve_amp', 'True')  # to disable amp sols
    add_default('LOFAR_ddserial', 'target_dir', '')  # ra,dec
    add_default('LOFAR_ddserial', 'manual_dd_cal', '')
    add_default('LOFAR_ddparallel', 'use_shm', 'False')  # use /dev/shm for temporary files, if available
    add_default('LOFAR_ddserial', 'develop', 'False')  # if true make more debug images (slower)
    # add_default('LOFAR_ddserial', 'solve_tec', 'False') # per default, solve each dd for scalarphase. if solve_tec==True, solve for TEC instead.
    # extract
    add_default('LOFAR_extract', 'max_niter', '10')
    add_default('LOFAR_extract', 'subtract_region',
                '')  # Sources inside extract-reg that should still be subtracted! Use this e.g. for individual problematic sources in a large extractReg
    add_default('LOFAR_extract', 'ph_sol_mode', 'phase')  # tecandphase, phase
    add_default('LOFAR_extract', 'amp_sol_mode', 'diagonal')  # diagonal, fulljones
    add_default('LOFAR_extract', 'beam_cut', '0.3')  # up to which distance a pointing will be considered
    add_default('LOFAR_extract', 'no_selfcal',
                'False')  # just extract the data, do not perform selfcal - use this if u want to use e.g. Reinout van Weeren's facet_seflcal script
    add_default('LOFAR_extract', 'ampcal', 'auto')
    add_default('LOFAR_extract', 'extractRegion', 'target.reg')
    # quality
    add_default('LOFAR_quality', 'ddparallel_dir', 'ddparallel')
    add_default('LOFAR_quality', 'ddserial_dir', 'ddserial')
    # virgo
    add_default('LOFAR_virgo', 'cal_dir', '')
    add_default('LOFAR_virgo', 'data_dir', './')
    # m87
    add_default('LOFAR_m87', 'data_dir', './')
    add_default('LOFAR_m87', 'updateweights', 'False')
    add_default('LOFAR_m87', 'skipmodel', 'False')
    add_default('LOFAR_m87', 'model_dir', '')


    return config
