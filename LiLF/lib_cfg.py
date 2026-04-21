# Configuration loader for the LiLF pipeline.
# Reads an optional LiLF.conf file from the current or parent directory and
# fills in all pipeline sections with their default values.

import glob
import os
from configparser import ConfigParser

from LiLF.lib_log import logger


def getParset(parsetFile=''):
    """
    Locate and parse the LiLF configuration file, then populate every pipeline
    section with default values for any options not already set.

    Parameters
    ----------
    parsetFile : str, optional
        Explicit path to a config file.  When omitted (or empty), the function
        searches for 'LiLF.conf*' in the current and parent directories.

    Returns
    -------
    ConfigParser
        Fully-populated configuration object.
    """
    def add_default(section, option, val):
        """Set `option` in `section` only when it is not already present."""
        if not config.has_option(section, option):
            config.set(section, option, val)

    # Auto-discover config file when none is given explicitly.
    if not parsetFile:
        matched_conf_files = (
            glob.glob('[Ll][Ii][Ll][Ff].conf*') +
            glob.glob('../[Ll][Ii][Ll][Ff].conf*')
        )
        if len(matched_conf_files) > 1:
            raise LookupError(f'Found more than one configuration file: {matched_conf_files}')
        elif len(matched_conf_files) == 1:
            parsetFile = matched_conf_files[0]
            logger.info(f'Found config file: {parsetFile}')

    config = ConfigParser(defaults=None)
    # Only attempt to read a file when one was actually found; config.read('')
    # would silently try to open a file with an empty name.
    if parsetFile:
        config.read(parsetFile)

    # Register each parset sub-directory as a pipeline section and record its
    # default parset_dir so pipelines can locate their parset templates.
    parsets_root = os.path.join(os.path.dirname(__file__), '..', 'parsets')
    for pipeline_path in glob.glob(os.path.join(parsets_root, '*')):
        pipeline_name = os.path.basename(pipeline_path)
        if not config.has_section(pipeline_name):
            config.add_section(pipeline_name)
        if not config.has_option(pipeline_name, 'parset_dir'):
            config.set(pipeline_name, 'parset_dir',
                       os.path.join(parsets_root, pipeline_name))

    # Add shared sections not tied to a specific pipeline.
    if not config.has_section('flag'):
        config.add_section('flag')
    if not config.has_section('model'):
        config.add_section('model')
    if not config.has_section('PiLL'):
        config.add_section('PiLL')

    # --- PiLL (top-level pipeline) ---
    add_default('PiLL', 'working_dir', os.getcwd())
    add_default('PiLL', 'redo_cal', 'False')        # re-do the calibrator even when it is in the archive
    add_default('PiLL', 'download_file', '')         # html.txt file to use instead of staging
    add_default('PiLL', 'project', '')
    add_default('PiLL', 'target', '')
    add_default('PiLL', 'obsid', '')                 # unique observation ID
    add_default('PiLL', 'minmaxhrs', '0,9999')      # min and max hours for an obs to be selected
    add_default('PiLL', 'logfile', '')               # logfile for PiLL
    # --- preprocess ---
    add_default('LOFAR_preprocess', 'fix_table', 'True')               # fix bug present in some old observations
    add_default('LOFAR_preprocess', 'renameavg', 'True')
    add_default('LOFAR_preprocess', 'keep_IS', 'True')
    add_default('LOFAR_preprocess', 'backup_full_res', 'False')
    add_default('LOFAR_preprocess', 'demix_sources', '')                # demix sources in these patches, e.g. [VirA,TauA]; default: no demix
    add_default('LOFAR_preprocess', 'demix_skymodel', '')               # non-default demix skymodel; empty = use built-in
    add_default('LOFAR_preprocess', 'demix_field_skymodel', 'LOTSS-DR3') # custom target skymodel; set to '' to skip target contribution
    add_default('LOFAR_preprocess', 'run_aoflagger', 'False')          # run aoflagger on individual sub-bands (only when not done by the observatory)
    add_default('LOFAR_preprocess', 'tar', 'True')                      # tar MS files at the end
    add_default('LOFAR_preprocess', 'data_dir', '')
    # --- cal ---
    add_default('LOFAR_cal', 'data_dir', 'data-bkp/')
    add_default('LOFAR_cal', 'skymodel', '')               # default: calib-simple.skydb (LBA) or calib-hba.skydb (HBA)
    add_default('LOFAR_cal', 'imaging', 'False')
    add_default('LOFAR_cal', 'fillmissingedges', 'True')
    add_default('LOFAR_cal', 'less_aggressive_flag', 'False') # softer flagging for data with alternating SBs or many flagged points
    add_default('LOFAR_cal', 'develop', 'False')           # if True, skip file deletion (useful for debugging)
    add_default('LOFAR_cal', 'use_GNSS', 'False')          # use GNSS satellite data to pre-correct TEC and Faraday rotation
    add_default('LOFAR_cal', 'use_shm', 'False')           # use /dev/shm for temporary files when available
    # add_default('LOFAR_cal', 'beam_model', 'hamaker')
    # --- timesplit ---
    add_default('LOFAR_timesplit', 'data_dir', 'data-bkp/')
    add_default('LOFAR_timesplit', 'cal_dir', '')          # default: check repository, then ../obsid_3[c|C]*
    add_default('LOFAR_timesplit', 'ngroups', '1')
    add_default('LOFAR_timesplit', 'initc', '0')
    add_default('LOFAR_timesplit', 'fillmissingedges', 'True')
    add_default('LOFAR_timesplit', 'apply_fr', 'False')    # also transfer rotation-measure solutions (e.g. nearby calibrator)
    add_default('LOFAR_timesplit', 'no_aoflagger', 'False') # skip aoflagger (e.g. for A-Team source observations)
    add_default('LOFAR_timesplit', 'ateam_clip', '')        # e.g. [CygA,CasA] — clip listed A-Team sources not demixed by observatory
    add_default('LOFAR_timesplit', 'use_GNSS', 'False')    # use GNSS satellite data to pre-correct TEC and Faraday rotation
    # --- ddparallel ---
    add_default('LOFAR_ddparallel', 'maxIter', '2')
    add_default('LOFAR_ddparallel', 'subfield', '')         # optional ds9 box region for sub-field; '' = auto-detect via subfield_min_flux
    add_default('LOFAR_ddparallel', 'subfield_min_flux', '20') # minimum flux (Jy) within the calibration sub-field
    add_default('LOFAR_ddparallel', 'ph_sol_mode', 'phase') # phase solution mode: 'phase' or 'tecandphase'
    add_default('LOFAR_ddparallel', 'remove3c', 'True')
    add_default('LOFAR_ddparallel', 'fulljones', 'False')
    add_default('LOFAR_ddparallel', 'min_facets', '')
    add_default('LOFAR_ddparallel', 'max_facets', '')
    add_default('LOFAR_ddparallel', 'develop', 'False')     # if True, produce extra debug images (slower)
    add_default('LOFAR_ddparallel', 'data_dir', '')
    add_default('LOFAR_ddparallel', 'use_shm', 'False')     # use /dev/shm for temporary files when available
    # --- ddserial ---
    add_default('LOFAR_ddserial', 'maxIter', '1')
    add_default('LOFAR_ddserial', 'minCalFlux60', '0.8')
    add_default('LOFAR_ddserial', 'solve_amp', 'True')      # set False to skip amplitude solutions
    add_default('LOFAR_ddserial', 'use_shm', 'False')       # use /dev/shm for temporary files when available
    add_default('LOFAR_ddserial', 'use_shm_ddcal', 'True')  # use /dev/shm for temporary dd-cal files when available
    add_default('LOFAR_ddserial', 'target_dir', '')          # target direction as 'ra,dec'
    add_default('LOFAR_ddserial', 'manual_dd_cal', '')
    add_default('LOFAR_ddserial', 'develop', 'False')        # if True, produce extra debug images (slower)
    # add_default('LOFAR_ddserial', 'solve_tec', 'False')   # default: scalarphase; set True to solve for TEC instead
    # --- extract ---
    add_default('LOFAR_extract', 'max_niter', '10')
    add_default('LOFAR_extract', 'subtract_region', '')     # sources inside extract-reg to keep subtracting (e.g. bright sources within a large region)
    add_default('LOFAR_extract', 'ph_sol_mode', 'phase')   # phase solution mode: 'phase' or 'tecandphase'
    add_default('LOFAR_extract', 'amp_sol_mode', 'diagonal') # amplitude solution mode: 'diagonal' or 'fulljones'
    add_default('LOFAR_extract', 'beam_cut', '0.3')         # maximum pointing distance (fraction of primary beam) to include
    add_default('LOFAR_extract', 'no_selfcal', 'False')     # if True, extract only without self-calibration (e.g. for external facet_selfcal scripts)
    add_default('LOFAR_extract', 'ampcal', 'auto')
    add_default('LOFAR_extract', 'extractRegion', 'target.reg')
    # --- quality ---
    add_default('LOFAR_quality', 'ddparallel_dir', 'ddparallel')
    add_default('LOFAR_quality', 'ddserial_dir', 'ddserial')
    # --- virgo ---
    add_default('LOFAR_virgo', 'cal_dir', '')
    add_default('LOFAR_virgo', 'data_dir', './')
    # --- m87 ---
    add_default('LOFAR_m87', 'data_dir', './')
    add_default('LOFAR_m87', 'updateweights', 'False')
    add_default('LOFAR_m87', 'skipmodel', 'False')
    add_default('LOFAR_m87', 'model_dir', '')

    # --- General (shared across pipelines) ---

    # flag
    add_default('flag', 'stations', '')  # LOFAR station flags
    # model
    add_default('model', 'sourcedb', '')
    add_default('model', 'fits_model', '')
    add_default('model', 'apparent', 'False')
    add_default('model', 'userReg', '')

    return config


