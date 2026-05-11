The LiLF parset
========================

LiLF is able to use a parset file. The default parset file can be found in ``LiLF.lib_util`` and is shown below as well:

.. code-block::
 
  [LOFAR_cal]
  parset_dir = /path/to/parset/dir/LOFAR_cal
  data_dir = data-bkp/
  skymodel = 
  imaging = False
  fillmissingedges = True
  less_aggressive_flag = False
  develop = False
  use_gnss = False
  use_shm = False
  
  [LOFAR_3c]
  parset_dir = /path/to/parset/dir/LOFAR_3c
  
  [LOFAR_timesplit]
  parset_dir = /path/to/parset/dir/LOFAR_timesplit
  data_dir = data-bkp/
  cal_dir = 
  ngroups = 1
  initc = 0
  fillmissingedges = True
  apply_fr = False
  no_aoflagger = False
  ateam_clip = 
  use_gnss = False
  
  [LOFAR_ddparallel]
  parset_dir = /path/to/parset/dir/LOFAR_ddparallel
  maxiter = 2
  subfield = 
  subfield_min_flux = 20
  ph_sol_mode = phase
  remove3c = True
  fulljones = False
  min_facets = 
  max_facets = 
  develop = False
  data_dir = 
  use_shm = False
  
  [LOFAR_quick-self]
  parset_dir = /path/to/parset/dir/LOFAR_quick-self
  
  [LOFAR_ateam]
  parset_dir = /path/to/parset/dir/LOFAR_ateam
  
  [LOFAR_m87]
  parset_dir = /path/to/parset/dir/LOFAR_m87
  data_dir = ./
  updateweights = False
  skipmodel = False
  model_dir = 
  
  [LOFAR_quality]
  parset_dir = /path/to/parset/dir/LOFAR_quality
  ddparallel_dir = ddparallel
  ddserial_dir = ddserial
  
  [deprecated]
  parset_dir = /path/to/parset/dir/deprecated
  
  [LOFAR_extract]
  parset_dir = /path/to/parset/dir/LOFAR_extract
  max_niter = 10
  subtract_region = 
  ph_sol_mode = phase
  amp_sol_mode = diagonal
  beam_cut = 0.3
  no_selfcal = False
  ampcal = auto
  extractregion = target.reg
  
  [LOFAR_preprocess]
  parset_dir = /path/to/parset/dir/LOFAR_preprocess
  fix_table = True
  renameavg = True
  keep_is = True
  backup_full_res = False
  demix_sources = 
  demix_skymodel = 
  demix_field_skymodel = LOTSS-DR3
  run_aoflagger = False
  tar = True
  data_dir = 
  
  [LOFAR_virgo]
  parset_dir = /path/to/parset/dir/LOFAR_virgo
  cal_dir = 
  data_dir = ./
  
  [LOFAR_ddserial]
  parset_dir = /path/to/parset/dir/LOFAR_ddserial
  maxiter = 1
  mincalflux60 = 0.8
  solve_amp = True
  use_shm = False
  target_dir = 
  manual_dd_cal = 
  develop = False
  
  [flag]
  stations = 
  antennas = 
  
  [model]
  sourcedb = 
  fits_model = 
  apparent = False
  userreg = 
  
  [PiLL]
  working_dir = /data2/user/scratch/
  redo_cal = False
  download_file = 
  project = 
  target = 
  obsid = 
  minmaxhrs = 0,9999
  logfile = 
