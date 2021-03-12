# LiLF
Library for Low Frequencies

LiLF is a set of functions linked together in pipelines for the reduction of low-frequency interferometric data in radio astronomy. It is built upon LOFAR software. LiLF can be used on both LOFAR and uGMRT data.


LOFAR: http://lofar.org/
uGMRT: http://www.ncra.tifr.res.in/ncra/gmrt

### Files:
- ~/.stagingrc >> with the login and pass for LTA archive
- ~/.surveys >> with the pass for the surveydb

### Container:
check here if you want a container that can run the pipeline:
LiLF/container/docker_build.sh

# lilf.config specifications:

### PiLL
working_dir: dir path [./] # working directory

### LOFAR_download
fix_table: bool [True] # fix bug in some old observations

renameavg: bool [True]

flag_elev: bool [True] # flag anything below 30 deg elevation

keep_IS: bool [False] # keep international stations
    
### LOFAR_demix
data_dir: dir path [../cals-bkp/]

demix_model: skaydb file [/home/fdg/scripts/model/demix_all.skydb]

### LOFAR_cal
imaging: bool [False] # perform test imaging of the calibrator data

skymodel: skydb file ["LiLF_dir"/models/calib-simple.skydb']

data_dir: dir path [../cals-bkp/]
    
### LOFAR_timesplit
data_dir: dir path [../tgts-bkp/]

cal_dir: dir path [../cals/]

ngroups: int [1] # number of frequency groups to create

initc: int [0] # init number for time chunk numbering

### LOFAR_self

### LOFAR_dd-serial
maxIter: int [2] # iterations

minCalFlux60: float [1.5] # apparent flux [Jy] at 60 MHz to be considered a calibrator

removeExtendedCutoff: float [0.0001] # remove extended sources from possible DD-calibrator list

### LOFAR_facet_self
maxniter: int [10]
