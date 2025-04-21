# LiLF
Library for Low Frequencies

LiLF is a set of functions linked together in pipelines for the reduction of low-frequency interferometric data in radio astronomy. It is built upon LOFAR software. LiLF can be used on both LOFAR and uGMRT data.

- LOFAR: <https://www.astron.nl/telescopes/lofar/>
- uGMRT: <http://www.ncra.tifr.res.in/ncra/gmrt>

### Files:
- ~/.stagingrc >> with the login and pass for LTA archive, for example:
```
user = <username>
password = <password>
```
- ~/.surveys >> with the pass for the surveydb

- ~/.awe/Environment.cfg >> if you want to use the LiLF/scripts/LOFAR_stager.py to stage and download, your LTA credentials should be also in this file, for example:
```
[global]
database_user : <your username>
database_password : <your password>
```

### Container:
Check here if you want a container that can run the pipeline:
LiLF/container/docker_build.sh

# LBA data reduction How-To:
To calibrate LOFAR LBA data is possible to use the script PILL.py (in the pipeline dir) that automatically does all the different pipeline steps, or to run the single steps by hand to check the intermediate results.

If you use these scripts, please cite:
- [de Gasperin+ 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...622A...5D/abstract) (DIE calibration)
- [de Gasperin+ 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..85D/abstract) (DDE calibration)

Information on the ionosphere systematic effects can be found here:
- [deGasperin+ 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...615A.179D/abstract)

If you demixed A-team sources, here is the paper describing the models:
- [de Gasperin+ 2020b](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.150D/abstract)

### Environment
The preferred way is to use the singularity container as described above, enter in the singularity and go to your working directory.
Check the path for LiLF, let's say it is /opt/LiLF/, and include it in your ~/.bashrc as follow:

`export PYTHONPATH=/opt/LiLF:$PYTHONPATH`

`export PATH=/opt/LiLF/scripts:$PATH`

I also reccommended to set:

`ulimit -n 8192`

as the pipeline requires to open many files at the same time.

### Run the pipelines:

* To use PILL run

`python3 /opt/LiLF/pipelines/PiLL.py`

* To run the pipeline step-by-step follow these commands:
    1. On you working directory create a `Download` directory and put here the html.txt files obtained from a data staging request on Long Term Archive (LTA).  Then run: `python3 /opt/LiLF/pipelines/LOFAR_preprocess.py` to download the data from LTA, unpack, average to 4 chan/sb and 4 sec and finally arrange the data in sub-directories that you can find in `Download/mss`. The subdirectories are called `id000_CAL` and `id000_TARGET`, where 000 is the id of your observation and CAL and TARGET are the name of the calibrator and target. If your observation is split in more than one night, you will have a calibrator and target directory for every observation. Inside that directories you will find all the ms files. Copy them in a directory called data-bkp (Don't change the name otherwise the pipeline doesn't find the ms files). So to summarize you will have `Download/mss/id000_CAL/data-bkp` and `Download/mss/id000_TARGET/data-bkp`.

    2. In your cal directory `Download/mss/id000_CAL/` run the calibrator pipeline that will estimate the contribution of systematic effects on your observations: `python3 /opt/LiLF/pipelines/LOFAR_cal.py`. Do it for every observations if you have more than one.

    3. In your target directory `Download/mss/id000_TARGET` run the split pipeline to apply the calibrator solutions and split the data to MS of 1h to run the next steps in parallel: `python3 /opt/LiLF/pipelines/LOFAR_timesplit.py`. You can find the new MS in `Download/mss/id000_TARGET/mss` named as TC00.MS, TC01.MS etc., one for every hour of observation. Now if you have more than one observation, copy all the MS files in a single directory as `Download/mss/TARGET/mss`, pay attention to rename the files with an increasing number, as for each observation the name start from T00.MS.

    4. In your target directory, run the self pipeline that performs the direction-indipendent calibration: `python3 /opt/LiLF/pipelines/LOFAR_ddparallel.py`. In the self directory you can find some plots useful to understand the quality of the ionosphere and calibration solutions, you can find some examples in the papers mentioned below. The images are in the `img` directory.

    5. In your target directory, run the dd-serial pipeline that performs the direction-dependent calibration: `python3 /opt/LiLF/pipelines/LOFAR_ddserial.py`. The pipeline selects the DD-calibrators, calibrates them one after the other (check the plots in `ddcal/` and the images of the varius steps of self-calibration), then it transfers the solutions to the associated facet and with DDFacet creates an image of the widefield. It performs two major cycles called c00 and c01. In the `/ddcal/c00/skymodels/` directory you can find the `all-c00.reg` file that indicates all the source selected as DD-calibrators, load it on the image you obtain from the previous step to check which sources it selects, they are indicated with a red circle and named ddcal000. To have a good image of your source is important that it is selected as a calibrator. The images of the single calibrators are named as `ddcalM-c01-ddcal0059-cdd00-MFS-image.fits, ddcalM-c01-ddcal0059-cdd01-MFS-image.fits` where cdd are the different steps of selfcal. You can use the best of them as final image of your source. The image of the widefield instead is `wideDD-c01.app.restored.fits`. You can re-image it using DDFacet with your prefered parameters (for example try to use --Deconv-Mode SSD).


# Extraction of LBA data:

Usage: `python LiLF/pipelines/LOFAR_extract.py -p [/path/to/observation] --radec [RA and DEC in deg]`

You can extract a target of interest to improve selfcalibration and try to fix ionospheric effects.
If you wish to extract only one target, simply run the command above indicating the path to the directory of
the observation [-p], which must contain the subdirectories /ddcal and /mss-avg obtained from the calibration pipeline (e.g. if one
has /example1/example2/myobs/ddcal, it will be [-p /example1/example2]), and RA and DEC where to center the extraction (--radec).
Optionally one can add redshift [--z] and target name [--name].
The pipeline is able to process multiple pointings of the same target altogether: it will look in the path specified with -p for directories
calibrated with LiLF (see previous steps), check if the input coordinates are covered by a specific observation, and move on to the next ones.
All observations for which the beam sensitivity is above 30% at the coordinates position will be used for the extraction.
The code will then create an extraction region based on the flux density within a given area (larger if low flux and vice-versa), run the selfcalibration and produce images. Outputs will be stored in the /img subdirectory within the target directory, while extracted .MS files can be found in the /mss-extract subdirectory. A final, nominal-resolution image will be produced (`extractM-final-MFS-image.fits`), as well as high-resolution, low-resolution and source-subtracted images. The source subtraction is still on beta version, please check it carefully before using the relative image.
Multiple default parameters can be changed through the command line - see below for a brief description.

If you wish to extract more than one target, prepare a .fits file with a minimum of 4 columns: Name, RA, DEC, z (the column names must exactly match these ones, but the order can be different). Run the extraction script through:

`python LiLF/pipelines/LOFAR_extract.py -p [/path/to/observation] --l [/path/to/fits/file]`

Coordinates must be in degrees, while the purpose of the name is just to create a subdirectory in which all the results will be stored. Redshift is mandatory only if one wants to perform the subtraction of compact sources, otherwise just put -99.
An optional column can also be added with the name EXTREG. This must be a .reg file that the pipeline uses as extraction region. The script is usually able to create a good one by itself, use this option only if needed - always try a run without it first.
Another optional .reg file can be provided in a column named 'MASKREG'; the script will use it as cleaning mask.
Finally, one can specify a subtraction region under the column 'SUBREG'; every source within this region will be subtracted.
The extraction will be run for each object in a different directory named after the 'Name' column.  

### Command line parameters

`-p`, `--path`: Path where to look for observations. It must lead to a directory where subdirectories contain /ddcal and /mss-avg derived from calibration.

`--radec`: RA/DEC where to center the extraction in deg. Use if you wish to extract only one target.

`--z`: Redshift of the target. Not necessary unless one wants to perform compact source subtraction.

`--name`: Name of the target. Will be used to create the directory containing the extracted data.

`-l`,`--list`: Name of .fits file which lists Name, RA, DEC and z. Optionally an extraction region and a mask region can be added. Use only for multiple extractions.

`--beamcut`: Beam sensitivity threshold. Default is 0.3.

`--noselfcal`: Do not perform selfcalibration.

`--extreg`: Provide an optional extraction region. If not, one will be created automatically.

`--maskreg`: Provide an optional user mask for cleaning.

`--ampcal`: Perform amplitude calibration. Can be set to True, False or auto. Default is auto.

`--ampsol`: How to solve for amplitudes. Can be set to diagonal or fulljones. Default is diagonal.

`--phsol`: How to solve for phases. Can be set to tecandphase or phase. Default is tecandphase.

`--maxniter`: Maximum number of selfcalibration cycles to perform. Default is 10.

`--subreg`: Provide an optional mask for sources one wishes to completely subtract


# lilf.config specifications:

One can create a `lilf.config` (or `lilf.conf`) file to pass specific settings to the pipeline. This file must be saved either in the working directory or in the one directly above it. A non-extensive list (for now) of parameters can be found below, they are also listed at the beginning of every script of the pipeline. For the config file to work one needs to include as first row the name of the step (for example, to put some specific options to the preprocess, the first row of the file must be `[LOFAR_preprocess]`). A few small example files are reported at the end of the documentation.

### PiLL
working_dir: dir path [./] # working directory

download_file: str # html.txt file for downloading if no project/target/obsid provided

project: str # LOFAR project

target: str # observation target

obsid: str # unique observation ID

### flag
stations: str [DE*;FR*;SE*;UK*;IE*;PL*] # For LOFAR, by default flag all IS.

antennas: str, # uGMRT

### model
sourcedb: str # model sourcedb

apparent: bool [False] # is model in apparent sky?

userReg: str # user defined region for cleaning

### LOFAR_preprocess
fix_table: bool [True] # fix bug in some old observations

renameavg: bool [True] # rename and average the data

keep_IS: bool [True] # keep the LOFAR international stations?

demix_skymodel: path/to/file.skydb

demix_sources: name(s) of the source(s) to demix # if there are two or more sources, the names must be written inside square brackets, for example [CygA, TauA]

### LOFAR_demix
This step has been implemented in the preprocess
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

### LOFAR_ddparallel

### LOFAR_ddserial
maxIter: int [2] # iterations

minCalFlux60: float [1.5] # apparent flux [Jy] at 60 MHz to be considered a calibrator

removeExtendedCutoff: float [0.0001] # remove extended sources from possible DD-calibrator list


### lilf.config examples
An example of a small lilf.config file is the following:

```
[LOFAR_preprocess]
demix_sources = CygA
demix_skymodel = /path/to/LiLF/models/demix_all.skymodel
```
this config file here adds a source to demix and passes in the skymodel.

```
[flag]
stations = RS106
```
this config file works only for the LOFAR_cal and LOFAR_timesplit steps and is applied when running one of the two. Stations allows to select the stations to flag.
