# LiLF
Library for Low Frequencies

LiLF is a set of functions linked together in pipelines for the reduction of low-frequency interferometric data in radio astronomy. It is built upon LOFAR software. LiLF can be used on both LOFAR and uGMRT data.

- LOFAR: http://lofar.org/
- uGMRT: http://www.ncra.tifr.res.in/ncra/gmrt

### Files:
- ~/.stagingrc >> with the login and pass for LTA archive, for example:
```
user = <username>
password = <password>
```
- ~/.surveys >> with the pass for the surveydb

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

`ulimit -n 4000`

as the pipeline requires to open many files at the same time.

### Run the pipelines:

* To use PILL run

`python3 /opt/LiLF/pipelines/PILL.py`

* To run the pipeline step-by-step follow these commands:
    1. On you working directory create a `Download` directory and put here the html.txt files obtained from a data staging request on Long Term Archive (LTA).  Then run: `python3 /opt/LiLF/pipelines/LOFAR_preprocess.py` to download the data from LTA, unpack, averaged to 4 chan/sb and 4 sec and finally arrange the data in sub-directories that you can find in `Download/mss` The subdirectories are called `id000_CAL` and `id000_TARGET`, where 000 is the id of your observation and CAL and TARGET are the name of the calibrator and target. If your observation is split in more than one night, you will have a calibrator and target directory for every observation. Inside that directories you will find all the ms files. Copy them in a directory called data-bkp (Don't change the name otherwise the pipeline doesn't find the ms files). So to summarize you will have `Download/mss/id000_CAL/data-bkp` and `Download/mss/id000_TARGET/data-bkp`.
    
    2. In your cal directory `Download/mss/id000_CAL/` run the calibrator pipeline that will estimate the contribution of systematic effects on your observations: `python3 /opt/LiLF/pipelines/LOFAR_cal.py`. Do it for every observations if you have more than one.
    
    3. In your target directory `Download/mss/id000_TARGET` run the split pipeline to apply the calibrator solutions and split the data to MS of 1h to run the next steps in parallel: `python3 /opt/LiLF/pipelines/LOFAR_timesplit.py`. You can find the new MS in `Download/mss/id000_TARGET/mss` named as TC00.MS, TC01.MS etc., one for every hour of observation. Now if you have more than one observation, copy all the MS files in a single directory as `Download/mss/TARGET/mss`, pay attention to rename the files with an increasing number, as for each observation the name start from T00.MS.
    
    4. In your target directory, run the self pipeline that performs the direction-indipendent calibration: `python3 /opt/LiLF/pipelines/LOFAR_self.py`. In the self directory you can find some plots useful to understand the quality of the ionosphere and calibration solutions, you can find some examples in the papers mentioned below. The images are in the `img` directory.
    
    5. In your target directory, run the dd-serial pipeline that performs the direction-dependent calibration: `python3 /opt/LiLF/pipelines/LOFAR_dd-serial.py`. The pipeline selects the DD-calibrators, calibrates them one after the other (check the plots in `ddcal/` and the images of the varius steps of self-calibration), then it transfers the solutions to the associated facet and with DDFacet creates an image of the widefield. It performs two major cycles called c00 and c01. In the `/ddcal/c00/skymodels/` directory you can find the `all-c00.reg` file that indicates all the source selected as DD-calibrators, load it on the image you obtain from the previous step to check which sources it selects, they are indicated with a red circle and named ddcal000. To have a good image of your source is important that it is selected as a calibrator. The images of the single calibrators are named as `ddcalM-c01-ddcal0059-cdd00-MFS-image.fits, ddcalM-c01-ddcal0059-cdd01-MFS-image.fits` where cdd are the different steps of selfcal. You can use the best of them as final image of your source. The image of the widefield instead is `wideDD-c01.app.restored.fits`. You can re-image it using DDFacet with your prefered parameters (for example try to use --Deconv-Mode SSD).

            
# Extraction of LBA Survey data:

Usage: python LiLF/pipelines/Cluster_extraction.py -l [.fits file name]

You can extract a source visibilities to improve selfcalibration and try to get rid of ionospheric effects.
To do so, prepare a .fits file with a minimum of 4 columns: Name, RA, DEC, z (the column names must exactly match these ones, but the order can be different). Coordinates must be in degrees, while the purpose of the name is just to create a subdirectory in which all the results will be stored. The coordinates denote where the phase center will be shifted - it is not necessary to point them at the center of the target of interest. An optional column can also be added with the name extreg. This must be a .reg file that the pipeline uses as extraction region. The script is usually able to create a good one by itself, use this option just if needed - always try a run without it first. The .fits file can contain more than one entry, the extraction will be run for each object in a different directory named after the 'Name' column.  
The pipeline is able to process multiple pointings of the same target altogether: it will look in the working directory for directories calibrated with LiLF (see previous steps), check if the provided coordinates are covered by that specific pointing, and move on to the next ones. All pointings for which the beam sensitivity does not drop below 30% at the coordinates position will be used for the extraction. The code will then create an extraction region based on the flux density within the same region (larger if low flux and vice-versa), run the selfcalibration and produce images. Outputs will be stored in the /img subdirectory within the target directory, while extracted .MS files can be found in the /mss-extract subdirectory. A final, nominal-resolution image will be produced, as well as high-resolution, low-resolution and source-subtracted images. The source subtraction is still on beta version, please check it carefully before using it.
Specific passages can be repeated without re-running the whole extraction by deleting the relative entry from the pipeline-extract.walker file within the target directory. The extraction of a single target can also be repeated by deleting the relative entry from the cluster-extraction.walker file within the working directory. IMPORTANT: in case of failures/errors, always delete the relative entry from the cluster-extraction.walker, otherwise the code will skip it.
Some parameters can be tuned by adding them in the lilf.config file:

- Mask region:  
`[model]`     
`userReg = userReg.reg`


- Force/avoid amplitude calibration (keep it to auto as per default unless needed):

  `[ampcal]`         
`ampcal = true/false/auto`


# lilf.config specifications:

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

flagelev: bool [True] # flag anything below 30 deg elevation

keep_IS: bool [False] # keep the LOFAR international stations?
    
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

### LOFAR_extract
maxniter: int [10] # maximum number of iterations to perform

extractRegion: str [target.reg] # ds9 region defining the extract target

phSolMode: str [tecandphase] # mode to use for phase solutions, either 'tecandphase' (to reduce number of free parameters) or 'phase' (to allow for more accurate solutions)
