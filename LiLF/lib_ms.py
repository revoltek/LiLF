#!/usr/bin/python

import os, sys, shutil
import numpy as np
import pyregion
from pyregion.parser_helper import Shape
from casacore import tables
from astropy.coordinates import get_sun, get_body, SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u

from LiLF import lib_util
from LiLF.lib_log import logger

# remove ires warning
from astropy.utils import iers
iers.conf.auto_download = False  

class AllMSs(object):

    def __init__(self, pathsMS, scheduler, check_flags=True, check_sun=False, min_sun_dist=10, check_consistency=False):
        """
        pathsMS:    list of MS paths
        scheduler:  scheduler object
        check_flag: if true ignore fully flagged ms
        check_sun: if true check sun distance
        min_sun_dist: if check_sun and distance from the sun < than this deg, skip
        check_consistency: it quits if the MSs are not all of the same antenna_mode and time/freq resolution
        """
        self.scheduler = scheduler

        # sort them, useful for some concatenating steps
        if len(pathsMS) == 0:
            logger.error('Cannot find MS files.')
            raise('Cannot find MS files.')

        self.mssListObj = []
        for pathMS in sorted(pathsMS):
            ms = MS(pathMS)
            if check_flags and ms.isAllFlagged(): 
                logger.warning('Skip fully flagged ms: %s' % pathMS)
            elif check_sun and ms.get_sun_dist() < min_sun_dist:
                logger.warning('Skip too close to sun (%.0f deg) ms: %s' % (ms.get_sun_dist(), pathMS))
            else:
                self.mssListObj.append(ms)
        
        if len(self.mssListObj) == 0:
            raise('ALL MS files flagged.')

        # check that antenna mode and time/freq resolutions are the same for all datasets
        if check_consistency:
            antenna_modes = [ms.getAntennaSet() for ms in self.mssListObj]
            if len(set(antenna_modes)) > 1:
                logger.error('Mixed antenna modes in AllMSs:', antenna_modes)
                sys.exit()
            time_ints = [round(ms.getTimeInt()) for ms in self.mssListObj]
            if len(set(time_ints)) > 1:
                logger.error('Mixed time intervals in AllMSs:', time_ints)
                sys.exit()
            freq_chan = [ms.getNchan() for ms in self.mssListObj]
            if len(set(freq_chan)) > 1:
                logger.error('Mixed nchan in AllMSs:', freq_chan)
                sys.exit()

        self.mssListStr = [ms.pathMS for ms in self.mssListObj]
        self.resolution = self.mssListObj[0].getResolution(check_flags=False)

        if len(self.mssListObj) > 500:
            logger.warning('Many MSs detected, using only the first to determine antenna set (HBA/LBA) and presence of IS.')
            self.isLBA = 'LBA' in self.mssListObj[0].getAntennaSet()
            self.isHBA = 'HBA' in self.mssListObj[0].getAntennaSet()
            self.hasIS = self.mssListObj[0].getMaxBL(check_flags=False) > 150e3
        else:
            self.isLBA = all(['LBA' in ms.getAntennaSet() for ms in self.mssListObj])
            self.isHBA = all(['HBA' in ms.getAntennaSet() for ms in self.mssListObj])
            self.hasIS = any([ms.getMaxBL(check_flags=False) > 150e3 for ms in self.mssListObj])


    def getListObj(self):
        """
        """
        return self.mssListObj


    def getListStr(self):
        """
        """
        return self.mssListStr

    def getNThreads(self, maxProcs=None):
        """
        Return the max number of threads per process assuming all MSs run at the same time
        with a max parallel runs of maxProcs
        """
        NumMSs = len(self.mssListStr)
        if self.scheduler.max_cpucores < NumMSs:
            NThreads = 1
        else:
            if maxProcs != None:
                Nprocs = min(maxProcs, NumMSs)
            else:
                Nprocs = NumMSs

            NThreads = round( self.scheduler.max_cpucores/Nprocs )

        return NThreads

    def getStrWsclean(self):
        """
        Return a string with all MS paths, useful for wsclean
        """
        return ' '.join(self.mssListStr)


    def getFreqs(self):
        """
        Return a list of freqs per chan per SB
        """
        freqs = [ list(ms.getFreqs()) for ms in self.mssListObj ]
        return np.array([item for sublist in freqs for item in sublist]).flatten()


    def getBandwidth(self):
        """
        Return the total span of frequency covered by this MS set
        """
        freqs = self.getFreqs()
        return freqs.max() - freqs.min()

    def getChout(self, size=4.e6):
        """
        Returns the channels-out parameter for wsclean
        size: size of each out channel in Hz
        """
        return int(round(self.getBandwidth()/(size)))

    def getMaxBL(self, check_flags=True, dutch_only=False, uvw=False):
        """
        Return length of longest baseline
        Parameters
        ----------
        check_flags: bool, check flags? default = True
        dutch_only: bool, check only dutch-dutch BL? default = False
        uvw: bool, consider also w values? default = False

        Returns
        -------
        max_bl: float
        """
        return np.max([ms.getMaxBL(check_flags=check_flags, dutch_only=dutch_only, uvw=uvw) for ms in self.getListObj()])


    def meanFractionalFlag(self):
        """
        Return the mean fractional flags (assuming all MSs have the same size)
        """
        return np.mean([ms.fractionalFlag() for ms in self.getListObj()])


    def run(self, command, log, commandType=None, maxProcs=None):
        """
        Run command 'command' of type 'commandType', and use 'log' for logger,for each MS of AllMSs.
        The command and log file path can be customised for each MS using keywords (see: 'MS.concretiseString()').
        Beware: depending on the value of 'Scheduler.max_threads' (see: lib_util.py), the commands are run in parallel.
        """
        # add max num of threads given the total jobs to run
        # e.g. in a 64 processors machine running on 16 MSs, would result in numthreads=4
        # TODO we need to figure out how to set numthreads for the slurm case...
        if commandType == 'DP3':
            if self.scheduler.backend == 'slurm':
                command += ' numthreads=16'
            else:
                command += ' numthreads='+str(self.getNThreads(maxProcs))

        for MSObject in self.mssListObj:
            commandCurrent = MSObject.concretiseString(command)
            logCurrent     = MSObject.concretiseString(log)

            self.scheduler.add(cmd = commandCurrent, log = logCurrent, commandType = commandType)

            # Provide debug output.
            #lib_util.printLineBold("commandCurrent:")
            #print (commandCurrent)
            #lib_util.printLineBold("logCurrent:")
            #print (logCurrent)

        self.scheduler.run(check = True)

    def addcol(self, newcol, fromcol, usedysco='auto', log='$nameMS_addcol.log', overwrite=True):
        """
        # TODO: it might be that if col exists and is dysco, forcing no dysco will not work. Maybe force TiledColumnStMan in such cases?
        Use DP3 to add a new data column using values from an existing column.
        Parameters
        ----------
        newcol: string, name of new column
        fromcol: string, name of existing column
        usedysco: bool or string, if bool: use dysco? if 'auto', use dysco if fromcol uses dysco.
        log: string, logfile name
        """

        check=0
        for ms_file in self.mssListStr:
            with tables.table(ms_file, ack=False) as t:
                if newcol in t.colnames() and overwrite==False:
                    logger.info(f'Column {newcol} already exists in {ms_file} and overwrite=False. Skipping..')
                    continue
                else:
                    check += 1
                    break

        if check != 0:
            sm = '' # storagemanager
            if usedysco == 'auto': # if col is dysco compressed in first MS, assume it is for all MSs
                with tables.table(self.mssListStr[0], ack=False) as t:
                    if t.getdminfo(fromcol)['TYPE'] == 'DyscoStMan':
                        sm = 'dysco'
            elif usedysco:
                sm = 'dysco'

            self.run(f'DP3 msin=$pathMS msin.datacolumn={fromcol} msout=. msout.datacolumn={newcol} \
                     msout.storagemanager={sm} steps=[]', log=log, commandType="DP3")
            
    def deletecol(self, col):
        """
        Use DP3 to delete a column.
        Parameters
        ----------
        col: string, name of column to delete
        log: string, logfile name
        """
        for ms_file in self.mssListStr:
            with tables.table(ms_file, ack=False, readonly=False) as t:
                if col not in t.colnames():
                    logger.info(f'deletecol: Column {col} does not exist in {ms_file}. Skipping..')
                    continue
                else:
                    logger.debug(f'Deleting column {col} from {ms_file}....')
                    t.removecols(col)

    def run_Blsmooth(self, incol='DATA', outcol='SMOOTHED_DATA', ionf='auto',  notime=False, nofreq=False, logstr='smooth'):
        """
        Execute BLsmooth incol-> outcol on a group of MSs in a way that tries to maximise resource efficiency.

        Parameters
        ----------
        ionf: float, ionofactor. Indicator of ionosphere strength,
            default='auto' -> 0.2e-3 for IS, else 1e-3
        incol: str, input column name. Default: 'DATA'
        outcol: str, output column name. Default: 'SMOOTHED_DATA'
        notime: bool, do not smooth in time?
        nofreq: bool, do not smooth in freq?
        logstr: str, logfile name suffix. Default: 'smooth'.
        """

        if ionf == 'auto': ionf = .2e-3 if self.hasIS else 1e-3

        # if multiple MSs - parallelize on MSs before adding chunks
        n_ms = len(self.getListObj())
        maxProcs = min([n_ms, 8])
        # possibly, we need to reduce the maxProcs for IS observations? Let's see.
        # if self.hasIS:
        #     maxProcs = 1

        # calculate the "size" of a single MS (~times*freq*BL). We assume that all MSs have the same size here.
        ms_size = self.mssListObj[0].getNtime() # N_times
        ms_size *= self.mssListObj[0].getNchan() # N_chan
        ms_size *= len(self.mssListObj[0].getAntennas())*(len(self.mssListObj[0].getAntennas())-1)/2 # N_BL

        # normalize by 1h dutch BL with 8chan/122SB/4s
        reference_size = 900 * 976 * 38*37/2
        # of such a ref MS, we can run 8 threads / 4 chunks in parallel
        # TODO: If this runs out of memory, we need to increase the prefactor (4) below
        chunks = 4 * ms_size / reference_size
        # if we have less than 8 threads, we can also reduce the number of chunks
        chunks *= maxProcs/8
        # make sure chunks >= 1 and integer
        if chunks < 1: chunks = 1
        chunks = int(np.round(chunks))

        ncpu = round(self.scheduler.max_cpucores / maxProcs)  # cpu max_proc / threads

        extra_flags = ''
        if notime: extra_flags += ' -t'
        if nofreq: extra_flags += ' -q'

        logger.info(f'BL-smooth: chunks={chunks}; ncpu={ncpu}; max processes={maxProcs}...')
        self.run(f'BLsmooth.py -c {chunks} -n {ncpu} -f {ionf} -r -i {incol} -o {outcol} {extra_flags} $pathMS',
                log=f'$nameMS_{logstr}.log', commandType='python', maxProcs=maxProcs)

    def print_HAcov(self, png=None):
        """
        some info on the MSs
        """
        has = []; elevs = []
        for ms in self.mssListObj:
            logger.info('%s (%s): Hour angle: %.1f hrs - Elev: %.2f (Sun distance: %.0f; Jupiter distance: %.0f)' % \
                        (ms.nameMS,ms.get_time().iso,ms.get_hour_angle(),ms.get_elev(),ms.get_sun_dist(),ms.get_jupiter_dist()))
            has.append(ms.get_hour_angle())
            elevs.append(ms.get_elev())

        if png is not None:
            import matplotlib.pyplot as pl
            pl.figure(figsize=(6,6))
            ax1 = pl.gca()
            ax1.plot(has, elevs, 'ko')
            ax1.set_xlabel('HA [hrs]')
            ax1.set_ylabel('elevs [deg]')
            logger.debug('Save plot: %s' % png)
            pl.savefig(png)


class MS(object):

    def __init__(self, pathMS):
        """
        pathMS:        path of the MS
        pathDirectory: path of the parent directory of the MS
        nameMS:        name of the MS, without parent directories and extension (which is assumed to be ".MS" always)
        """
        if pathMS[-1] == "/": pathMS = pathMS[:-1]
        self.setPathVariables(pathMS)
        
    def get_telescope_coords(self):
        """
        Return astropy coords of telescope location
        """
        telescope = self.getTelescope()
        if telescope == 'LOFAR':
            telescope_coords = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m)
        elif telescope == 'GMRT':
            telescope_coords = EarthLocation(lat=19.0948*u.deg, lon=74.0493*u.deg, height=0*u.m)
        else:
            raise('Unknown Telescope.')
        return telescope_coords

    def get_time(self):
        """
        Return mean time of the observation in mjd (time obj)
        """
        time = np.mean(self.getTimeRange())
        time = Time( time/(24*3600.), format='mjd')
        time.delta_ut1_utc = 0. # no need to download precise table for leap seconds
        return time

    def get_elev(self):
        """
        Return mean elevation
        """
        coord = self.getPhaseCentre(skycoordobj=True)
        return coord.transform_to(AltAz(obstime=self.get_time(),location=self.get_telescope_coords())).alt.deg

    def get_sun_dist(self):
        """
        Return sun distance from the pointing centre in deg
        """
        coord_sun = get_sun(self.get_time())
        coord_sun = SkyCoord(ra=coord_sun.ra,dec=coord_sun.dec) # fix transformation issue
        coord = self.getPhaseCentre(skycoordobj=True)
        return coord.separation(coord_sun).deg
    
    def get_jupiter_dist(self):
        """
        Return Jupiter distance from the pointing centre in deg
        """
        coord_jupiter = get_body('jupiter', self.get_time(), self.get_telescope_coords())
        coord_jupiter = SkyCoord(ra=coord_jupiter.ra,dec=coord_jupiter.dec) # fix transformation issue
        coord = self.getPhaseCentre(skycoordobj=True)
        return coord.separation(coord_jupiter).deg
    
    def get_hour_angle(self):
        """
        Return the hour angle in hrs of the phase centre
        """
        coord = self.getPhaseCentre(skycoordobj=True)
        lst = self.get_time().sidereal_time('mean', self.get_telescope_coords().lon)
        ha = lst - coord.ra # hour angle
        return ha.deg/15.


    def get_hist(self):
        """
        Return the Hystory subtable cell 1 (containing preprocessin ginfo) as a list
        """

        with tables.table(self.pathMS + "/HISTORY", ack=False) as table:
            col = table.col('APP_PARAMS')
            try:
                return col.getcell(1)
            except:
                return []


    def print_ateam_demix(self):
        """
        Some debug logs on the ateam
        """
        hist = self.get_hist()
        for h in hist:
            #print(h)
            if 'demixer.subtractsources' in h:
                logger.debug(f'{self.nameMS}: DEMIX - Sources = {h.split('=')[1]}')
            if 'demixer.ignoretarget' in h:
                logger.debug(f'{self.nameMS}: DEMIX - Ignoretarget = {h.split('=')[1]}')


    def get_ateam_demix(self):
        """
        Return a list of the demixed ateam sources
        """
        hist = self.get_hist()
        demixed = None
        for h in hist:
            if 'demixer.subtractsources' in h:
                demixed = h.split('=')[1].replace('\'','').replace(' ','')
        try:
            return demixed[1:-1].split(',')
        except:
            return []

    def distBrightSource(self, name):
        """
        Get the distance in deg from some bright sources
        """
        ateam={'CygA':{'ra':299.8679167, 'dec':40.7338889},
                'CasA':{'ra':350.8583333, 'dec':58.8000000},
                'TauA':{'ra':83.6333333, 'dec':22.0144444},
                'VirA':{'ra':187.7058333, 'dec':12.3911111},
                '3C338':{'ra':247.160333, 'dec':39.551556},
                '3C380':{'ra':277.382420,'dec':48.746156}
        }
        
        if name not in ateam.keys():
            logger.error('Unknown source for distance: %s' % name)
            logger.error('Use: '+' '.join(ateam.keys()))
            raise

        coord_bright = SkyCoord(ra=ateam[name]['ra']*u.deg, dec=ateam[name]['dec']*u.deg)
        ra, dec = self.getPhaseCentre()
        return coord_bright.separation(SkyCoord(ra*u.deg, dec*u.deg)).deg

    def setPathVariables(self, pathMS):
        """
        Set logistical variables.
        """
        self.pathMS        = pathMS

        indexLastSlash     = self.pathMS.rfind('/')

        if '/' in self.pathMS: self.pathDirectory = self.pathMS[ : indexLastSlash]
        else: self.pathDirectory = './'

        self.nameMS        = self.pathMS[indexLastSlash + 1 : -3]

    def move(self, pathMSNew, overwrite=False, keepOrig=False):
        """
        Move (or rename) the MS to another locus in the file system.
        """
        logger.debug('Move: '+self.pathMS+' -> '+pathMSNew)
        if overwrite == True:
            lib_util.check_rm(pathMSNew)
        if not os.path.exists(pathMSNew):
            if keepOrig:
                shutil.copytree(self.pathMS, pathMSNew)
            else:
                shutil.move(self.pathMS, pathMSNew)

            self.setPathVariables(pathMSNew)


    def setNameField(self, nameField):
        """
        Set field name.
        """
        pathFieldTable = self.pathMS + "/FIELD"
        tables.taql("update $pathFieldTable set NAME=$nameField")

    def setObsID(self, ObsID):
        """
        Set ObsID.
        """
        pathObsIDTable = self.pathMS + "/OBSERVATION"
        tables.taql("update $pathObsIDTable set LOFAR_OBSERVATION_ID=$ObsID")

    def setSpwName(self, SpwName):
        """
        Set Spw Name.
        """
        pathSpwTable = self.pathMS + "/SPECTRAL_WINDOW"
        tables.taql("update $pathSpwTable set NAME=$SpwName")

    def setCode(self, code):
        """
        Set Observation Code.
        """
        pathcodeTable = self.pathMS + "/FIELD"
        tables.taql("update $pathcodeTable set CODE=$code")


    def getNameField(self, checkCalName=False, updateCalName=False):
        """
        Retrieve field name.
        checkCalNAme: if true check if the target is a calibrator although the fieldname says otherwise
        updateCalName: if it turns out to be a calibrator, update the name in FIELD table
        """
        pathFieldTable = self.pathMS + "/FIELD"
        nameField      = (tables.taql("select NAME from $pathFieldTable")).getcol("NAME")[0]

        # If the field name is not a recognised calibrator name, one of two scenarios is true:
        # 1. The field is not a calibrator field;
        # 2. The field is a calibrator field, but the name was not properly set.
        # The following lines correct the field name if scenario 2 is the case.
        if (checkCalName and not self.isCalibrator()):
            calibratorDistanceThreshold = 0.1 # in degrees
            if (self.getCalibratorDistancesSorted()[0] < calibratorDistanceThreshold):
                nameField = self.getCalibratorNamesSorted()[0]
                if updateCalName: self.setNameField(nameField)

        return nameField


    def getCalibratorDistancesSorted(self):
        """
        Returns a list of distances (in degrees) to known calibrators, sorted by distance from small to large.
        """
        myRA, myDec                                    = self.getPhaseCentre()
        calibratorRAs, calibratorDecs, calibratorNames = lib_util.getCalibratorProperties()
        calibratorDistances                            = lib_util.distanceOnSphere(myRA, myDec, calibratorRAs, calibratorDecs)

        calibratorDistancesSorted                      = np.sort(calibratorDistances)
        return calibratorDistancesSorted


    def getCalibratorNamesSorted(self):
        """
        Returns a list of names of known calibrators, sorted by distance from small to large.
        """
        myRA, myDec                                    = self.getPhaseCentre()
        calibratorRAs, calibratorDecs, calibratorNames = lib_util.getCalibratorProperties()
        calibratorDistances                            = lib_util.distanceOnSphere(myRA, myDec, calibratorRAs, calibratorDecs)

        calibratorNamesSorted                          = calibratorNames[np.argsort(calibratorDistances)]
        return calibratorNamesSorted


    def isCalibrator(self):
        """
        Returns whether the field is a known calibrator field or not.
        """
        calibratorRAs, calibratorDecs, calibratorNames = lib_util.getCalibratorProperties()
        return (self.getNameField() in calibratorNames)


    def concretiseString(self, stringOriginal):
        """
        Returns a concretised version of the string 'stringOriginal', with keywords filled in.
        More keywords (which start with '$') and their conversions can be added below.
        """
        stringCurrent = stringOriginal.replace("$pathMS",        self.pathMS)
        stringCurrent = stringCurrent.replace( "$pathDirectory", self.pathDirectory)
        stringCurrent = stringCurrent.replace( "$nameMS",        self.nameMS)
        stringCurrent = stringCurrent.replace( "$nameField",     self.getNameField())

        return stringCurrent


    def getFreqs(self):
        """
        Get chan frequencies in Hz
        """
        with tables.table(self.pathMS + "/SPECTRAL_WINDOW", ack = False) as t:
            freqs = t.getcol("CHAN_FREQ")

        return freqs[0]


    def getNchan(self):
        """
        Find number of channels
        """
        with tables.table(self.pathMS + "/SPECTRAL_WINDOW", ack = False) as t:
            nchan = t.getcol("NUM_CHAN")
        assert (nchan[0] == nchan).all() # all SpWs have same channels?

        #logger.debug("%s: channel number: %i", self.pathMS, nchan[0])
        return nchan[0]


    def getChanband(self):
        """
        Find bandwidth of a channel in Hz
        """
        with tables.table(self.pathMS + "/SPECTRAL_WINDOW", ack = False) as t:
            chan_w = t.getcol("CHAN_WIDTH")[0]
        assert all(x == chan_w[0] for x in chan_w) # all chans have same width

        #logger.debug("%s: channel width (MHz): %f", self.pathMS, chan_w[0] / 1.e6)
        return chan_w[0]


    def getTimeRange(self):
        """
        Return the time interval of this observation
        """
        with tables.table(self.pathMS, ack = False) as t:
            return ( t.getcol("TIME")[0], t.getcol("TIME")[-1] )


    def getNtime(self):
        """
        Returns the number of time slots in this MS
        """
        with tables.table(self.pathMS, ack = False) as t:
            return len(np.unique(t.getcol("TIME")))


    def getTimeInt(self):
        """
        Get time interval in seconds
        """
        with tables.table(self.pathMS, ack = False) as t:
            nTimes = len(set(t.getcol("TIME")))
        t_init, t_end = self.getTimeRange()
        deltaT = (t_end - t_init) / nTimes

        #logger.debug("%s: time interval (seconds): %f", self.pathMS, deltaT)
        return deltaT


    def getPhaseCentre(self, skycoordobj=False):
        """
        Get the phase centre (in degrees) of the first source (is it a problem?) of an MS.
        skycoordobj: if True return a skycoord object, otherwise ra,dec in deg
        """
        field_no = 0
        ant_no   = 0
        with tables.table(self.pathMS + "/FIELD", ack = False) as field_table:
            direction = field_table.getcol("PHASE_DIR")
        RA        = direction[ant_no, field_no, 0]
        Dec       = direction[ant_no, field_no, 1]

        if (RA < 0):
            RA += 2 * np.pi

        if skycoordobj:
            return SkyCoord(np.degrees(RA)*u.deg, np.degrees(Dec)*u.deg)
        else:
            #logger.debug("%s: phase centre (degrees): (%f, %f)", self.pathMS, np.degrees(RA), np.degrees(Dec))
            return (np.degrees(RA), np.degrees(Dec))

    def getTelescope(self):
        """
        Return telescope name such as "LOFAR" or "GMRT"
        """
        with tables.table(self.pathMS+'/OBSERVATION', ack = False) as t:
            return t.getcell("TELESCOPE_NAME",0)

    def getAntennaSet(self):
        """
        If LBA observation, return obs mode: INNER, OUTER, SPARSE_EVEN, SPARSE_ODD
        """
        if self.getTelescope() != 'LOFAR':
            raise("Only LOFAR has Antenna Sets.")

        with tables.table(self.pathMS+'/OBSERVATION', ack = False) as t:
            return t.getcell("LOFAR_ANTENNA_SET",0)

    def getObsID(self):
        """
        Return LOFAR observation ID
        """
        with tables.table(self.pathMS+'/OBSERVATION', ack = False) as t:
            return int(t.getcell("LOFAR_OBSERVATION_ID",0))

    def getFWHM(self, freq='mid', elliptical=False):
        """
        Return the expected FWHM in degree (fwhm_maj, fwhm_min) where maj is N-S and min E-W
        freq: min,max,med - which frequency to use to estimate the beam size
        """
        # get minimum freq as it has the largest FWHM
        if freq == 'min':
            beamfreq = np.min(self.getFreqs()) 
        elif freq == 'max':
            beamfreq = np.max(self.getFreqs()) 
        elif freq == 'mid':
            beamfreq = np.mean(self.getFreqs()) 
        else:
            raise "Wrong freq value for beam size, use: min|mid|max."

        if self.getTelescope() == 'LOFAR':

            # Following numbers are based at 60 MHz (old.astron.nl/radio-observatory/astronomers/lofar-imaging-capabilities-sensitivity/lofar-imaging-capabilities/lofa)
            scale = 60e6/beamfreq

            if 'OUTER' in self.getAntennaSet():
                fwhm = 3.88*scale
            elif 'SPARSE' in self.getAntennaSet():
                fwhm = 4.85*scale
            elif 'INNER' in self.getAntennaSet():
                fwhm = 9.77*scale
            else:
                logger.info('Unknown antenna configuration, assuming SPARSE...')
                fwhm = 4.85*scale #same configuration of SPARSE
            __ra, dec = self.getPhaseCentre()

            if elliptical:
                return np.array([fwhm/np.cos(np.deg2rad(53-dec)), fwhm])
            else:
                return fwhm
                
        elif self.getTelescope() == 'GMRT':
            # equation from http://gmrt.ncra.tifr.res.in/gmrt_hpage/Users/doc/manual/Manual_2013/manual_20Sep2013.pdf    
            return (85.2/60) * (325.e6 / beamfreq)

        else:
            raise('Only LOFAR or GMRT implemented.')

    def getAvgFactors(self, keep_IS):
        """
        Get the time and frequency averaging factor to arrive at the standard LiLF widefield processing resolution
        depending on SPARSE/OUTER, frequency coverage and Dutch/IS.
        keep_IS: compute resolution for IS data
        """
        nchan_per_sb = round(0.195312e6/self.getChanband())
        logger.debug(f'nchan_per_sb: {nchan_per_sb}')
        timeint = self.getTimeInt()
        minfreq = np.min(self.getFreqs())
        if nchan_per_sb == 1:
            avg_factor_f = 1
        # elif nchan % 2 == 0 and MSs.isHBA: # case HBA
        #    avg_factor_f = int(nchan / 4)  # to 2 ch/SB
        elif nchan_per_sb % 8 == 0 and minfreq < 40e6:
            avg_factor_f = int(nchan_per_sb / 8)  # to 8 ch/SB
        elif nchan_per_sb % 8 == 0 and 'SPARSE' in self.getAntennaSet():
            avg_factor_f = int(nchan_per_sb / 8)  # to 8 ch/SB
        elif nchan_per_sb % 4 == 0:
            avg_factor_f = int(nchan_per_sb / 4)  # to 4 ch/SB
        elif nchan_per_sb % 5 == 0:
            avg_factor_f = int(nchan_per_sb / 5)  # to 5 ch/SB
        else:
            logger.error('Channels should be a multiple of 4 or 5.')
            sys.exit(1)

        if keep_IS:
            avg_factor_f = int(nchan_per_sb / 32)  # to have the full FoV in LBA we need 32 ch/SB
        if avg_factor_f < 1: avg_factor_f = 1

        avg_factor_t = int(np.round(2 / timeint)) if keep_IS else int(np.round(4 / timeint))  # to 4 sec (2 for IS)
        if avg_factor_t < 1: avg_factor_t = 1
        return avg_factor_t, avg_factor_f

    def makeBeamReg(self, outfile, pb_cut=None, to_pbval=0.5, freq='mid'):
        """
        Create a ds9 region of the beam to FWHM by default
        outfile : str
            output file
        pb_cut : float, optional
            diameter of the beam in deg
        to_pbval: float, optional
            make to this value of the beam (e.g. 0 to arrive to the null), default is FWHM (0.5)
        freq: min,max,med 
            which frequency to use to estimate the beam size
        """
        logger.debug('Making PB region: '+outfile)
        ra, dec = self.getPhaseCentre()

        if pb_cut is None:
            if to_pbval < 0.1: to_pbval=0.1 # it's a gaussian, not a real beam, so 0 would be inf
            fwhm = self.getFWHM(freq=freq, elliptical=True)
            sigma = fwhm/(2*np.sqrt(2*np.log(2)))
            def inv_gaus(y, sigma, x0=0):
                x = x0 + np.sqrt(-2 * np.log(y) * sigma ** 2)
                return x
            radius = inv_gaus(to_pbval, sigma)
        else:
            radius = np.array([pb_cut/(2.*np.cos(np.deg2rad(53-dec))),pb_cut/2.])

        s = Shape('ellipse', None)
        s.coord_format = 'fk5'
        s.coord_list = [ ra, dec, radius[1], radius[0], 0.0 ] # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '2', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color=red text="beam"'

        regions = pyregion.ShapeList([s])
        lib_util.check_rm(outfile)
        regions.write(outfile)

    def getMaxBL(self, check_flags=True, dutch_only=False, uvw=True):
        """
        Return the max BL length in meters
        dutch_only: bool, check only the dutch stations, default=False
        uvw: bool, get max uvw length, not just uv.
        """
        if dutch_only:
            with tables.taql(f"SELECT FROM {self.pathMS} WHERE ANTENNA1 IN [SELECT ROWID() FROM ::ANTENNA WHERE NAME \
                         ~p/[CR]S*/] && ANTENNA2 in [SELECT ROWID() FROM ::ANTENNA WHERE NAME ~p/[CR]S*/]") as t:
                if check_flags:
                    t = t.query('not all(FLAG)')
                col = t.getcol('UVW')
        else:
            with tables.table(self.pathMS, ack=False) as t:
                if check_flags:
                    t = t.query('not all(FLAG)')
                col = t.getcol('UVW')
        if uvw:
            maxdist = np.nanmax( np.sqrt(col[:,0] ** 2 + col[:,1] ** 2 + col[:,2] ** 2) )
        else:
            maxdist = np.nanmax( np.sqrt(col[:,0] ** 2 + col[:,1] ** 2) )
        return maxdist

    def getResolution(self, check_flags=True):
        """
        Return the expected resolution (in arcsec) of the MS
        Completely flagged lines can be removed
        """
        c = 299792458. # in metres per second

        with tables.table(self.pathMS+'/SPECTRAL_WINDOW', ack = False) as t:
            wavelength = c / np.max(t.getcol("CHAN_FREQ")[0]) # in metres
        
        maxdist = self.getMaxBL(check_flags)

        return float('%.1f'%(wavelength / maxdist * (180 / np.pi) * 3600)) # in arcsec

    def getPixelScale(self, check_flags=True):
        """
        Return a reasonable pixel scale
        """
        res = self.getResolution(check_flags)
        return round(res*2/4) # reasonable value (4" for Dutch LBA)

    def getAntennas(self):
        """
        Return a list of antenna names
        """
        pathAntennaTable = self.pathMS + "/ANTENNA"
        antennas = (tables.taql('select NAME from $pathAntennaTable')).getcol('NAME')
        return antennas

    def isAllFlagged(self):
        """
        Is the dataset fully flagged?
        """
        try:
            with tables.table(self.pathMS, ack = False) as t:
                return np.all(t.getcol('FLAG'))
        except MemoryError: # can happen e.g. for full MS in timesplit (with IS and/or HBA)
            logger.warning('Caugt MemoryError in checking for fully flagged MS! This can happen when working with large '
                           'measurement sets. You might want to manually inspect the flags. Trying to proceed...')
            
    def fractionalFlag(self):
        """
        Return the fraction of flagged data
        """
        with tables.table(self.pathMS, ack = False) as t:
            f = t.getcol('FLAG')
            return np.sum(f)/f.size
        
