#!/usr/bin/python

import os, sys, shutil

from casacore import tables
import numpy as np
import pyregion
from pyregion.parser_helper import Shape
from LiLF import lib_util

from LiLF.lib_log import logger

class AllMSs(object):

    def __init__(self, pathsMS, scheduler, check_flags=True):
        """
        pathsMS:    list of MS paths
        scheduler:  scheduler object
        check_flag: if true ignore fully flagged ms
        """
        self.scheduler = scheduler

        # sort them, useful for some concatenating steps
        if len(pathsMS) == 0:
            raise('Cannot find MS files.')

        self.mssListObj = []
        for pathMS in sorted(pathsMS):
            ms = MS(pathMS)
            if check_flags and ms.isAllFlagged(): 
                logger.warning('Skip fully flagged ms: %s' % pathMS)
            else:
                self.mssListObj.append(MS(pathMS))

        if len(self.mssListObj) == 0:
            raise('ALL MS files flagged.')

        self.mssListStr = [ms.pathMS for ms in self.mssListObj]


    def getListObj(self):
        """
        """
        return self.mssListObj


    def getListStr(self):
        """
        """
        return self.mssListStr

    def getNThreads(self):
        """
        Return the max number of threads in one machine assuming all MSs run at the same time
        """
        if self.scheduler.max_processors < len(self.mssListStr): NThreads = 1
        else:
            NThreads = int(np.rint( self.scheduler.max_processors/len(self.mssListStr) ))

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
        return [item for sublist in freqs for item in sublist] # flatten


    def getBandwidth(self):
        """
        Return the total span of frequency covered by this MS set
        """
        freqs = self.getFreqs()
        return max(freqs) - min(freqs)


    def run(self, command, log, commandType='', maxThreads=None):
        """
        Run command 'command' of type 'commandType', and use 'log' for logger,
        for each MS of AllMSs.
        The command and log file path can be customised for each MS using keywords (see: 'MS.concretiseString()').
        Beware: depending on the value of 'Scheduler.max_threads' (see: lib_util.py), the commands are run in parallel.
        """
        # add max num of threads given the total jobs to run
        # e.g. in a 64 processors machine running on 16 MSs, would result in numthreads=4
        if commandType == 'DPPP': command += ' numthreads='+str(self.getNThreads())

        for MSObject in self.mssListObj:
            commandCurrent = MSObject.concretiseString(command)
            logCurrent     = MSObject.concretiseString(log)

            self.scheduler.add(cmd = commandCurrent, log = logCurrent, commandType = commandType)

            # Provide debug output.
            #lib_util.printLineBold("commandCurrent:")
            #print (commandCurrent)
            #lib_util.printLineBold("logCurrent:")
            #print (logCurrent)

        self.scheduler.run(check = True, maxThreads = maxThreads)

    def plot_HAcov(self, plotname='HAcov.png'):
        """
        Show the coverage in HA
        """
        from astropy.coordinates import get_sun, SkyCoord, EarthLocation, AltAz
        from astropy.time import Time
        from astropy import units as u

        telescope = self.mssListObj[0].getTelescope()
        if telescope == 'LOFAR':
            telescope_coords = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m)
        elif telescope == 'GMRT':
            telescope_coords = EarthLocation(lat=19.0948*u.deg, lon=74.0493*u.deg, height=0*u.m)
        else:
            raise('Unknown Telescope.')
        
        for ms in self.mssListObj:
            time = np.mean(ms.getTimeRange())
            time = Time( time/86400, format='mjd')
            time.delta_ut1_utc = 0. # no need to download precise table for leap seconds
            coord_sun = get_sun(time)
            ra, dec = ms.getPhaseCentre()
            coord = SkyCoord(ra*u.deg, dec*u.deg)
            elev = coord.transform_to(AltAz(obstime=time,location=telescope_coords)).alt
            sun_dist = coord.separation(coord_sun)
            lst = time.sidereal_time('mean', telescope_coords.lon)
            ha = lst - coord.ra # hour angle
            logger.info('%s (%s): Hour angle: %.1f hrs - Elev: %.2f (Sun distance: %.0f)' % (ms.nameMS,time.iso,ha.deg/15.,elev.deg,sun_dist.deg))


class MS(object):

    def __init__(self, pathMS):
        """
        pathMS:        path of the MS, without '/' at the end!
        pathDirectory: path of the parent directory of the MS
        nameMS:        name of the MS, without parent directories and extension (which is assumed to be ".MS" always)
        """
        self.setPathVariables(pathMS)

        # If the field name is not a recognised calibrator name, one of two scenarios is true:
        # 1. The field is not a calibrator field;
        # 2. The field is a calibrator field, but the name was not properly set.
        # The following lines correct the field name if scenario 2 is the case.
        calibratorDistanceThreshold = 0.5 # in degrees
        if (not self.isCalibrator()):
            if (self.getCalibratorDistancesSorted()[0] < calibratorDistanceThreshold):
                nameFieldOld = self.getNameField()
                nameFieldNew = self.getCalibratorNamesSorted()[0]
                #logger.warning("Although the field name '" + nameFieldOld + "' is not recognised as a known calibrator name, " +
                #                "the phase centre coordinates suggest that this scan is a calibrator scan. Changing field name into '" +
                #                nameFieldNew + "'...")
                self.setNameField(nameFieldNew)


    def setPathVariables(self, pathMS):
        """
        Set logistical variables.
        """
        self.pathMS        = pathMS

        indexLastSlash     = self.pathMS.rfind('/')

        self.pathDirectory = self.pathMS[ : indexLastSlash]
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


    def getNameField(self):
        """
        Retrieve field name.
        """
        pathFieldTable = self.pathMS + "/FIELD"
        nameField      = (tables.taql("select NAME from $pathFieldTable")).getcol("NAME")[0]
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

        logger.debug("%s: channel number: %i", self.pathMS, nchan[0])
        return nchan[0]


    def getChanband(self):
        """
        Find bandwidth of a channel in Hz
        """
        with tables.table(self.pathMS + "/SPECTRAL_WINDOW", ack = False) as t:
            chan_w = t.getcol("CHAN_WIDTH")[0]
        assert all(x == chan_w[0] for x in chan_w) # all chans have same width

        logger.debug("%s: channel width (MHz): %f", self.pathMS, chan_w[0] / 1.e6)
        return chan_w[0]


    def getTimeInt(self):
        """
        Get time interval in seconds
        """
        with tables.table(self.pathMS, ack = False) as t:
            nTimes = len(set(t.getcol("TIME")))
        t_init, t_end = self.getTimeRange()
        deltaT = (t_end - t_init) / nTimes

        logger.debug("%s: time interval (seconds): %f", self.pathMS, deltaT)
        return deltaT


    def getTimeRange(self):
        """
        Return the time interval of this observation
        """
        with tables.table(self.pathMS, ack = False) as t:
            return ( t.getcol("TIME")[0], t.getcol("TIME")[-1] )


    def getPhaseCentre(self):
        """
        Get the phase centre (in degrees) of the first source (is it a problem?) of an MS.
        """
        field_no = 0
        ant_no   = 0
        with tables.table(self.pathMS + "/FIELD", ack = False) as field_table:
            direction = field_table.getcol("PHASE_DIR")
        RA        = direction[ant_no, field_no, 0]
        Dec       = direction[ant_no, field_no, 1]

        if (RA < 0):
            RA += 2 * np.pi

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
        
    def getFWHM(self, freq='mid'):
        """
        Return the expected FWHM in degree
        freq: min,max,med - which frequency to use to estimate the beam size
        """
        # get minimum freq as it has the largest FWHM
        if freq == 'min':
            beamfreq = np.min(self.getFreqs()) 
        elif freq == 'max':
            beamfreq = np.max(self.getFreqs()) 
        elif freq == 'mid':
            beamfreq = np.mean(self.getFreqs()) 

        if self.getTelescope() == 'LOFAR':

            # Following numbers are based at 60 MHz (old.astron.nl/radio-observatory/astronomers/lofar-imaging-capabilities-sensitivity/lofar-imaging-capabilities/lofa)
            scale = 60e6/beamfreq 

            if 'OUTER' in self.getAntennaSet():
                return 3.88*scale
            elif 'SPARSE' in self.getAntennaSet():
                return 4.85*scale
            elif 'INNER' in self.getAntennaSet():
                return 9.77*scale
                
        elif self.getTelescope() == 'GMRT':
            # equation from http://gmrt.ncra.tifr.res.in/gmrt_hpage/Users/doc/manual/Manual_2013/manual_20Sep2013.pdf    
            return (85.2/60) * (325.e6 / beamfreq)

        else:
            raise('Only LOFAR or GMRT implemented.')

    def makeBeamReg(self, outfile, pb_cut=None, to_null=False, freq='mid'):
        """
        Create a ds9 region of the beam
        outfile : str
            output file
        pb_cut : float, optional
            diameter of the beam
        to_null : bool, optional
            arrive to the first null, not the FWHM
        freq: min,max,med 
            which frequency to use to estimate the beam size
        """
        logger.debug('Making PB region: '+outfile)
        ra, dec = self.getPhaseCentre()

        if pb_cut is None:
            radius = self.getFWHM(freq=freq)/2.
        else:
            radius = pb_cut/2.

        if to_null: radius *= 2 # rough estimation

        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ ra, dec, radius ] # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '2', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color=red text="beam"'

        regions = pyregion.ShapeList([s])
        lib_util.check_rm(outfile)
        regions.write(outfile)

    def getResolution(self):
        """
        Return the expected resolution (in arcsec) of the MS
        Completely flagged lines are removed
        """
        c = 299792458. # in metres per second

        with tables.table(self.pathMS, ack = False).query('not all(FLAG)') as t:
            col = t.getcol('UVW')

        with tables.table(self.pathMS+'/SPECTRAL_WINDOW', ack = False) as t:
            wavelength = c / t.getcol('REF_FREQUENCY')[0]             # in metres
        #print 'Wavelength:', wavelength,'m (Freq: '+str(t.getcol('REF_FREQUENCY')[0]/1.e6)+' MHz)'

        maxdist = np.nanmax( np.sqrt(col[:,0] ** 2 + col[:,1] ** 2) )

        #return int(round(wavelength / maxdist * (180 / np.pi) * 3600)) # in arcseconds
        return float('%.1f'%(wavelength / maxdist * (180 / np.pi) * 3600)) # in arcsec

    def isAllFlagged(self):
        """
        Is the dataset fully flagged?
        """
        with tables.table(self.pathMS, ack = False) as t:
            return np.all(t.getcol('FLAG'))
