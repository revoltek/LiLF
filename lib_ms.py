#!/usr/bin/python

import os, sys, logging

from casacore import tables
import numpy as np


#logger = logging.getLogger("PiLL")


class AllMSs(object):

    def __init__(self, pathsMS, scheduler):
        """
        pathsMS:   list of MS paths
        scheduler: scheduler object
        """
        self.scheduler    = scheduler

        self.mss_list_str = sorted(pathsMS)

        self.mss_list_obj = []
        for pathMS in self.mss_list_str:
            self.mss_list_obj.append(MS(pathMS))


    def get_list_obj(self):
        return self.mss_list_obj


    def get_list_str(self):
        return self.mss_list_str


    def get_str_wsclean(self):
        """
        Return a string with all MS paths, useful for wsclean
        """
        return ' '.join(self.mss_list_str)


    def run(self, command, commandType, log):
        for pathMS in self.mss_list_str:
            commandCurrent = command.replace("$ms", pathMS)
            logCurrent     = log.replace("$ms", pathMS)
            self.scheduler.add(commandCurrent, logCurrent, commandType)

        self.scheduler.run(check = True)


class MS(object):

    def __init__(self, pathMS):
        """
        pathMS:        path to the MS
        pathDirectory: path to the parent directory of the MS
        nameMS:        name of the MS, without parent directories and extension
        """
        self.pathMS        = pathMS

        indexLastSlash     = self.pathMS.rfind('/')

        self.pathDirectory = self.pathMS[ : indexLastSlash]
        self.nameMS        = self.pathMS[indexLastSlash + 1 : -3]


    def getName(self):
        """
        Retrieve source name.
        """
        pathFieldTable = self.pathMS + "/FIELD"
        name           = (tables.taql("select NAME from $pathFieldTable")).getcol("NAME")[0]
        return name


    def isCalibrator(self):
        """
        Returns whether the source is a calibrator or not.
        """
        return (self.getName() in ["CygA", "3C48", "3C147", "3C196", "3C286", "3C295", "3C380"]) # This list should be expanded to include all potential calibrators.


    def find_nchan(self):
        """
        Find number of channels
        """
        with tables.table(self.pathMS + "/SPECTRAL_WINDOW", ack = False) as t:
            nchan = t.getcol("NUM_CHAN")
        assert (nchan[0] == nchan).all() # all SpWs have same channels?

        logging.debug("%s: channel number (1): %i", self.pathMS, nchan[0])
        #logger.debug("%s: channel number (1): %i", self.pathMS, nchan[0])
        return nchan[0]


    def find_chanband(self):
        """
        Find bandwidth of a channel in Hz
        """
        with tables.table(self.pathMS + "/SPECTRAL_WINDOW", ack = False) as t:
            chan_w = t.getcol("CHAN_WIDTH")[0]
        assert all(x == chan_w[0] for x in chan_w) # all chans have same width

        logging.debug("%s: channel width (MHz): %f", self.pathMS, chan_w[0] / 1.e6)
        #logger.debug("%s: channel width (MHz): %f", self.pathMS, chan_w[0] / 1.e6)
        return chan_w[0]


    def find_timeint(self):
        """
        Get time interval in seconds
        """
        with tables.table(self.pathMS, ack = False) as t:
            nTimes = len(set(t.getcol("TIME")))
        with tables.table(self.pathMS + "/OBSERVATION", ack = False) as t:
            deltat = (t.getcol("TIME_RANGE")[0][1] - t.getcol("TIME_RANGE")[0][0]) / nTimes

        logging.debug("%s: time interval (s): %f", self.pathMS, deltat)
        #logger.debug("%s: time interval (s): %f", self.pathMS, deltat)
        return deltat


    def get_phase_centre(self):
        """
        Get the phase centre of the first source (is it a problem?) of an MS
        values in deg
        """
        field_no = 0
        ant_no   = 0
        with tables.table(self.pathMS + "/FIELD", ack = False) as field_table:
            direction = field_table.getcol("PHASE_DIR")
            ra        = direction[ant_no, field_no, 0]
            dec       = direction[ant_no, field_no, 1]

        if (ra < 0):
            ra += 2 * np.pi

        logging.debug("%s: phase centre (deg): %f - %f", self.pathMS, np.degrees(ra), np.degrees(dec))
        #logger.debug("%s: phase centre (deg): %f - %f", self.pathMS, np.degrees(ra), np.degrees(dec))
        return (np.degrees(ra), np.degrees(dec))


    '''
    import lib_util

    def get_calname(self):
        """
        Check if MS is of a calibrator and return patch name
        """
        RA, Dec                 = self.get_phase_centre()

        # The following list should be expanded to include all calibrators that are possibly used.
        calibratorRAs           = np.array([24.4220808, 85.6505746, 277.3824204, 212.835495, 123.4001379, 299.8681525, 202.784479167]) # in degrees
        calibratorDecs          = np.array([33.1597594, 49.8520094, 48.7461556,  52.202770,  48.2173778,  40.7339156,  30.509088])     # in degrees
        calibratorNames         = np.array(["3C48",     "3C147",    "3C380",     "3C295",    "3C196",     "CygA",      "3C286"])

        distances               = lib_util.distanceOnSphere(RA, Dec, calibratorRAs, calibratorDecs)

        distanceCutoff          = 1                                                                                                       # in degrees

        calibratorIsNear        = (distances < distanceCutoff)
        numberOfCalibratorsNear = np.sum(calibratorIsNear)

        if (numberOfCalibratorsNear == 0):
            logger.info("Error: unknown calibrator.")
            sys.exit()
        elif (numberOfCalibratorsNear == 1):
            calibratorNameCurrent = calibratorNames[calibratorIsNear][0]
            logger.info("Calibrator found: %s." % calibratorNameCurrent)
            return calibratorNameCurrent
        else:
            logger.info("Error: multiple calibrators were nearby.")
            sys.exit()
    '''