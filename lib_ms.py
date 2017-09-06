#!/usr/bin/python

import logging, os, shutil, sys

from casacore import tables
import numpy as np

import lib_util
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
        """
        """
        return self.mss_list_obj


    def get_list_str(self):
        """
        """
        return self.mss_list_str


    def get_str_wsclean(self):
        """
        Return a string with all MS paths, useful for wsclean
        """
        return ' '.join(self.mss_list_str)


    def run(self, command, commandType, log):
        """
        """
        for MSObject in self.mss_list_obj:
            commandCurrent = MSObject.concretiseString(command)
            logCurrent     = MSObject.concretiseString(log)

            self.scheduler.add(cmd = commandCurrent, log = logCurrent, commandType = commandType)

            # Provide debug output.
            lib_util.printLineBold("commandCurrent:")
            print (commandCurrent)
            lib_util.printLineBold("logCurrent:")
            print (logCurrent)

        self.scheduler.run(check = True)


class MS(object):

    def __init__(self, pathMS):
        """
        pathMS:        path to the MS, without '/' at the end!
        pathDirectory: path to the parent directory of the MS
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
                logging.warning("Although the field name '" + nameFieldOld + "' is not recognised as a known calibrator name, " +
                                "the phase centre coordinates suggest that this scan is a calibrator scan. Changing field name into '" +
                                nameFieldNew + "'...")
                self.setNameField(nameFieldNew)


    def setPathVariables(self, pathMS):
        """
        Set logistical variables.
        """
        self.pathMS        = pathMS

        indexLastSlash     = self.pathMS.rfind('/')

        self.pathDirectory = self.pathMS[ : indexLastSlash]
        self.nameMS        = self.pathMS[indexLastSlash + 1 : -3]


    def move(self, pathMSNew):
        """
        Move (or rename) the MS to another locus in the file system.
        """
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
        myRA, myDec                                    = self.get_phase_centre()
        calibratorRAs, calibratorDecs, calibratorNames = lib_util.getCalibratorProperties()
        calibratorDistances                            = lib_util.distanceOnSphere(myRA, myDec, calibratorRAs, calibratorDecs)

        calibratorDistancesSorted                      = np.sort(calibratorDistances)
        return calibratorDistancesSorted


    def getCalibratorNamesSorted(self):
        """
        Returns a list of names of known calibrators, sorted by distance from small to large.
        """
        myRA, myDec                                    = self.get_phase_centre()
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
        stringCurrent = stringOriginal.replace("$pathMS",    self.pathMS)
        stringCurrent = stringCurrent.replace( "$nameMS",    self.nameMS)
        stringCurrent = stringCurrent.replace( "$nameField", self.getNameField())

        return stringCurrent


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
        return chan_w[0]


    def find_timeint(self):
        """
        Get time interval in seconds
        """
        with tables.table(self.pathMS, ack = False) as t:
            nTimes = len(set(t.getcol("TIME")))
        with tables.table(self.pathMS + "/OBSERVATION", ack = False) as t:
            deltat = (t.getcol("TIME_RANGE")[0][1] - t.getcol("TIME_RANGE")[0][0]) / nTimes

        logging.debug("%s: time interval (seconds): %f", self.pathMS, deltat)
        return deltat


    def get_phase_centre(self):
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

        logging.debug("%s: phase centre (degrees): (%f, %f)", self.pathMS, np.degrees(RA), np.degrees(Dec))
        return (np.degrees(RA), np.degrees(Dec))