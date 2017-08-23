#!/usr/bin/python

import os, sys
import numpy as np

import lib_util

import pyrap.tables as tb
from casacore import tables

import logging
logger = logging.getLogger('PiLL')

class AllMss(object):

    def __init__(self, mss, s):
        """
        mss : list of MS filenames
        s : scheduler obj
        """
        self.s = s
        self.mss_list_str = sorted(mss)
        self.mss_list_obj = []
        for ms in sorted(mss):
            mss_list_obj.append(Ms(ms))

    def get_list_obj(self):
        return self.mss_list_obj

    def get_list_str(self):
        return self.mss_liat_str

    def get_str_wsclean(self):
        """
        Return a string with all mss names,
        useful for wsclean
        """
        return ' '.join(self.mss_liat_str)

    def run(self, cmd, log, cmd_type):
        for ms in self.mss_list_str:
            cmd = cmd.replace('$ms', ms)
            log = log.replace('$ms', ms)
            self.s.add(cmd, log, cmd_type)
        self.s.run(check = True)

    
class Ms(object):

    def __init__(self, filename):
        self.ms = filename
        
    def getName(self):
        """
        Retrieve field name.
        """
        pathFieldTable = self.ms + "/FIELD"
        name           = (tables.taql("select NAME from $pathFieldTable")).getcol("NAME")[0]
        return name
    
    def isCalibrator(self):
        """
        Returns whether the source is a calibrator or not.
        """
        if (self.getName() in ["CygA", "3C48", "3C147", "3C196", "3C286", "3C295", "3C380"]):
            return True
        else:
            return False
    
    '''
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
    
    def find_nchan(self):
        """
        Find number of channels
        """
        with tb.table(self.ms + '/SPECTRAL_WINDOW', ack = False) as t:
            nchan = t.getcol('NUM_CHAN')
        assert (nchan[0] == nchan).all() # all spw have same channels?
        logger.debug('%s: Number of channels: %i' (self.ms, nchan[0]))
        return nchan[0]
    
    
    def find_chanband(self):
        """+
        Find bandwidth of a channel in Hz
        """
        with tb.table(self.ms + '/SPECTRAL_WINDOW', ack = False) as t:
            chan_w = t.getcol('CHAN_WIDTH')[0]
        assert all(x == chan_w[0] for x in chan_w) # all chans have same width
        logger.debug('%s: Chan-width: %f MHz' (self.ms, chan_w[0] / 1.e6))
        return chan_w[0]
    
    
    def find_timeint(self):
        """
        Get time interval in seconds
        """
        with tb.table(self.ms, ack = False) as t:
            Ntimes = len(set(t.getcol('TIME')))
        with tb.table(self.ms + '/OBSERVATION', ack = False) as t:
            deltat = (t.getcol('TIME_RANGE')[0][1] - t.getcol('TIME_RANGE')[0][0]) / Ntimes
        logger.debug('%s: Time interval: %f s' (self.ms, deltat))
        return deltat
    
    
    def get_phase_centre(self):
        """
        Get the phase centre of the first source (is it a problem?) of an MS
        values in deg
        """
        field_no = 0
        ant_no   = 0
        with tb.table(self.ms + "/FIELD", ack = False) as field_table:
            direction = field_table.getcol("PHASE_DIR")
            ra        = direction[ ant_no, field_no, 0 ]
            dec       = direction[ ant_no, field_no, 1 ]
        #logger.debug("%s: Phase centre: %f deg - %f deg" (self.ms, ra*180/np.pi, dec*180/np.pi))
        if (ra < 0):
            ra += 2 * np.pi
        return (np.degrees(ra), np.degrees(dec))

'''
def find_nchan(ms):
    """
    Find number of channel in this ms
    """
    with tb.table(ms+'/SPECTRAL_WINDOW', ack=False) as t:
        nchan = t.getcol('NUM_CHAN')
    assert (nchan[0] == nchan).all() # all spw have same channels?
    logger.debug('Channel in '+ms+': '+str(nchan[0]))
    return nchan[0]


def find_chanband(ms):
    """
    Find bandwidth of a channel
    """
    with tb.table(ms+'/SPECTRAL_WINDOW', ack=False) as t:
        chan_w = t.getcol('CHAN_WIDTH')[0]
    assert all(x==chan_w[0] for x in chan_w) # all chans have same width
    logger.debug('Channel width in '+ms+': '+str(chan_w[0]/1e6)+' MHz')
    return chan_w[0]


def find_timeint(ms):
    """
    Get time interval in seconds
    """
    with tb.table(ms, ack=False) as t:
        Ntimes = len(set(t.getcol('TIME')))
    with tb.table(ms+'/OBSERVATION', ack=False) as t:
        deltat = (t.getcol('TIME_RANGE')[0][1]-t.getcol('TIME_RANGE')[0][0])/Ntimes
    logger.debug('Time interval for '+ms+': '+str(deltat))
    return deltat


def get_phase_centre(ms):
    """
    Get the phase centre of the first source (is it a problem?) of an MS
    values in deg
    """
    field_no = 0
    ant_no = 0
    with tb.table(ms + '/FIELD', ack=False) as field_table:
        direction = field_table.getcol('PHASE_DIR')
        ra = direction[ ant_no, field_no, 0 ]
        dec = direction[ ant_no, field_no, 1 ]
    return (ra*180/np.pi, dec*180/np.pi)
'''