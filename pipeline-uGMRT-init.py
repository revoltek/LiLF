#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Tammo Jan Dijkema & Reinout van Weeren

This pipeline chunk
1. Downloads uGMRT data.
2. Converts downloaded data into MS format.
3. Sets up a series of directories.
4. Splits up the MS on a scan-by-scan basis and puts the MSs thus created in the right directories.
5. Modifies the MSs to LOFAR format.
'''

import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *

from casacore import tables
import shutil


def printLineBold(line):
    '''
    Print a line of text 'line' using boldface.
    '''

    boldStart = "\033[1m"
    boldEnd   = "\033[0;0m"
    print (boldStart + line + boldEnd)


def getScanIDs(pathMS, verbose = False):
    '''
    Return the different scan IDs in a MS.
    '''

    scanIDs = (tables.taql("select distinct SCAN_NUMBER from $pathMS")).getcol("SCAN_NUMBER")

    if (verbose):
        printLineBold("scanIDs:")
        print (scanIDs)

    return scanIDs



# 1. Download data from NCRA server



# 2. Convert into MS
# Temporary:
pathDirectoryMS = "/disks/strw3/oei/"
nameMS          = "uGMRTCosmosCut.ms"
pathMS          = pathDirectoryMS + nameMS

if (not os.path.isdir(pathMS)):
    print ("'" + pathMS + "' is not a valid path.")
    sys.exit()



# 3. Set-up directory structure
pathDirectoryMain = pathDirectoryMS + nameMS[ : -3] + "-PiLF/"

if (os.path.isdir(pathDirectoryMain)):
    shutil.rmtree(pathDirectoryMain)
os.mkdir(pathDirectoryMain)
os.mkdir(pathDirectoryMain + "fieldsTarget/")
os.mkdir(pathDirectoryMain + "fieldsCalibrator/")



# 4. Split up the MS, and place the fragments in the right position.
scanIDs       = getScanIDs(pathMS)
numberOfScans = len(scanIDs)

printLineBold("pathMS:")
print (pathMS)

printLineBold("- Split-up of original MS '" + pathMS + "' -")
for scanID, i in zip(scanIDs, np.arange(numberOfScans) + 1):
    pathMSNew = pathDirectoryMain + "scanID" + str(scanID) + ".MS"
    tables.taql("SELECT from $pathMS where SCAN_NUMBER = $scanID giving $pathMSNew as plain")

    printLineBold(str(i) + " / " + str(numberOfScans) + ". Created MS '" + pathMSNew + "'.")