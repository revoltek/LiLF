#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Reinout van Weeren, Tammo Jan Dijkema and Andre Offringa

This pipeline chunk
1. Downloads uGMRT data.
2. Converts downloaded data into MS format.
3. Sets up a series of directories.
4. Splits up the MS on a scan-by-scan basis.
5. Modifies the MSs to LOFAR format.
6. Puts the MSs in the right place in the file tree.

Notes:
Paths to directories do not end with a '/'.
'''

import argparse, logging, os, shutil, sys

from casacore import tables
import numpy as np

import lib_ms, lib_util



def getScanIDs(pathMS, verbose = False):
    '''
    Return the different scan IDs in a MS.
    '''

    scanIDs = (tables.taql("select distinct SCAN_NUMBER from $pathMS")).getcol("SCAN_NUMBER")

    if (verbose):
        lib_util.printLineBold("scanIDs:")
        print (scanIDs)

    return scanIDs


def pipeline_uGMRT_init(pathDirectoryMS, nameMS, verbose = False):

    # 1. Download data from NCRA server




    # 2. Convert into MS




    # 3. Set-up directory structure
    pathMS                        = pathDirectoryMS + '/' + nameMS

    if (not os.path.isdir(pathMS)):
        print ("'" + pathMS + "' is not a valid path.")
        sys.exit()

    pathDirectoryMain             = pathDirectoryMS + '/' + nameMS[ : -3] + "-PiLF"
    pathDirectoryFieldsTarget     = pathDirectoryMain + "/fieldsTarget"
    pathDirectoryFieldsCalibrator = pathDirectoryMain + "/fieldsCalibrator"
    pathDirectoryLogs             = pathDirectoryMain + "/logs"

    if (os.path.isdir(pathDirectoryMain)):
        shutil.rmtree(pathDirectoryMain)
    os.mkdir(pathDirectoryMain)
    os.mkdir(pathDirectoryFieldsTarget)
    os.mkdir(pathDirectoryFieldsCalibrator)
    os.mkdir(pathDirectoryLogs)

    # Initialise logging settings.
    nameFileLog        = "pipeline_uGMRT_init.log"
    pathFileLog        = pathDirectoryLogs + '/' + nameFileLog

    # Initialise logging.
    lib_util.printLineBold("Starting log at '" + pathFileLog + "'...")
    logging.basicConfig(filename = pathFileLog, level = logging.DEBUG)
    logging.info("Started 'pipeline_uGMRT_init.py'!")



    # 4. Split up the MS.
    columnNameVisibilitiesNew = "DATA"
    columnNameFlagsNew        = "FLAG"
    columnNameWeightsNew      = "WEIGHT_SPECTRUM"

    scanIDs                   = getScanIDs(pathMS)
    numberOfScans             = len(scanIDs)


    logging.debug("pathMS:")
    logging.debug(pathMS)

    # Create temporary paths for the sub-MSs.
    pathsMSNew    = np.empty(len(scanIDs), dtype = object)
    i             = 0
    for scanID in scanIDs:
        pathMSNew      = pathDirectoryMain + "/scanID" + str(scanID) + ".MS"
        pathsMSNew[i]  = pathMSNew
        i             += 1
    logging.debug("pathsMSNew:")
    logging.debug(pathsMSNew)


    logging.info("- Split-up of original MS '" + pathMS + "' -")
    for scanID, pathMSNew, i in zip(scanIDs, pathsMSNew, np.arange(numberOfScans) + 1):
        tables.taql("SELECT from $pathMS where SCAN_NUMBER = $scanID giving $pathMSNew as plain")

        logging.info(str(i) + " / " + str(numberOfScans) + ". Created MS '" + pathMSNew + "'.")


    # 6. Move the sub-MSs from a temporary place to their right locus in the file tree.
    for pathMSNew in pathsMSNew:

        MSObject    = lib_ms.MS(pathMSNew)

        if (MSObject.isCalibrator()):
            pathDirectoryCalibrator = pathDirectoryFieldsCalibrator + '/' + MSObject.nameMS
            pathMSFinal             = pathDirectoryCalibrator + '/' + MSObject.nameMS + ".MS"
            if (not os.path.isdir(pathDirectoryCalibrator)):
                os.mkdir(pathDirectoryCalibrator)
                os.mkdir(pathDirectoryCalibrator + "/plots")
                os.mkdir(pathDirectoryCalibrator + "/solutions")
        else:
            pathDirectoryTarget = pathDirectoryFieldsTarget + '/' + MSObject.getNameField()
            pathMSFinal         = pathDirectoryTarget + "/MSs/" + MSObject.nameMS + ".MS"
            if (not os.path.isdir(pathDirectoryTarget)):
                os.mkdir(pathDirectoryTarget)
                os.mkdir(pathDirectoryTarget + "/plots")
                os.mkdir(pathDirectoryTarget + "/MSs")
                os.mkdir(pathDirectoryTarget + "/images")

        logging.info("Moving MS at '" + pathMSNew + "' to '" + pathMSFinal + "'...")

        MSObject.move(pathMSFinal)



if (__name__ == "__main__"):

    # If the program is run from the command line, parse arguments.
    parser                    = argparse.ArgumentParser(description = "Pipeline step 1: Retrieval and reformatting of uGMRT data.")
    parser.add_argument("pathDirectoryMS",                  help = "Path of the directory containing target MS.")
    parser.add_argument("nameMS",                           help = "Name of the target MS, including extension.")
    parser.add_argument("-v", "--verbose", default = False, help = "Whether or not to provide extensive textual diagnostic output. Default: False")
    arguments                 = parser.parse_args()

    # Temporary!
    arguments.pathDirectoryMS = "/disks/strw3/oei"
    arguments.nameMS          = "uGMRTCosmosCut.ms"

    lib_util.printLineBold("Parameters to use:")
    print (arguments)

    # Run the program with appropriate input.
    pipeline_uGMRT_init(arguments.pathDirectoryMS, arguments.nameMS, arguments.verbose)
