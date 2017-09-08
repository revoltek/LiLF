#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Reinout van Weeren, Tammo Jan Dijkema and Andre Offringa

This pipeline chunk works exclusively on calibrator scans, and
1. Creates model data.
2. Creates complex gains.
3. Generates amplitude bandpasses.
4. Generates phase bandpasses (which are more accurate than single delay parameters) and TEC solutions.

Notes:
Paths to directories do not end with a '/'.
'''

import argparse, logging, os

import numpy as np

import lib_ms, lib_util


def pipeline_uGMRT_bandpass(pathsMS, pathDirectoryLogs, pathDirectoryParSets = "./parsets", verbose = False):

    # Initialise parameter set settings.
    nameParSetPredict = "DPPP_uGMRT_predict.parset"
    nameParSetSolve   = "DPPP_uGMRT_sol.parset"
    pathParSetPredict = pathDirectoryParSets + '/' + nameParSetPredict
    pathParSetSolve   = pathDirectoryParSets + '/' + nameParSetSolve

    # Initialise logging settings.
    nameFileLog        = "pipeline_uGMRT_bandpass.log"
    pathFileLog        = pathDirectoryLogs + '/' + nameFileLog

    # Initialise logging.
    lib_util.printLineBold("Starting log at '" + pathFileLog + "'...")
    logging.basicConfig(filename = pathFileLog, level = logging.DEBUG)
    logging.info("Started 'pipeline_uGMRT_bandpass.py'!")

    # Initialise processing objects.
    scheduler          = lib_util.Scheduler(dry = False, log_dir = pathDirectoryLogs)
    MSs                = lib_ms.AllMSs(pathsMS, scheduler)


    # Add model data column (for predict), and corrected data column (for gaincal).
    for MSObject in MSs.get_list_obj():

        lib_util.columnAddSimilar(MSObject.pathMS, "MODEL_DATA",     "DATA", "TiledMODEL_DATAMartijn",
                                  overwrite = False, fillWithZeros = True, comment = "", verbose = True)
        lib_util.columnAddSimilar(MSObject.pathMS, "CORRECTED_DATA", "DATA", "TiledCORRECTED_DATAMartijn",
                                  overwrite = False, fillWithZeros = True, comment = "", verbose = True)

        # Test functionality of class MS.
        print (MSObject.find_nchan())
        print (MSObject.find_chanband())
        print (MSObject.pathDirectory)
        print (MSObject.nameMS)

    '''
    # Set model data column. Instead of predicting 'on the fly' whilst calculating gains, we predict and store in MODEL_DATA.
    # This is a disk space versus computing time trade-off.
    logging.info("Predicting calibrator data...")
    sourceDB = "./models/calib-simple.skydb"
    MSs.run(command = "DPPP " + pathParSetPredict + " msin=$pathMS predict.sourcedb=" + sourceDB + " predict.sources=$nameField",
            commandType = "DPPP", log = "bandpass_$nameMS.log")


    # Calculate complex gains and store in ParmDB format.
    logging.info("Calculating complex gains...")
    for pathMS in MSs.get_list_str():
        print (pathMS)
        lib_util.check_rm(pathMS + "/instrument")
    MSs.run(command = "DPPP " + pathParSetSolve + " msin=$pathMS gaincal.parmdb=$pathMS/instrument",
            commandType = "DPPP", log = "bandpass_$nameMS.log")


    # As long as the transition from ParmDB to H5Parm is incomplete, the following conversion step remains.
    logging.info("Converting ParmDB to H5Parm...")

    # Delete H5Parm files if already existing.
    for MSObject in MSs.get_list_obj():
        pathH5Parm = MSObject.pathDirectory + "/solutions/gainsRaw.h5"
        if (os.path.isfile(pathH5Parm)):
            logging.info("Removing old H5Parm file at '" + pathH5Parm + "'...")
            os.remove(pathH5Parm)

    # Create H5Parm files.
    MSs.run(command = "H5parm_importer.py $pathDirectory/solutions/gainsRaw.h5 $pathMS", commandType = "python", log = "bandpass_$nameMS.log")
    '''

    # Determine and store amplitude and phase bandpass (as well as calibrator TEC solutions).
    logging.info("Calculating amplitude bandpass, phase bandpass and calibrator TEC solutions...")

    MSs.run(command = "dedicated_uGMRT_bandpass.py $pathDirectory/solutions/gainsRaw.h5", commandType = "python", log = "bandpass_$nameMS.log")


if (__name__ == "__main__"):

    # If the program is run from the command line, parse arguments.
    parser                      = argparse.ArgumentParser(description = "Pipeline step 3: Generation of bandpasses.")
    parser.add_argument("pathsMS", help = "Paths to the MSs to act upon.")
    parser.add_argument("pathDirectoryLogs", help = "Directory containing log files.")
    parser.add_argument("-p", "--pathDirectoryParSets", default = "./parsets", help = "Directory containing parameter sets.")
    parser.add_argument("-v", "--verbose", default = False, help = "Whether or not to provide extensive textual diagnostic output. Default: False")
    arguments                   = parser.parse_args()

    # Temporary!
    arguments.pathsMS           = ["/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsCalibrator/scanID1/scanID1.MS"]
    arguments.pathDirectoryLogs =  "/disks/strw3/oei/uGMRTCosmosCut-PiLF/logs"
    arguments.verbose           = True


    lib_util.printLineBold("Parameters to use:")
    print (arguments)

    # Run the program with appropriate input.
    pipeline_uGMRT_bandpass(arguments.pathsMS, arguments.pathDirectoryLogs, arguments.pathDirectoryParSets, arguments.verbose)