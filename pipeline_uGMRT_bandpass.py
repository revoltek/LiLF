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

import argparse, logging

import lib_ms, lib_util


def pipeline_uGMRT_bandpass(pathsMS, pathDirectoryLogs, pathDirectoryParSets = "./parsets"):

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

    # TEMPORARY:
    for MSObject in MSs.get_list_obj():
        print (MSObject)
        print (type(MSObject))
        print (MSObject.isCalibrator())
        print (MSObject.getNameField())

    # Set model data column. Instead of predicting 'on the fly' whilst calculating gains, we predict and store in MODEL_DATA.
    # This is a disk space versus computing time trade-off.
    logging.info("Predicting calibrator data...")
    sourcedb = "./models/calib-simple.skydb"
    MSs.run("DPPP " + pathParSetPredict + " msin=$pathMS predict.sourcedb=" + sourcedb + " predict.sources=$nameField", log = "bandpass_$nameMS.log", commandType = "DPPP")


    # Calculate complex gains and store in ParmDB format.
    logging.info("Calculating complex gains...")
    for pathMS in MSs.get_list_str():
        print (pathMS)
        lib_util.check_rm(pathMS + "/instrument")
    MSs.run("DPPP " + pathParSetSolve + " msin=$pathMS gaincal.parmdb=$pathMS/instrument", log = "bandpass_$nameMS.log", commandType = "DPPP")


    # As long as the transition from ParmDB to H5Parm is incomplete, the following conversion step remains.
    logging.info("Converting ParmDB to H5Parm...")
    MSs.run("H5parm_importer.py $nameMS.h5 $pathMS", log = "bandpass_$nameMS.log", commandType = "python")


    # Determine and store amplitude and phase bandpass (as well as calibrator TEC solutions).
    logging.info("Calculating amplitude bandpass, phase bandpass and calibrator TEC solutions...")
    MSs.run("dedicated_uGMRT_bandpass.py $nameMS.h5", log = "bandpass_$nameMS.log", commandType = "python")


if (__name__ == "__main__"):

    # If the program is run from the command line, parse arguments.
    parser                      = argparse.ArgumentParser(description = "Pipeline step 3: Generation of bandpasses.")
    parser.add_argument("pathsMS", help = "Paths to the MSs to act upon.")
    parser.add_argument("pathDirectoryLogs", help = "Directory containing log files.")
    parser.add_argument("-p", "--pathDirectoryParSets", default = "./parsets", help = "Directory containing parameter sets.")
    arguments                   = parser.parse_args()

    # Temporary!
    arguments.pathsMS           = ["/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsCalibrator/scanID1/scanID1.MS"]
    arguments.pathDirectoryLogs =  "/disks/strw3/oei/uGMRTCosmosCut-PiLF/logs"


    lib_util.printLineBold("Parameters to use:")
    print (arguments)

    # Run the program with appropriate input.
    pipeline_uGMRT_bandpass(arguments.pathsMS, arguments.pathDirectoryLogs, arguments.pathDirectoryParSets)