#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Reinout van Weeren, Tammo Jan Dijkema and Andre Offringa
'''

import argparse

from losoto import h5parm

import lib_ms, lib_util


def pipeline_uGMRT_transfer(pathsMS, pathCalibratorH5Parm, pathDirectoryLogs, pathDirectoryParSets = "./parsets"):

    # Initialise parameter set settings.
    nameParSetSolve   = "DPPP_uGMRT_sol_dummy.parset"
    pathParSetSolve   = pathDirectoryParSets + '/' + nameParSetSolve

    # Initialise processing objects.
    scheduler          = lib_util.Scheduler(dry = False, log_dir = pathDirectoryLogs)
    MSs                = lib_ms.AllMSs(pathsMS, scheduler)

    # Add model column, and fill with ones.
    for MSObject in MSs.get_list_obj():
        lib_util.columnAddSimilar(MSObject.pathMS, "MODEL_DATA", "DATA", "TiledMODEL_DATAMartijn",
                                  overwrite = False, fillWithOnes = True, comment = "", verbose = True)

    # Create ParmDBs with dummy values.
    MSs.run(command = "DPPP " + pathParSetSolve + " msin=$pathMS gaincal.parmdb=$pathMS/instrument",
            commandType = "DPPP", log = "transfer_$nameMS.log")

    # Create H5Parm files.
    MSs.run(command = "H5parm_importer.py $pathDirectory/$nameMS.h5 $pathMS", commandType = "python", log = "transfer_$nameMS.log")

    # Open H5Parm files, and fill with bandpass from 'pathCalibratorH5Parm'.
    # 1. Load calibrator data.
    # 2. Fill target H5Parms with bandpass solutions.
    for MSObject in MSs.get_list_obj():
        objectH5Parm = h5parm.h5parm(MSObject.pathDirectory + "/" + MSObject.nameMS + ".h5", readonly = False)
        # Do stuff.
        print ("Under construction!")
        objectH5Parm.close()
    # 3. Save H5Parms.

    # Apply solutions to target fields.



if (__name__ == "__main__"):

    # If the program is run from the command line, parse arguments.
    parser                         = argparse.ArgumentParser(description = "Pipeline step 4: Transfer of solutions.")
    parser.add_argument("pathsMS",              help = "Paths to the MSs to transfer solutions to.")
    parser.add_argument("pathCalibratorH5Parm", help = "Path to the calibrator H5Parm file whose solutions are to be applied.")
    parser.add_argument("pathDirectoryLogs",    help = "Directory containing log files.")
    arguments                      = parser.parse_args()

    # Temporary!
    arguments.pathsMS              = ["/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsTarget/P149.7+03.4/MSs/scanID2.MS",
                                      "/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsTarget/P149.7+03.4/MSs/scanID14.MS"]
    arguments.pathCalibratorH5Parm = "/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsCalibrator/scanID1/solutions/bandpassesTECs.h5"
    arguments.pathDirectoryLogs    = "/disks/strw3/oei/uGMRTCosmosCut-PiLF/logs"

    pipeline_uGMRT_transfer(arguments.pathsMS, arguments.pathCalibratorH5Parm, arguments.pathDirectoryLogs)