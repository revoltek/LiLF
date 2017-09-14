#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Reinout van Weeren, Tammo Jan Dijkema and Andre Offringa
'''

import argparse

import lib_ms, lib_util


def pipeline_uGMRT_transfer(pathsMS, pathCalibratorH5Parm, pathDirectoryLogs, pathDirectoryParSets = "./parsets"):

    # Initialise parameter set settings.
    nameParSetSolve   = "DPPP_uGMRT_sol_dummy.parset"
    pathParSetSolve   = pathDirectoryParSets + '/' + nameParSetSolve

    # Initialise processing objects.
    scheduler          = lib_util.Scheduler(dry = False, log_dir = pathDirectoryLogs)
    MSs                = lib_ms.AllMSs(pathsMS, scheduler)

    # Add model column, and fill with ones. In this way, we avoid a predict.
    for MSObject in MSs.get_list_obj():
        lib_util.columnAddSimilar(MSObject.pathMS, "MODEL_DATA", "DATA", "TiledMODEL_DATAMartijn",
                                  overwrite = False, fillWithOnes = True, comment = "", verbose = True)

    # Create ParmDBs with meaningless values.
    MSs.run(command = "DPPP " + pathParSetSolve + " msin=$pathMS gaincal.parmdb=$pathMS/instrument",
            commandType = "DPPP", log = "transfer_$nameMS.log")

    # Create H5Parm files.
    MSs.run(command = "H5parm_importer.py $pathDirectory/$nameMS.h5 $pathMS", commandType = "python", log = "transfer_$nameMS.log")

    # Open H5Parm files, and fill with bandpass from 'pathCalibratorH5Parm'.
    # 1. Load calibrator data.
    # 2. Fill target H5Parms with bandpass solutions.
    for MSObject in MSs.get_list_obj():
        # Do stuff.
    # 3. Save H5Parms.

    # Apply solutions to target fields.



if (__name__ == "__main__"):

    pipeline_uGMRT_transfer()