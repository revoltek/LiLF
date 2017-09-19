#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Reinout van Weeren, Tammo Jan Dijkema and Andre Offringa
'''

import argparse

import numpy as np
from losoto import h5parm

import lib_ms, lib_util


def pipeline_uGMRT_transfer(pathsMS, pathCalibratorH5Parm, pathDirectoryLogs, pathDirectoryParSets = "./parsets"):

    # Initialise parameter set settings.
    nameParSetSolve                 = "DPPP_uGMRT_sol_dummy.parset"
    nameParSetApply                 = "DPPP_uGMRT_apply.parset"
    pathParSetSolve                 = pathDirectoryParSets + "/" + nameParSetSolve
    pathParSetApply                 = pathDirectoryParSets + "/" + nameParSetApply

    # Initialise processing objects.
    scheduler                       = lib_util.Scheduler(dry = False, log_dir = pathDirectoryLogs)
    MSs                             = lib_ms.AllMSs(pathsMS, scheduler)

    # Add model column, and fill with ones.
    for MSObject in MSs.get_list_obj():
        lib_util.columnAddSimilar(MSObject.pathMS, "MODEL_DATA", "DATA", "TiledMODEL_DATAMartijn",
                                  overwrite = False, fillWithOnes = True, comment = "", verbose = True)

    # Create ParmDBs with dummy values.
    #MSs.run(command = "DPPP " + pathParSetSolve + " msin=$pathMS gaincal.parmdb=$pathMS/instrument",
    #        commandType = "DPPP", log = "transfer_$nameMS.log")

    # Create H5Parm files.
    #MSs.run(command = "H5parm_importer.py $pathDirectory/$nameMS.h5 $pathMS", commandType = "python", log = "transfer_$nameMS.log")

    # Open H5Parm files, and fill with bandpass from 'pathCalibratorH5Parm'.
    # 1. Load calibrator data.
    objectH5Parm                    = h5parm.h5parm(pathCalibratorH5Parm, readonly = True)
    objectSolSet                    = objectH5Parm.getSolset("sol000")
    objectSolTabBandpassesAmplitude = objectSolSet.getSoltab("bandpassAmplitude")
    objectSolTabBandpassesPhase     = objectSolSet.getSoltab("bandpassPhase")
    bandpassesAmplitude             = objectSolTabBandpassesAmplitude.getValues(retAxesVals = False, weight = False)
    bandpassesPhase                 = objectSolTabBandpassesPhase.getValues(    retAxesVals = False, weight = False)
    objectH5Parm.close()

    # Reshape (2, 30, 2048) array to (2, 1, 30, 2048, 24)... How? Using numpy.tile? numpy.broadcast_to? numpy.reshape...?
    print (bandpassesAmplitude.shape)
    bandpassesAmplitudeReshaped     = np.expand_dims(bandpassesAmplitude,         axis = 1)
    bandpassesAmplitudeReshaped     = np.expand_dims(bandpassesAmplitudeReshaped, axis = 4)

    bandpassesPhaseReshaped         = np.expand_dims(bandpassesPhase,             axis = 1)
    bandpassesPhaseReshaped         = np.expand_dims(bandpassesPhaseReshaped,     axis = 4)

    bandpassesAmplitudeTiled = np.tile(bandpassesAmplitudeReshaped, (1, 1, 1, 1, 24))
    bandpassesPhaseTiled     = np.tile(bandpassesPhaseReshaped,     (1, 1, 1, 1, 24))
    print (bandpassesAmplitudeTiled.shape)
    print (bandpassesPhaseTiled.shape)
    from matplotlib import pyplot
    from matplotlib import cm
    for ant in [0, 1]:
        pyplot.imshow(bandpassesAmplitudeTiled[0, 0, ant, :, :], interpolation = "none", aspect = "auto")
        pyplot.savefig("/disks/strw3/oei/uGMRTCosmosCut-PiLF/test" + str(ant) + ".pdf")
        pyplot.close()

        pyplot.imshow(bandpassesPhaseTiled[0, 0, ant, :, :], interpolation = "none", cmap = cm.hsv, vmin = -180, vmax = 180)
        pyplot.savefig("/disks/strw3/oei/uGMRTCosmosCut-PiLF/testFase" + str(ant) + ".pdf")
        pyplot.close()



    # 2. Fill target H5Parms with bandpass solutions.
    for MSObject in MSs.get_list_obj():
        objectH5Parm               = h5parm.h5parm(MSObject.pathDirectory + "/" + MSObject.nameMS + ".h5", readonly = False)
        objectSolSet               = objectH5Parm.getSolset("sol000")
        objectSolTabGainAmplitudes = objectSolSet.getSoltab("amplitude000")
        objectSolTabGainPhases     = objectSolSet.getSoltab("phase000")
        gainAmplitudes             = objectSolTabGainAmplitudes.getValues(retAxesVals = False, weight = False)
        gainPhases                 = objectSolTabGainPhases.getValues(    retAxesVals = False, weight = False)
        numberOfTimeStamps         = gainAmplitudes.shape[4]
        print (gainAmplitudes.shape)
        # Fill 'gainAmplitudes' and 'gainPhases' with values from 'bandpassesAmplitude' and 'bandpassesPhase'.
        # Make 'numberOfTimeStamps' copies.
        #gainAmplitudes             =
        #gainPhases                 =
        weights                    = np.logical_not(np.isnan(gainAmplitudes))
        # Fill existing SolTabs with 'gainAmplitudes', 'gainPhases' and 'weights'.
        objectH5Parm.close()
    # 3. Save H5Parms.


    # Apply solutions to target fields.
    #MSs.run(command = "DPPP " + pathParSetApply + " msin=$pathMS " +
    #        "applyBandpassAmplitude.parmdb=$pathDirectory/$nameMS.h5 applyBandpassPhase.parmdb=$pathDirectory/$nameMS.h5",
    #        commandType = "DPPP", log = "transfer_$nameMS.log")



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