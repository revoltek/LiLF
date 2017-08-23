#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Francesco de Gasperin & Martijn Oei, 2017
In collaboration with: Tammo Jan Dijkema and Reinout van Weeren

This pipeline chunk
1. Downloads uGMRT data.
2. Converts downloaded data into MS format.
3. Sets up a series of directories.
4. Splits up the MS on a scan-by-scan basis.
5. Modifies the MSs to LOFAR format.
6. Puts the MSs in the right place in the file tree.
'''

import sys, os, shutil

from casacore import tables
import numpy as np

import lib_ms
import lib_util



def columnExists(tableObject, columnName):
    '''
    Check whether a column with name 'columnName' exists in  table 'tableObject'.
    '''
    columnNames = tableObject.colnames()
    
    return (columnName in columnNames)


def getScanIDs(pathMS, verbose = False):
    '''
    Return the different scan IDs in a MS.
    '''

    scanIDs = (tables.taql("select distinct SCAN_NUMBER from $pathMS")).getcol("SCAN_NUMBER")

    if (verbose):
        lib_util.printLineBold("scanIDs:")
        print (scanIDs)

    return scanIDs


def getTimes(pathMS, verbose = False):
    '''
    Return the different time stamp centres in a MS.
    '''

    times = (tables.taql("select distinct TIME from $pathMS")).getcol("TIME")

    if (verbose):
        lib_util.printLineBold("times:")
        print (times)

    return times


# 1. Download data from NCRA server




# 2. Convert into MS




# 3. Set-up directory structure
# Temporary:
pathDirectoryMS   = "/disks/strw3/oei/"
nameMS            = "uGMRTCosmosCut.ms"
pathMS            = pathDirectoryMS + nameMS

if (not os.path.isdir(pathMS)):
    print ("'" + pathMS + "' is not a valid path.")
    sys.exit()

pathDirectoryMain = pathDirectoryMS + nameMS[ : -3] + "-PiLF/"

if (os.path.isdir(pathDirectoryMain)):
    shutil.rmtree(pathDirectoryMain)
os.mkdir(pathDirectoryMain)
os.mkdir(pathDirectoryMain + "fieldsTarget/")
os.mkdir(pathDirectoryMain + "fieldsCalibrator/")
os.mkdir(pathDirectoryMain + "logs/")




# 4. Split up the MS.
verbose                   = False
columnNameVisibilitiesNew = "DATA"
columnNameFlagsNew        = "FLAG"
columnNameWeightsNew      = "WEIGHT_SPECTRUM"

scanIDs                   = getScanIDs(pathMS)
numberOfScans             = len(scanIDs)

lib_util.printLineBold("pathMS:")
print (pathMS)

# Create temporary paths for the sub-MSs.
pathsMSNew    = np.empty(len(scanIDs), dtype = object)
i             = 0
for scanID in scanIDs:
    pathMSNew      = pathDirectoryMain + "scanID" + str(scanID) + ".MS"
    pathsMSNew[i]  = pathMSNew
    i             += 1
lib_util.printLineBold("pathsMSNew:")
print (pathsMSNew)


lib_util.printLineBold("- Split-up of original MS '" + pathMS + "' -")
for scanID, pathMSNew, i in zip(scanIDs, pathsMSNew, np.arange(numberOfScans) + 1):
    tables.taql("SELECT from $pathMS where SCAN_NUMBER = $scanID giving $pathMSNew as plain")
    
    lib_util.printLineBold(str(i) + " / " + str(numberOfScans) + ". Created MS '" + pathMSNew + "'.")




# 5. Convert MSs from uGMRT to LOFAR format.
for pathMSNew in pathsMSNew:
    lib_util.printLineBold("Starting work on MS at '" + pathMSNew + "'...")


    # 5.1
    # Removal of columns can give errors when executing LOFAR command 'msoverview'.
    lib_util.printLineBold("- Removal of unnecessary data columns -")

    t                        = tables.table(pathMSNew, readonly = False)

    # Old array: ["EXPOSURE", "SIGMA_SPECTRUM"]. The column "EXPOSURE" seems necessary for LOFAR's msoverview and for DPPP.
    for columnName in ["SIGMA_SPECTRUM"]: # This list could possibly be expanded.
        if (columnExists(t, columnName)):
            t.removecols(columnName)

    t.close()


    # 5.2
    # CORR_TYPE    column description comment: 'The polarization type for each correlation product, as a Stokes enum.'
    # CORR_PRODUCT column description comment: 'Indices describing receptors of feed going into correlation'
    lib_util.printLineBold("- Adaptation of polarisation metadata -")

    t                        = tables.table(pathMSNew + "/POLARIZATION", readonly = False)

    correlationTypesNew      = np.array([[5, 6, 7, 8]])
    correlationProductsNew   = np.array([[[0, 0], [0, 1], [1, 0], [1, 1]]])
    numberOfCorrelationsNew  = 4

    t.putcol("CORR_TYPE",    correlationTypesNew)
    t.putcol("CORR_PRODUCT", correlationProductsNew)
    t.putcol("NUM_CORR",     numberOfCorrelationsNew)

    t.close()


    # 5.3
    lib_util.printLineBold("- Adaptation of frequency metadata -")

    t                       = tables.table(pathMSNew + "/SPECTRAL_WINDOW", readonly = False)
    frequencies             = t.getcol("CHAN_FREQ")
    frequenciesNew          = np.fliplr(frequencies)
    t.putcol("CHAN_FREQ",   frequenciesNew)
    t.close()

    if (verbose):
        lib_util.printLineBold("frequencies:")
        print (frequencies)
        lib_util.printLineBold("frequencies (updated):")
        print (frequenciesNew)


    # 5.4
    lib_util.printLineBold("- Adaptation of field information -")

    pathMSNewField = pathMSNew + "/FIELD"

    # Remove metadata of other fields in the FIELD subtable
    tables.taql("delete from $pathMSNewField where rownr() not in (select distinct FIELD_ID from $pathMSNew)")

    # Set 'SOURCE_ID' to 0 in the FIELD subtable
    tables.taql("update $pathMSNewField set SOURCE_ID=0")

    # Set 'FIELD_ID' to 0 in the main table
    tables.taql("update $pathMSNew set FIELD_ID=0")


    # 5.5
    lib_util.printLineBold("- Adaptation of intervals -")

    times           = getTimes(pathMSNew)
    intervalPrecise = times[1] - times[0]

    # Open the MS as a table in a way that changes can be made.
    t               = tables.table(pathMSNew, readonly = False)
    intervals       = t.getcol("INTERVAL")
    intervalsNew    = np.ones_like(intervals) * intervalPrecise
    t.putcol("INTERVAL", intervalsNew)
    t.close()

    if (verbose):
        lib_util.printLineBold("Time intervals (should be equal):")
        print (times[1:] - times[:-1])


    # 5.6
    lib_util.printLineBold("- Change existing or create alternative columns for data, flags and weights -")

    # Open the MS as a table in a way that changes can be made.
    t                        = tables.table(pathMSNew, readonly = False)

    # First, change the 'direction' of the frequency axis.
    # If the visibilities were sorted in descending frequency order (id est given for high frequencies first),
    # they are converted to ascending frequency order. Note: applying this operation twice returns the MS to its old state!
    if (verbose):
        lib_util.printLineBold("Loading visibilities...")
    visibilities             = t.getcol("DATA")
    if (verbose):
        lib_util.printLineBold("Visibilities loaded!")

    visibilities             = np.fliplr(visibilities)
    visibilitiesNew          = np.zeros((visibilities.shape[0], visibilities.shape[1], 4), dtype = np.complex128)
    visibilitiesNew[:, :, 0] = visibilities[:, :, 0]
    visibilitiesNew[:, :, 3] = visibilities[:, :, 1]

    keywordNames             = t.colkeywordnames("DATA")
    columnDescription        = t.getcoldesc("DATA")
    dataManagerInfo          = t.getdminfo("DATA")

    if (verbose):
        lib_util.printLineBold("keywordNames:")
        print (keywordNames)
        lib_util.printLineBold("columnDescription:")
        print (columnDescription)
        lib_util.printLineBold("dataManagerInfo:")
        print (dataManagerInfo)

    dataManagerInfo["NAME"]                                  = "TiledDATAMartijn"
    dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, visibilities.shape[1], visibilities.shape[0]], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, visibilities.shape[1]], dtype = np.int32)

    if (verbose):
        lib_util.printLineBold("dataManagerInfo (updated):")
        print (dataManagerInfo)

    if (verbose):
        lib_util.printLineBold("Removing column '" + columnNameVisibilitiesNew + "', if it exists...")
    if (columnExists(t, columnNameVisibilitiesNew)):
        t.removecols(columnNameVisibilitiesNew)

    if (verbose):
        lib_util.printLineBold("Adding column '" + columnNameVisibilitiesNew + "'...")
    t.addcols(tables.makecoldesc(columnNameVisibilitiesNew, columnDescription), dataManagerInfo)

    if (verbose):
        lib_util.printLineBold("Filling column '" + columnNameVisibilitiesNew + "'...")
    t.putcol(columnNameVisibilitiesNew, visibilitiesNew)

    if (verbose):
        lib_util.printLineBold("Visibilities flipped along frequency axis and placeholder polarisations added!")



    if (verbose):
        lib_util.printLineBold("Loading flags...")
    flags                    = t.getcol("FLAG")
    if (verbose):
        lib_util.printLineBold("Flags loaded!")

    flags                    = np.fliplr(flags)
    flagsNew                 = np.zeros((flags.shape[0], flags.shape[1], 4), dtype = np.bool_)
    flagsNew[:, :, 0]        = flags[:, :, 0]
    flagsNew[:, :, 1]        = flags[:, :, 0] # Take over flags from LL correlation
    flagsNew[:, :, 2]        = flags[:, :, 0] # Take over flags from LL correlation
    flagsNew[:, :, 3]        = flags[:, :, 1]

    keywordNames             = t.colkeywordnames("FLAG")
    columnDescription        = t.getcoldesc("FLAG")
    dataManagerInfo          = t.getdminfo("FLAG")

    if (verbose):
        lib_util.printLineBold("keywordNames:")
        print (keywordNames)
        lib_util.printLineBold("columnDescription:")
        print (columnDescription)
        lib_util.printLineBold("dataManagerInfo:")
        print (dataManagerInfo)

    dataManagerInfo["NAME"]                                  = "TiledFlagMartijn"
    dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, flags.shape[1], flags.shape[0]], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, flags.shape[1]], dtype = np.int32)

    if (verbose):
        lib_util.printLineBold("dataManagerInfo (updated):")
        print (dataManagerInfo)

    if (verbose):
        lib_util.printLineBold("Removing column '" + columnNameFlagsNew + "', if it exists...")
    if (columnExists(t, columnNameFlagsNew)):
        t.removecols(columnNameFlagsNew)

    if (verbose):
        lib_util.printLineBold("Adding column '" + columnNameFlagsNew + "'...")
    t.addcols(tables.makecoldesc(columnNameFlagsNew, columnDescription), dataManagerInfo)

    if (verbose):
        lib_util.printLineBold("Filling column '" + columnNameFlagsNew + "'...")
    t.putcol(columnNameFlagsNew, flagsNew)

    if (verbose):
        lib_util.printLineBold("Flags flipped along frequency axis and placeholder polarisations added!")



    if (verbose):
        lib_util.printLineBold("Loading weights...")
    weights                  = t.getcol("WEIGHT_SPECTRUM")
    if (verbose):
        lib_util.printLineBold("Weights loaded!")

    weights                  = np.fliplr(weights)
    weightsNew               = np.zeros((weights.shape[0], weights.shape[1], 4), dtype = np.float64)
    weightsNew[:, :, 0]      = weights[:, :, 0]
    weightsNew[:, :, 3]      = weights[:, :, 1]

    keywordNames             = t.colkeywordnames("WEIGHT_SPECTRUM")
    columnDescription        = t.getcoldesc("WEIGHT_SPECTRUM")
    dataManagerInfo          = t.getdminfo("WEIGHT_SPECTRUM")
    if (verbose):
        lib_util.printLineBold("keywordNames:")
        print (keywordNames)
        lib_util.printLineBold("columnDescription:")
        print (columnDescription)
        lib_util.printLineBold("dataManagerInfo:")
        print (dataManagerInfo)

    dataManagerInfo["NAME"]                                  = "TiledWgtSpectrumMartijn"
    dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, weights.shape[1], weights.shape[0]], dtype = np.int32)
    dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, weights.shape[1]], dtype = np.int32)

    if (verbose):
        lib_util.printLineBold("dataManagerInfo (updated):")
        print (dataManagerInfo)

    if (verbose):
        lib_util.printLineBold("Removing column '" + columnNameWeightsNew + "', if it exists...")
    if (columnExists(t, columnNameWeightsNew)):
        t.removecols(columnNameWeightsNew)

    if (verbose):
        lib_util.printLineBold("Adding column '" + columnNameWeightsNew + "'...")
    t.addcols(tables.makecoldesc(columnNameWeightsNew, columnDescription), dataManagerInfo)

    if (verbose):
        lib_util.printLineBold("Filling column '" + columnNameWeightsNew + "'...")
    t.putcol(columnNameWeightsNew, weightsNew)

    if (verbose):
        lib_util.printLineBold("Weights flipped along frequency axis and placeholder polarisations added!")


    t.close()

    lib_util.printLineBold("Finished work on MS at '" + pathMSNew + "'!")




# 6. Move the sub-MSs from a temporary place to their right locus in the file tree.
for pathMSNew in pathsMSNew:
    
    MSObject    = lib_ms.Ms(pathMSNew)
    
    if (MSObject.isCalibrator()):
        pathDirectoryCalibrator = pathDirectoryMain + "fieldsCalibrator/" + pathMSNew[len(pathDirectoryMain) : -3]
        pathMSFinal             = pathDirectoryCalibrator + "/" + pathMSNew[len(pathDirectoryMain) : ]
        if (not os.path.isdir(pathDirectoryCalibrator)):
            os.mkdir(pathDirectoryCalibrator)
            os.mkdir(pathDirectoryCalibrator + "/plots")
            os.mkdir(pathDirectoryCalibrator + "/solutions")
    else:
        pathDirectoryTarget = pathDirectoryMain + "fieldsTarget/" + MSObject.getName()
        pathMSFinal         = pathDirectoryTarget + "/MSs/" + pathMSNew[len(pathDirectoryMain) : ]
        if (not os.path.isdir(pathDirectoryTarget)):
            os.mkdir(pathDirectoryTarget)
            os.mkdir(pathDirectoryTarget + "/plots")
            os.mkdir(pathDirectoryTarget + "/MSs")
            os.mkdir(pathDirectoryTarget + "/images")

    lib_util.printLineBold("Moving MS at '" + pathMSNew + "' to '" + pathMSFinal + "'...")
    shutil.move(pathMSNew, pathMSFinal)
    MSObject.ms = pathMSFinal
    
    print (MSObject.getName())
    print (MSObject.isCalibrator())