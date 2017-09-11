#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse, logging

from losoto import h5parm
from matplotlib import cm
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy

import lib_util


def wrapPhasesZeroCentred(phases, unitDegree = True):
        '''
        This method assumes that 'phases' is expressed on a -180 to 180 degree scale (or -pi to pi if 'unitDegree' is False),
        with some of the phases lying out of bounds - which is the motivation to call this function.
        It correctly phase wraps these phases, and returns them expressed on the same scale.
        '''
        if (unitDegree):
            # Move the phases to a 0 to 360 degrees scale (with some out of bounds).
            phases += 180

            # Next, put all phases between 0 and 360 degrees.
            # Note: 'numpy.mod(...)' correctly handles negative numbers. So e.g. 'numpy.mod([-10, -365], 360) == [350, 355]'.
            phases  = numpy.mod(phases, 360)

            # Move back to a -180 to 180 degree scale (with none out of bounds).
            phases -= 180
        else:
            phases += numpy.pi
            phases  = numpy.mod(phases, 2 * numpy.pi)
            phases -= numpy.pi

        return phases


def plotAmplitudes2D(amplitudes, antennaeWorking, pathDirectoryPlots, namePolarisation, nameField, timeStart, timeRange,
                     frequencyStart = 300, frequencyRange = 200, nameH5Parm = "?", nameTelescope = "uGMRT"):
    """
    Generate time-frequency plots of antenna-based gain amplitudes, of one specific polarisation.
    """

    # Establish data properties.
    numberOfAntennae, numberOfChannels, numberOfTimeStamps = amplitudes.shape

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            logging.info("Starting gain amplitudes visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + "...")

            # Create 2D antenna-based gain amplitudes plot.
            figure       = pyplot.figure(figsize = (12, 6))
            image        = pyplot.imshow(amplitudes[i], aspect = "auto", interpolation = "nearest", cmap = cm.viridis, vmin = 0, vmax = 1)

            pyplot.xlabel("time (s)")
            pyplot.ylabel("frequency (MHz)")
            pyplot.xticks(numpy.linspace(0, numberOfTimeStamps - 1, num = 5, endpoint = True),
                          numpy.linspace(timeStart, timeStart + timeRange, num = 5, endpoint = True).astype(int))
            pyplot.yticks(numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
                          numpy.linspace(frequencyStart, frequencyStart + frequencyRange, num = 5, endpoint = True).astype(int))
            pyplot.title("antenna-based gain amplitudes of uncalibrated calibrator visibilities\ndata set: "
                         + nameH5Parm + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
                         + nameField + " | polarisation: " + namePolarisation)

            colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
            colorBar     = pyplot.colorbar(image, cax = colorBarAxis, ticks = [0, 0.2, 0.4, 0.6, 0.8, 1])
            colorBar.ax.set_ylabel("gain amplitude $(1)$")

            pyplot.subplots_adjust(left = .06, right = .94, bottom = 0.08, top = 0.91)
            pyplot.savefig(pathDirectoryPlots + "/amplitudes2D_ant" + str(i) + "_pol" + namePolarisation + ".pdf")
            pyplot.close()
        else:
            logging.info("Skipping gain amplitudes visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + ": all data are flagged.")


def plotPhases2D(phases, antennaeWorking, pathDirectoryPlots, namePolarisation, nameField, timeStart, timeRange,
                 frequencyStart = 300, frequencyRange = 200, nameH5Parm = "?", nameTelescope = "uGMRT"):
    """
    Generate time-frequency plots of antenna-based gain phases, of one specific polarisation.
    """

    # Establish data properties.
    numberOfAntennae, numberOfChannels, numberOfTimeStamps = phases.shape

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            logging.info("Starting gain phases visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + "...")

            # Create 2D antenna-based gain phases plot.
            figure       = pyplot.figure(figsize = (12, 6))
            image        = pyplot.imshow(phases[i], aspect = "auto", interpolation = "none", cmap = cm.hsv, vmin = -180, vmax = 180)

            pyplot.xlabel("time (s)")
            pyplot.ylabel("frequency (MHz)")
            pyplot.xticks(numpy.linspace(0, numberOfTimeStamps - 1, num = 5, endpoint = True),
                          numpy.linspace(timeStart, timeStart + timeRange, num = 5, endpoint = True).astype(int))
            pyplot.yticks(numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
                          numpy.linspace(frequencyStart, frequencyStart + frequencyRange, num = 5, endpoint = True).astype(int))
            pyplot.title("antenna-based gain phases of uncalibrated calibrator visibilities\ndata set: "
                         + nameH5Parm + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
                         + nameField + " | polarisation: " + namePolarisation)

            colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
            colorBar     = pyplot.colorbar(image, cax = colorBarAxis, ticks = [-180, -120, -60, 0, 60, 120, 180])
            colorBar.ax.set_ylabel("gain phase $(\degree)$")

            pyplot.subplots_adjust(left = .06, right = .94, bottom = 0.08, top = 0.91)
            pyplot.savefig(pathDirectoryPlots + "/phases2D_ant" + str(i) + "_pol" + namePolarisation + ".pdf")
            pyplot.close()
        else:
            logging.info("Skipping gain phases visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + ": all data are flagged.")


def dedicated_uGMRT_bandpass(pathH5Parm, referenceAntennaID = 0, verbose = False):

    # Load H5Parm file.
    objectH5Parm             = h5parm.h5parm(pathH5Parm)
    objectH5Parm.printInfo()

    # Load antenna-based gains.
    gainAmplitudes           = (objectH5Parm.H.root.sol000.amplitude000.val)   [ : , 0, : , : ]
    gainPhases               = (objectH5Parm.H.root.sol000.phase000.val)       [ : , 0, : , : ]

    # Load weights (generalised flags).
    weightsForAmplitudes     = (objectH5Parm.H.root.sol000.amplitude000.weight)[ : , 0, : , : ]
    weightsForPhases         = (objectH5Parm.H.root.sol000.phase000.weight)    [ : , 0, : , : ]

    # Split-up gains by polarisation.
    gainAmplitudesPol1       = gainAmplitudes[0]
    gainAmplitudesPol2       = gainAmplitudes[1]
    gainPhasesPol1           = gainPhases[0]
    gainPhasesPol2           = gainPhases[1]

    # Split-up weights by polarisation.
    weightsForAmplitudesPol1 = weightsForAmplitudes[0]
    weightsForAmplitudesPol2 = weightsForAmplitudes[1]
    weightsForPhasesPol1     = weightsForPhases[0]
    weightsForPhasesPol2     = weightsForPhases[1]

    # Establish data properties.
    numberOfAntennae, numberOfChannels, numberOfTimeStamps = gainAmplitudesPol1.shape

    # Make gain phases relative to reference antenna.
    gainPhasesPol1 -= numpy.tile(gainPhasesPol1[referenceAntennaID, : , : ], (numberOfAntennae, 1, 1))
    gainPhasesPol2 -= numpy.tile(gainPhasesPol2[referenceAntennaID, : , : ], (numberOfAntennae, 1, 1))
    gainPhasesPol1  = wrapPhasesZeroCentred(gainPhasesPol1, unitDegree = False)
    gainPhasesPol2  = wrapPhasesZeroCentred(gainPhasesPol2, unitDegree = False)

    # Convert gain phases from degrees to radians.
    gainPhasesPol1  = numpy.degrees(gainPhasesPol1)
    gainPhasesPol2  = numpy.degrees(gainPhasesPol2)
    print(numpy.amax(gainPhasesPol1), numpy.amin(gainPhasesPol2))

    print ((objectH5Parm.H.root.sol000.amplitude000.val).shape)

    #for valsThisTime, weights, coord, selection in getValuesIter(returnAxes = ["time",'freq'], weights = True):
    #    valsThisTime *= 2
    #    setValues(selection = selection)

    # These values can be taken from the MS, and perhaps also from the H5Parm file.
    # Temporary!
    pathDirectoryPlots = "/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsCalibrator/scanID1/plots"
    nameField          = "3C147"
    timeStart          = 0                                    # in seconds
    timeStampLength    = 8.05                                 # in seconds
    timeRange          = numberOfTimeStamps * timeStampLength # in seconds

    # Plot gain amplitudes.
    plotAmplitudes2D(gainAmplitudesPol1, [True] * numberOfAntennae, pathDirectoryPlots, "LL", nameField, timeStart, timeRange)
    plotAmplitudes2D(gainAmplitudesPol2, [True] * numberOfAntennae, pathDirectoryPlots, "RR", nameField, timeStart, timeRange)

    # Plot gain phases.
    plotPhases2D(    gainPhasesPol1,     [True] * numberOfAntennae, pathDirectoryPlots, "LL", nameField, timeStart, timeRange)
    plotPhases2D(    gainPhasesPol2,     [True] * numberOfAntennae, pathDirectoryPlots, "RR", nameField, timeStart, timeRange)


    # Output debug info.
    print (gainAmplitudes.shape)
    print (gainPhases.shape)
    print (gainAmplitudesPol1.shape)
    print (gainPhasesPol1.shape)
    print (weightsForAmplitudesPol2.shape)
    print (weightsForPhasesPol2.shape)
    print ("numberOfAntennae:", numberOfAntennae, "numberOfChannels:", numberOfChannels, "numberOfTimeStamps:", numberOfTimeStamps)

    import sys
    sys.exit()

    '''
    In this step, the amplitude bandpasses for each antenna are determined in two iterations.
    Iteration 1 is subdivided into several 'subiterations' that are meant to flag outliers.
    Iteration 2 takes the bandpass results from iteration 1 and performs median filtering and interpolation.
    '''

    # Create lists that store the normalised amplitude bandpasses and the bandpass normalisation factors.
    bandpassesAmplitudePol1          = []
    bandpassesAmplitudePol2          = []
    bandpassNormalisationFactorsPol1 = []
    bandpassNormalisationFactorsPol2 = []

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            print ("Starting amplitude bandpass calculation for antenna ID " + str(i) + "...")

            for subIteration in range(numberOfSubIterationsBandpassAmplitude):
                # Determine the amplitude bandpass (iteration 1).
                if (subIteration == numberOfSubIterationsBandpassAmplitude - 1): # If in the last iteration, determine the amplitude bandpass using the mean.
                    bandpassAmplitudePol1Iter1 = numpy.nanmean(gainAmplitudesPol1[i], axis = 1)
                    bandpassAmplitudePol2Iter1 = numpy.nanmean(gainAmplitudesPol2[i], axis = 1)
                else: # During earlier iterations, to avoid e.g. RFI signatures, we use the median.
                    bandpassAmplitudePol1Iter1 = numpy.nanmedian(gainAmplitudesPol1[i], axis = 1)
                    bandpassAmplitudePol2Iter1 = numpy.nanmedian(gainAmplitudesPol2[i], axis = 1)

                # Tile the bandpass into a grid.
                gridBandpassAmplitudePol1 = numpy.tile(bandpassAmplitudePol1Iter1, (numberOfTimeStamps, 1)).T
                gridBandpassAmplitudePol2 = numpy.tile(bandpassAmplitudePol2Iter1, (numberOfTimeStamps, 1)).T

                # Remove the bandpass from the data. Residuals will be centered around 1.
                gridResidualsPol1 = numpy.divide(gainAmplitudesPol1[i], gridBandpassAmplitudePol1)
                gridResidualsPol2 = numpy.divide(gainAmplitudesPol2[i], gridBandpassAmplitudePol2)

                # Calculate the standard deviation of the residuals.
                STDPol1 = numpy.nanstd(gridResidualsPol1)
                STDPol2 = numpy.nanstd(gridResidualsPol2)

                # Determine outliers and update the flags. This makes sure that the bandpass in the next iterations is better.
                gridIsOutlierPol1 = numpy.greater(numpy.absolute(gridResidualsPol1 - 1), flaggingThresholdFactorAmplitude * STDPol1)
                gridIsOutlierPol2 = numpy.greater(numpy.absolute(gridResidualsPol2 - 1), flaggingThresholdFactorAmplitude * STDPol2)
                flagsPol1[i] = numpy.logical_or(flagsPol1[i], gridIsOutlierPol1)
                flagsPol2[i] = numpy.logical_or(flagsPol2[i], gridIsOutlierPol2)

                # Set flagged data to 'numpy.nan'.
                gainAmplitudesPol1[i] = numpy.where(flagsPol1[i], numpy.nan, gainAmplitudesPol1[i])
                gainAmplitudesPol2[i] = numpy.where(flagsPol2[i], numpy.nan, gainAmplitudesPol2[i])
                gainPhasesPol1[i]     = numpy.where(flagsPol1[i], numpy.nan, gainPhasesPol1[i])
                gainPhasesPol2[i]     = numpy.where(flagsPol2[i], numpy.nan, gainPhasesPol2[i])
#
#                 # Median filter and interpolate the amplitude bandpass (iteration 2).
#                 bandpassAmplitudePol1Iter2 = fillGaps1D(filterMedian1D(bandpassAmplitudePol1Iter1, kernelSize = 7))
#                 bandpassAmplitudePol2Iter2 = fillGaps1D(filterMedian1D(bandpassAmplitudePol2Iter1, kernelSize = 7))
#
#                 '''
#                 # Use this piece of code to debug the second iteration of the amplitude bandpass.
#                 pyplot.plot(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter1, c = 'r', linestyle = "--")
#                 pyplot.plot(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter2, c = 'b', linestyle = ":")
#                 pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter1, s = 4, c = 'r')
#                 pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter2, s = 4, c = 'b')
#                 pyplot.show()
#                 sys.exit()
#                 '''
#
#                 # Normalise the amplitude bandpass for both polarisations.
#                 # We do not normalise by determining the factor that scales the amplitude bandpass peak to 1,
#                 # as this would result in wrong scaling behaviour in cases where the 'real' peak is absent in the
#                 # bandpass due to flagging.
#                 # We also do not normalise by determining the factor that would scale the mean of the bandpass
#                 # to 1, as this requires taking the mean of all channels except for those flagged. In such case,
#                 # when many channels at the lower and upper end of the frequency band are non-flagged, the scaling
#                 # is different than when many of those edge channels are flagged.
#                 # We conclude that the best thing to do is to normalise by determining the factor that would scale
#                 # the mean bandpass in a selected range of frequency channels to 1. These channels are chosen such
#                 # that the bandpass is approximately constant over them.
#                 # If we normalise bandpass iteration 2, then we could use 'numpy.mean' instead of 'numpy.nanmean' as well
#                 # in the determination of 'bandpassNormalisationFactorPolX'. After all, due to the linear interpolation
#                 # only the edges of the frequency band still lack a bandpass value. We still prefer using 'numpy.nanmean'
#                 # for antennae with a lot of flagged channels at the outer parts of the frequency band, as in such cases
#                 # a few 'numpy.nan's might lie into the domain '[frequencyNormalisationStart, frequencyNormalisationEnd]'.
#
#
#                 # Determine the normalisation factors.
#                 indexStart = numpy.argmax(gridFrequencies[ : , 0] > frequencyNormalisationStart)
#                 indexEnd   = numpy.argmax(gridFrequencies[ : , 0] > frequencyNormalisationEnd)
#                 bandpassNormalisationFactorPol1 = numpy.nanmean(bandpassAmplitudePol1Iter2[indexStart : indexEnd]) # 'numpy.mean' would work too, mostly
#                 bandpassNormalisationFactorPol2 = numpy.nanmean(bandpassAmplitudePol2Iter2[indexStart : indexEnd]) # 'numpy.mean' would work too, mostly
#
#                 # Divide the amplitude bandpasses by the normalisation factor.
#                 bandpassAmplitudePol1Iter1 = numpy.divide(bandpassAmplitudePol1Iter1, bandpassNormalisationFactorPol1)
#                 bandpassAmplitudePol2Iter1 = numpy.divide(bandpassAmplitudePol2Iter1, bandpassNormalisationFactorPol2)
#                 bandpassAmplitudePol1Iter2 = numpy.divide(bandpassAmplitudePol1Iter2, bandpassNormalisationFactorPol1)
#                 bandpassAmplitudePol2Iter2 = numpy.divide(bandpassAmplitudePol2Iter2, bandpassNormalisationFactorPol2)
#
#                 # Add the normalised amplitude bandpasses and normalisation factors to the appropriate lists.
#                 bandpassesAmplitudePol1.append(bandpassAmplitudePol1Iter2)
#                 bandpassesAmplitudePol2.append(bandpassAmplitudePol2Iter2)
#                 bandpassNormalisationFactorsPol1.append(bandpassNormalisationFactorPol1)
#                 bandpassNormalisationFactorsPol2.append(bandpassNormalisationFactorPol2)
#
#
#                 if (plot):
#                     print ("Starting amplitude bandpass visualisation for antenna ID " + str(i) + "...")
#
#                     # Create plot of amplitude bandpass (for both polarisations, iteration 1).
#                     pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter1, c = "navy", s = 16, lw = 0, label = "polarisation 1\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol1, 3)))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol2Iter1, c = "orangered", s = 16, lw = 0, label = "polarisation 2\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol2, 3)))
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.xlabel("frequency channel centre (MHz)")
#                     pyplot.ylabel("antenna-based gain amplitude (1)")
#                     pyplot.xlim(frequencyStart - plotFrequencyLimit, frequencyStart + frequencyRange + plotFrequencyLimit)
#                     pyplot.ylim(0, 2 + plotAmplitudeLimit)
#                     pyplot.title("normalised amplitude bandpass (iteration 1)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .07, right = .98, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "bandpassAmplitudeIter1_" + "Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     # Create plot of amplitude bandpass (for both polarisations, iteration 2).
#                     pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter2, c = "navy", s = 16, lw = 0, label = "polarisation 1\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol1, 3)))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol2Iter2, c = "orangered", s = 16, lw = 0, label = "polarisation 2\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol2, 3)))
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.xlabel("frequency channel centre (MHz)")
#                     pyplot.ylabel("antenna-based gain amplitude (1)")
#                     pyplot.xlim(frequencyStart - plotFrequencyLimit, frequencyStart + frequencyRange + plotFrequencyLimit)
#                     pyplot.ylim(0, 2 + plotAmplitudeLimit)
#                     pyplot.title("normalised amplitude bandpass (iteration 2)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .07, right = .98, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "bandpassAmplitudeIter2_" + "Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     if (plot3D):
#                         dataGrids = [gridBandpassAmplitudePol1, gridBandpassAmplitudePol2]
#                         for dataGrid, polarisationName in zip(dataGrids, polarisationNames):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             axes3D.plot_surface(gridFrequencies, gridTimes, numpy.where(numpy.isnan(dataGrid), 0, dataGrid), cmap = cm.Spectral_r, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("amplitude ($1$)")
#                             axes3D.set_zlim(0, 1)
#                             axes3D.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1])
#                             axes3D.set_title("unnormalised amplitude bandpass (iteration 1)\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "bandpassAmplitude3DIter1_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#
#                         # Create 3D plots of antenna-based gain amplitude residuals.
#                         dataGrids = [gridResidualsPol1, gridResidualsPol2]
#                         STDs = [STDPol1, STDPol2]
#                         for dataGrid, STD, polarisationName in zip(dataGrids, STDs, polarisationNames):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             axes3D.plot_surface(gridFrequencies, gridTimes, numpy.where(numpy.isnan(dataGrid), 1, dataGrid), cmap = cm.YlOrRd, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("amplitude (1)")
#                             axes3D.set_zlim(1 - .1, 1 + .1)
#                             axes3D.set_title("antenna-based gain amplitude residuals (iteration 1) | $\sigma$ = " + str(numpy.round(STD, 3)) + "\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "residualsAmplitudeIter1_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#             else:
#                 print ("Skipping amplitude bandpass calculation for antenna ID " + str(i) + ": all data are flagged.")
#
#         # Save the normalised amplitude bandpasses and their normalisation factors, after all working antennae have been looped over.
#         numpy.save(pathOutputBPsAmplitude, [bandpassesAmplitudePol1, bandpassesAmplitudePol2, bandpassNormalisationFactorsPol1, bandpassNormalisationFactorsPol2, antennaeWorking])
#

    # After amplitude and phase bandpass and TECs are created, save back to H5Parm file.
    # Write the TEC solutions to 'objectH5Parm.H.root.sol000.tec000.val'.


if (__name__ == "__main__"):
    # If the program is run from the command line, parse arguments.
    parser                      = argparse.ArgumentParser(description = "Pipeline step 3: Generation of bandpasses.")
    parser.add_argument("pathH5Parm", help = "Path of the H5Parm file to act upon.")
    arguments                   = parser.parse_args()

    lib_util.printLineBold("Parameters to use:")
    print (arguments)

    dedicated_uGMRT_bandpass(arguments.pathH5Parm)

# '''
# Martijn Oei, 2017
# This program:
#
# 1. Takes the measured and modelled visibilities of a single calibrator field from a Measurement Set,
# and calculates antenna-based gains for many time-frequency bins for both polarisations seperately.
# The resulting Calibration Table is stored and read-out again. A NumPy array containing the same data,
# but reordered in a more convenient format, is generated.
#
# 2. Visualises the gain amplitudes and phases by creating both 2D image and 3D landscape plots.
#
# 3. Finds the amplitude and phase bandpasses in multi-step processes, during which some flagging is performed.
# The various steps are visualised. The phase bandpass calculation also generates a function DeltaTEC(t) for each antenna.
# '''
# #!/usr/local/bin/python3.5 # This line might not be required on your system.
#
# import argparse, sys
#
# from matplotlib import cm
# from matplotlib import colors
# from matplotlib import pyplot
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import losoto, numpy
#
# #import libraryRadioAstronomy
#
#
# def dedicated_uGMRT_bandpass(pathH5Parm, verbose = False):
#
#     # Load the H5Parm file.
#     fileH5Parm = losoto.h5parm.h5parm(pathH5Parm)
#
#     if (verbose):
#         fileH5Parm.printInfo()
#
#     solutionSet = fileH5Parm.getSolset("sol000")
#     amplitudes  = solutionSet.amp000.val[ : , 0, : , : ]
#     phases      = solutionSet.phase000.val[ : , 0, : , : ]
#
#
#     # import H5 to numpy!!!!!!
#     # ParmDB to H5
#     # Losoto H5parm_importer.py
#     # pip install --allow-external --upgrade --user https://github.com/revoltek/losoto/archive/master.zip # Dit installt het in Python
#     # Dan, ergens in een directory dumpt-ie:
#     # git clone https://github.com/revoltek/losoto.git # Dit zorgt ervoor dat je het in de command line kunt gebruiken
#     # ../rvweeren/software/losoto/bin/losoto
#     #from losoto import h5parm
#     #H5 = losoto.h5parm.h5parm("DI.circ.H5")
#     #handy:
#     #H5.printInfo()
#     #sols = H5.getSolset("sol000")
#     #numpyArray = sols.phase000.val[:, 0, :, :] # would load for all polarisations, for the first (and hopefully only)
#     #H5parm_importer.py -v LXXXXX_3c295.h5 LXXXXX_3c295.MS
#     #for MSObject in MSs.get_list_obj():
#     #    MSs.
#     #losoto.H5parm_importer(
#     # copy amp from cd h5parm to iono h5parm
#     #s.add('H5parm_merge.py cal-cd.h5:sol000 cal-iono.h5:solcd', log='losoto-iono.log', log_append=True, cmd_type = 'python', processors='max')
#     #s.run(check = True)
#     #from losoto.h5parm import h5parm
#     #H = h5parm(h5parmFile, readonly = True)
#     #soltabs = H.getSoltabs(solset = solsetTo)
#     #numpyArray = sols.val(...)
#
#     # After amplitude and phase bandpass and TECs are created, save back to H5.
#     # Also write back the TEC solutions in sols.tec000.val
#
#     def filterMedian1D(array, kernelSize):
#         '''
#         Median filter a 1D array.
#         'kernelSize' must be an odd integer, and at least 3.
#         '''
#         kernelRadius = int((kernelSize - 1) / 2)
#
#         arrayNew = numpy.zeros_like(array)
#         arrayNew[ : kernelRadius] = array[ : kernelRadius]
#         arrayNew[-kernelRadius : ] = array[-kernelRadius : ]
#         for index in range(kernelRadius, len(array) - kernelRadius):
#             arrayNew[index] = numpy.nanmedian(array[index - kernelRadius : index + kernelRadius + 1])
#
#         return arrayNew
#
#
#     def fillGaps1D(array):
#         '''
#         Replaces 'numpy.nan' values wherever they lie between parts of 'array' that have measured values
#         via linear interpolation.
#         Occurences of 'numpy.nan' at the beginning or end of the array are not replaced, as it is unclear
#         how this should be done via interpolation.
#         '''
#         indexDatumLast = None # index of the last encountered valid value (not 'numpy.nan')
#         arrayNew = numpy.copy(array)
#
#         for index in range(len(array)):
#             if (not numpy.isnan(array[index])):
#
#                 if (not (indexDatumLast == None) and indexDatumLast < index - 1):
#                     interpolation = numpy.linspace(array[indexDatumLast], array[index], num = index - indexDatumLast + 1, endpoint = True)[1 : -1]
#                     arrayNew[indexDatumLast + 1 : index] = interpolation
#                 indexDatumLast = index
#
#         return arrayNew
#
#
#     def wrapPhasesCenter0(phaseData):
#         '''
#         This method assumes that 'phaseData' is expressed on a -180 to 180 degree scale, with
#         some of the phases lying out of bounds (which is the motivation to call this function).
#         It correctly phase wraps these phases, and returns them expressed on the same scale.
#         '''
#         # Move the phases to a 0 to 360 degrees scale (with some out of bounds).
#         phaseData += 180
#         # Next, put all phases between 0 and 360 degrees.
#         # Note: 'numpy.mod(...)' correctly handles negative numbers. So e.g. 'numpy.mod([-10, -365], 360) == [350, 355]'.
#         phaseData = numpy.mod(phaseData, 360)
#         # Move back to a -180 to 180 degree scale (with none out of bounds).
#         phaseData -= 180
#
#         return phaseData
#
#
#     def computeDerivative1D(array, stepSize, degree = 1, intermediateSampling = False):
#         '''
#         This function computes a numerical approximation to the first or second derivative, determined via 'degree'.
#
#         'intermediateSampling = False': if coordinates are t_0, t_1, ..., t_{N-2}, t_{N-1},
#         then the derivatives are given at matching points t_1, t_2, ..., t_{N-3}, t_{N-2}.
#
#         'intermediateSampling = True': if coordinats are t_0, t_1, ..., t_{N_2}, t_{N-1},
#         then the derivatives are given at intermediate points (t_0 + t_1) / 2, ..., (t_{N-2} + t_{N-1}) / 2.
#         '''
#         if (degree == 0):
#             arrayResult = array
#         elif (degree == 1):
#             if (intermediateSampling):
#                 arrayResult = (array[1 : ] - array[ : -1]) / stepSize
#             else:
#                 arrayResult = (array[2 : ] - array[ : -2]) / (2 * stepSize)
#         elif (degree == 2):
#             if (intermediateSampling):
#                 print ("Error: intermediate sampling is not yet featured for the second derivative.")
#                 sys.exit()
#             else:
#                 arrayResult = (array[4 : ] + array[ : -4] - 2 * array[2 : -2]) / (4 * stepSize ** 2)
#         else:
#             print ("Error: 'computeDerivative1D' does not support calculating derivatives of degree " + str(degree) + ".")
#             sys.exit()
#
#         return arrayResult
#
#
#     def computeDerivative2D(grid, stepSize, axis = 0, degree = 1, intermediateSampling = False):
#         '''
#         This function computes a numerical approximation to the first or second derivative.
#         Derivatives can either be taken along the rows ('axis = 0'), or along the columns ('axis = 1').
#
#         If the coordinates of the columns of 'grid' are t_0, t_1, t_2, ..., t_{N-1}, then
#         the derivatives are given at t_1, t_2, ..., t_{N-2}.
#
#         Sometimes, e.g. for numerical integration, it is desirable that derivatives are given at
#         the intermediate points: (t_0 + t_1) / 2, (t_1 + t_2) / 2, ..., (t_{N-2} + t_{N-1}) / 2.
#         In such cases, set 'intermediateSampling' to True. Beware that output array dimensions change!
#         '''
#         if (axis == 0):
#             if (degree == 0):
#                 gridResult = grid
#             elif (degree == 1):
#                 if (intermediateSampling):
#                     gridResult = (grid[1 : , : ] - grid[ : -1, : ]) / stepSize
#                 else:
#                     gridResult = (grid[2 : , : ] - grid[ : -2, : ]) / (2 * stepSize)
#             elif (degree == 2):
#                 if (intermediateSampling):
#                     print ("Error: intermediate sampling is not yet featured for the second derivative.")
#                     sys.exit()
#                 else:
#                     gridResult = (grid[4 : , : ] + grid[ : -4, : ] - 2 * grid[2 : -2, : ]) / (4 * stepSize ** 2)
#             else:
#                 print ("Error: 'computeDerivative2D' does not support calculating derivatives of degree " + str(degree) + ".")
#                 sys.exit()
#         else:
#             if (degree == 0):
#                 gridResult = grid
#             elif (degree == 1):
#                 if (intermediateSampling):
#                     gridResult = (grid[ : , 1 : ] - grid[ : , : -1]) / stepSize
#                 else:
#                     gridResult = (grid[ : , 2 : ] - grid[ : , : -2]) / (2 * stepSize)
#             elif (degree == 2):
#                 if (intermediateSampling):
#                     print ("Error: intermediate sampling is not yet featured for the second derivative.")
#                     sys.exit()
#                 else:
#                     gridResult = (grid[ : , 4 : ] + grid[ : , : -4] - 2 * grid[ : , 2 : -2]) / (4 * stepSize ** 2)
#             else:
#                 print ("Error: 'computeDerivative2D' does not support calculating derivatives of degree " + str(degree) + ".")
#                 sys.exit()
#
#         return gridResult
#
#
#
#     # Initialise settings specific to the MS and telescope used.
#     nameMS                      = "DDTB247-GWB_FULLPOL.MS" # e.g. "DDTB247-GWB_FULLPOL.MS" or "31_008_23feb2017_gwb4.MS"
#     pathMS                      = "/data1/MartijnOei/400MUGSPilotPilot/" + nameMS # e.g. "/data1/MartijnOei/400MUGSPilotPilot/" + ... or "/data1/MartijnOei/400MUGSPilot/" + ...
#     nameField                   = "3C286" # calibrator to use - e.g. "3C286" or "3C147"
#     nameScan                    = "89" # scans of this calibrator to use - e.g. "89" or "53", or "1"
#     pathOutputDirectory         = "../data/DDTB247/" # path to storage of non-CASA data products (bandpasses and DTECs) - e.g. "../data/DDTB247/" or "../data/31_008COSMOS/"
#     pathPlotDirectory           = "../figures/calibrationBandpass/DDTB247/" + nameField + "/" # e.g. "../figures/calibrationBandpass/DDTB247/" or "../figures/calibrationBandpass/31_008COSMOS/"
#
#     pathOutputCalTableCASA      = pathMS[ : -3]  + "-cal/bandpassGains"     + nameField + "_" + nameScan + ".table"
#     pathOutputCalTableNumPy     = pathOutputDirectory + "bandpassGains"     + nameField + "_" + nameScan + ".npy" # NumPy path to antenna-based gains
#     pathOutputBPsAmplitude      = pathOutputDirectory + "bandpassAmplitude" + nameField + "_" + nameScan + ".npy" # NumPy path to amplitude bandpasses
#     pathOutputBPsPhase          = pathOutputDirectory + "bandpassPhase"     + nameField + "_" + nameScan + ".npy" # NumPy path to phase     bandpasses
#     pathOutputDTECs             = pathOutputDirectory + "DTECs"             + nameField + "_" + nameScan + ".npy" # NumPy path to DTECs
#
#     nameTelescope               = "uGMRT"
#     numberOfFrequencyChannels   = 2048 # telescope property
#     timeStart                   = 0 # in seconds
#     frequencyRange              = 200 # in MHz
#     frequencyStart              = 300 # in MHz
#
#     # Initialise gain calculation parameters.
#     timeBinInterval             = 10 # user choice, in seconds
#     numberOfChannels       = 2048 # user choice, in 1
#     nameReferenceAntenna        = "C13" # user choice, e.g. "C00" or "C13"
#
#     # We normalise the frequency bandpasses using the following arbitrary criterion:
#     # the mean gain between 'frequencyNormalisationStart' and 'frequencyNormalisationEnd' must be 1.
#     frequencyNormalisationStart = 375 # in MHz
#     frequencyNormalisationEnd   = 450 # in MHz
#
#     # Initialise bandpass calculation algorithm settings.
#     numberOfSubIterationsBandpassAmplitude = 3
#     numberOfSubIterationsBandpassPhase     = 2
#     flaggingThresholdFactorAmplitude       = 3
#     flaggingThresholdFactorPhase           = 3
#
#
#     # Check whether the settings are valid.
#     if (numpy.mod(numberOfFrequencyChannels, numberOfChannels) == 0):
#         numberOfChannelsPerBin = numberOfFrequencyChannels / numberOfChannels
#     else:
#         print ("Error: the number of frequency channels per bin is not an integer!")
#         sys.exit()
#
#     # Sets plotting parameters.
#     polarisationNames      = ["1", "2"] # name of the two polarisation modes (influences plot titles and file names)
#     angleElevation         = 40 # in degrees (for 3D plotting)
#     angleAzimuthal         = -60 # in degrees (for 3D plotting)
#     plotTimeLimit          = 10 # in seconds
#     plotDerivTimeDTECLimit = .002 # in TECUs / s
#     plotDTECLimit          = .05 # in TECUs
#     plotFrequencyLimit     = 1 # in MHz
#     plotAmplitudeLimit     = .02 # in 1
#     plotPhaseLimit         = 5 # in degrees
#     plot                   = True # whether bandpass calibration plots should be made
#     plot3D                 = False # wether to create 3D plots; these require much more time to generate than 2D plots. Only applicable when 'plot' is True.
#
#
#     # Decide which steps to perform. Depending on the settings chosen, the bandpass calculation step can take a while.
#     performStepGenerateGains      = False
#     performStepPlotAmplitudes     = True
#     performStepPlotPhases         = False
#     performStepBandpassAmplitudes = False
#     performStepBandpassPhases     = False
#
#
#
#     # Calculate the bandpass using both time and frequency bins. Store the antenna-based complex gains in a NumPy array.
#     if (performStepGenerateGains):
#         # Load the packages here, instead of at the beginning, so that computers without CASANOVA can still run
#         # the plotting and bandpass calculation and visualisation procedures of this program.
#         import casac
#         tb = casac.casac.table()
#
#         # Calculate the bandpass using time-frequency bins. Store the result in a CASA Calibration Table.
#         '''
#         from casat import gaincal
#         gaincal = gaincal.gaincal
#         gaincal (vis = pathMS, caltable = pathOutputCalTableCASA, field = nameField, selectdata = True, scan = nameScan,
#                  solint = str(timeBinInterval) + "s, " + str(numberOfChannelsPerBin) + "ch", combine = "scan, field", refant = nameReferenceAntenna)
#         '''
#
#         # Running 'bandpass' with a finite time interval in 'solint' is actually equivalent to running 'gaincal'.
#         from casat import bandpass
#         bandpass = bandpass.bandpass
#         bandpass(vis = pathMS, caltable = pathOutputCalTableCASA, field = nameField, selectdata = True, scan = nameScan, solint = "int, " + str(numberOfChannelsPerBin) + "ch",
#                  combine = "scan, field", refant = nameReferenceAntenna)
#         #bandpass(vis = pathMS, caltable = pathOutputCalTableCASA, field = nameField, selectdata = True, scan = nameScan, solint = str(timeBinInterval) + "s, " + str(numberOfChannelsPerBin) + "ch",
#         #         combine = "scan, field", refant = nameReferenceAntenna)
#
#
#         # Open the Calibration Table, extract gains and flags, and close it again.
#         tb.open(pathOutputCalTableCASA)
#         gainsComplex = tb.getcol("CPARAM")
#         flags = tb.getcol("FLAG")
#         antennaIDs = tb.getcol("ANTENNA1")
#         tb.close()
#
#         # Infer the number of antennae from 'antennaIDs', by calculating the number of unique elements
#         # in the list with antennae IDs. Use CASA's 'browsetable()' to inspect the calibration table for more insight.
#         numberOfAntennae = len(numpy.unique(antenna1IDs))
#
#         # Infer the number of time bins from the data and the number of antennae.
#         numberOfTimeStamps = gainsComplex.shape[2] / numberOfAntennae
#
#         # Calculate the gain amplitudes and phases.
#         gainAmplitudes = numpy.absolute(gainsComplex) # in 1
#         gainPhases     = numpy.angle(gainsComplex, deg = True) # in degrees
#
#         # Divide up the data in two different polarizations.
#         gainAmplitudesPol1 = gainAmplitudes[0]
#         gainAmplitudesPol2 = gainAmplitudes[1]
#         gainPhasesPol1     = gainPhases[0]
#         gainPhasesPol2     = gainPhases[1]
#         flagsPol1          = flags[0]
#         flagsPol2          = flags[1]
#
#         # Reshape the gain amplitudes and phases and flags into a more practical fashion.
#         # (There is probably a way with fewer NumPy commands which I do not know.)
#         gainAmplitudesPol1Reshaped = []
#         gainAmplitudesPol2Reshaped = []
#         gainPhasesPol1Reshaped     = []
#         gainPhasesPol2Reshaped     = []
#         flagsPol1Reshaped          = []
#         flagsPol2Reshaped          = []
#
#         for i in range(numberOfAntennae):
#             # Create for the current antenna six 2D arrays, which will contain for each time - frequency pair either a gain amplitude or phase or a flag.
#             gainAmplitudesPol1Current = []
#             gainAmplitudesPol2Current = []
#             gainPhasesPol1Current     = []
#             gainPhasesPol2Current     = []
#             flagsPol1Current          = []
#             flagsPol2Current          = []
#
#             # Add a row to each of the arrays, which corresponds to a specific time.
#             for j in range(numberOfTimeStamps):
#                 gainAmplitudesPol1Current.append(gainAmplitudesPol1[ : , i + j * numberOfAntennae])
#                 gainAmplitudesPol2Current.append(gainAmplitudesPol2[ : , i + j * numberOfAntennae])
#                 gainPhasesPol1Current.append(    gainPhasesPol1    [ : , i + j * numberOfAntennae])
#                 gainPhasesPol2Current.append(    gainPhasesPol2    [ : , i + j * numberOfAntennae])
#                 flagsPol1Current.append(         flagsPol1         [ : , i + j * numberOfAntennae])
#                 flagsPol2Current.append(         flagsPol2         [ : , i + j * numberOfAntennae])
#
#             # NumPy-ify the arrays, transpose (flip time and frequency axes) and mirror the frequency axis.
#             gainAmplitudesPol1Current = numpy.flipud(numpy.transpose(numpy.array(gainAmplitudesPol1Current)))
#             gainAmplitudesPol2Current = numpy.flipud(numpy.transpose(numpy.array(gainAmplitudesPol2Current)))
#             gainPhasesPol1Current     = numpy.flipud(numpy.transpose(numpy.array(gainPhasesPol1Current)))
#             gainPhasesPol2Current     = numpy.flipud(numpy.transpose(numpy.array(gainPhasesPol2Current)))
#             flagsPol1Current          = numpy.flipud(numpy.transpose(numpy.array(flagsPol1Current)))
#             flagsPol2Current          = numpy.flipud(numpy.transpose(numpy.array(flagsPol2Current)))
#
#             # Add each of the final 2D arrays for this antenna to the right stack.
#             gainAmplitudesPol1Reshaped.append(gainAmplitudesPol1Current)
#             gainAmplitudesPol2Reshaped.append(gainAmplitudesPol2Current)
#             gainPhasesPol1Reshaped.append(    gainPhasesPol1Current)
#             gainPhasesPol2Reshaped.append(    gainPhasesPol2Current)
#             flagsPol1Reshaped.append(         flagsPol1Current)
#             flagsPol2Reshaped.append(         flagsPol2Current)
#
#         # NumPy-ify the arrays.
#         gainAmplitudesPol1Reshaped = numpy.array(gainAmplitudesPol1Reshaped)
#         gainAmplitudesPol2Reshaped = numpy.array(gainAmplitudesPol2Reshaped)
#         gainPhasesPol1Reshaped     = numpy.array(gainPhasesPol1Reshaped)
#         gainPhasesPol2Reshaped     = numpy.array(gainPhasesPol2Reshaped)
#         flagsPol1Reshaped          = numpy.array(flagsPol1Reshaped)
#         flagsPol2Reshaped          = numpy.array(flagsPol2Reshaped)
#
#         # Put together all data (for both polarisations), and save them to a NumPy file.
#         dataAll = numpy.array([gainAmplitudesPol1Reshaped, gainAmplitudesPol2Reshaped,
#                                gainPhasesPol1Reshaped,     gainPhasesPol2Reshaped,
#                                flagsPol1Reshaped,          flagsPol2Reshaped])
#         numpy.save(pathOutputCalTableNumPy, dataAll)
#
#         # Print diagnostic information.
#         print ("numberOfTimeStamps:\n", numberOfTimeStamps)
#         print ("gainsComplex.shape:\n", gainsComplex.shape)
#         print (gainAmplitudesPol1.shape,         gainAmplitudesPol2.shape)
#         print (gainAmplitudesPol1Reshaped.shape, gainAmplitudesPol2Reshaped.shape)
#
#
#
#     # Initialise data, flags and other variables for all subsequent steps.
#     if (performStepPlotAmplitudes or performStepPlotPhases or performStepBandpassAmplitudes or performStepBandpassPhases):
#         # Load data and flags from the NumPy file.
#         data = numpy.load(pathOutputCalTableNumPy)
#         gainAmplitudesPol1 = data[0]
#         gainAmplitudesPol2 = data[1]
#         gainPhasesPol1     = data[2] # in degrees, between -180 and 180
#         gainPhasesPol2     = data[3] # in degrees, between -180 and 180
#         flagsPol1          = data[4]
#         flagsPol2          = data[5]
#
#
#         # Flagged data should not be used in calculations.
#         # Masking using 'numpy.ma.masked_array(...)' is not always practical - the mask is lost during some NumPy operations.
#         # We choose to set flagged amplitudes and phases to 'numpy.nan'. This does not lead to desired colormap behaviour in
#         # 3D plotting. For that specific purpose, we use adapted grids where all flagged amplitudes and phases have value 0.
#         # If we don't alter the values of flagged amplitudes and phases, they will be crazy and this is annoying for plotting.
#         # We have the freedom to choose any value for flagged values, as they are meaningless anyways.
#         gainAmplitudesPol1 = numpy.where(flagsPol1, numpy.nan, gainAmplitudesPol1)
#         gainAmplitudesPol2 = numpy.where(flagsPol2, numpy.nan, gainAmplitudesPol2)
#         gainPhasesPol1     = numpy.where(flagsPol1, numpy.nan, gainPhasesPol1)
#         gainPhasesPol2     = numpy.where(flagsPol2, numpy.nan, gainPhasesPol2)
#
#         # 'numberOfAntennae' is defined in 'performStepGenerateGains', but might not exist here if that step has not been performed.
#         numberOfAntennae = data.shape[1]
#
#         # Create a list that indicates which antennae work.
#         antennaeWorking = []
#         for i in range(numberOfAntennae):
#             if (flagsPol1[i].all() and flagsPol2[i].all()):
#                 antennaeWorking.append(False)
#             else:
#                 antennaeWorking.append(True)
#
#         # Create 'gridFrequencies' and 'gridTimes', which are used for calculations and 3D plotting.
#         numberOfTimeStamps              = data.shape[3]
#         timeRange                     = timeBinInterval * numberOfTimeStamps # in seconds
#         gridIndexRow, gridIndexColumn = numpy.mgrid[0 : numberOfChannels, 0 : numberOfTimeStamps]
#         gridFrequencies               = numpy.divide(gridIndexRow, numberOfChannels - 1) * frequencyRange + frequencyStart # in MHz
#         gridTimes                     = numpy.divide(gridIndexColumn, numberOfTimeStamps - 1) * timeRange + timeStart # in seconds
#
#         # Print diagnostic information.
#         print ("data.shape:",                 data.shape) # E.g. (6, 30, 16, 67): 6 data types (amplitudes, phases, flags - for both pol.), 30 antennae, 16 frequency bins and 67 time bins
#         print ("gainAmplitudesPol1.shape:",   gainAmplitudesPol1.shape)
#         print ("gainPhasesPol1.shape:",       gainPhasesPol1.shape)
#         print ("gridIndexRow.shape:",         gridIndexRow.shape)
#         print ("antennaeWorking:",            antennaeWorking)
#         print ("number of working antennae:", numpy.sum(antennaeWorking))
#
#
#
#     if (performStepPlotAmplitudes):
#         # Create 2D and 3D plots of the antenna-based gain amplitudes for both polarisations.
#         for i in range(numberOfAntennae):
#             if (antennaeWorking[i]):
#                 print ("Starting gain amplitudes visualisation for antenna ID " + str(i) + "...")
#                 dataGrids = [gainAmplitudesPol1[i], gainAmplitudesPol2[i]]
#                 for dataGrid, polarisationName in zip(dataGrids, polarisationNames):
#                     # Create 2D antenna-based gain amplitudes plot.
#                     figure = pyplot.figure(figsize = (12, 6))
#                     image = pyplot.imshow(dataGrid,
#                                           aspect = "auto", interpolation = "nearest", cmap = cm.viridis, vmin = 0, vmax = 1)
#                     pyplot.xlabel("time (s)")
#                     pyplot.ylabel("frequency (MHz)")
#                     pyplot.xticks(numpy.linspace(0, numberOfTimeStamps - 1, num = 5, endpoint = True),
#                                   numpy.linspace(timeStart, timeStart + timeRange, num = 5, endpoint = True).astype(int))
#                     pyplot.yticks(numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
#                                   numpy.linspace(frequencyStart, frequencyStart + frequencyRange, num = 5, endpoint = True).astype(int))
#                     pyplot.title("antenna-based gain amplitudes of uncalibrated calibrator visibilities\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
#                     colorBar = pyplot.colorbar(image, cax = colorBarAxis, ticks = [0, 0.2, 0.4, 0.6, 0.8, 1])
#                     colorBar.ax.set_ylabel("gain amplitude ($1$)")
#                     pyplot.subplots_adjust(left = .06, right = .94, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "amplitudes2D_Ant" + str(i) + "_Pol" + polarisationName + ".pdf")
#                     pyplot.close()
#
#                     if (plot3D):
#                         # Create 3D antenna-based gain amplitudes plot. Because 3D plotting does not correctly
#                         # work with 'numpy.nan' (color artefacts appear), we replace all 'numpy.nan' with zeros.
#                         figure = pyplot.figure(figsize = (12, 6))
#                         axes3D = figure.add_subplot(111, projection = "3d")
#                         axes3D.plot_surface(gridFrequencies, gridTimes, numpy.where(numpy.isnan(dataGrid), 0, dataGrid),
#                                             vmin = 0, vmax = 1, cmap = cm.Spectral_r, alpha = .8, rstride = 1, cstride = 1)
#                         axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                         axes3D.set_xlabel("frequency (MHz)")
#                         axes3D.set_ylabel("time (s)")
#                         axes3D.set_zlabel("gain amplitude ($1$)")
#                         axes3D.set_zlim(0, 1)
#                         axes3D.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1])
#                         axes3D.set_title("antenna-based gain amplitudes of uncalibrated calibrator visibilities\ndata set: "
#                                      + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                      + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                         pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                         pyplot.savefig(pathPlotDirectory + "amplitudes3D_Ant" + str(i) + "_Pol" + polarisationName + ".png") # ".pdf" is possible here too - but 3D plots saved in PDF format are much bigger than those in PNG.
#                         pyplot.close()
#             else:
#                 print ("Skipping gain amplitudes visualisation for antenna ID " + str(i) + ": all data are flagged.")
#
#
#
#     if (performStepPlotPhases):
#         # Create 2D and 3D plots of the antenna-based gain phases for both polarisations.
#         for i in range(numberOfAntennae):
#             if (antennaeWorking[i]):
#                 print ("Starting gain phases visualisation for antenna ID " + str(i) + "...")
#                 dataGrids = [gainPhasesPol1[i], gainPhasesPol2[i]]
#                 for dataGrid, polarisationName in zip(dataGrids, polarisationNames):
#
#                     # Create 2D antenna-based gain phases plot.
#                     figure = pyplot.figure(figsize = (12, 6))
#                     image = pyplot.imshow(dataGrid,
#                                           aspect = "auto", interpolation = "none", cmap = cm.hsv, vmin = -180, vmax = 180)
#                     pyplot.xlabel("time (s)")
#                     pyplot.ylabel("frequency (MHz)")
#                     pyplot.xticks(numpy.linspace(0, numberOfTimeStamps - 1, num = 5, endpoint = True),
#                                   numpy.linspace(timeStart, timeStart + timeRange, num = 5, endpoint = True).astype(int))
#                     pyplot.yticks(numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
#                                   numpy.linspace(frequencyStart, frequencyStart + frequencyRange, num = 5, endpoint = True).astype(int))
#                     pyplot.title("antenna-based gain phases of uncalibrated calibrator visibilities\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "3%", pad = .05)
#                     colorBar = pyplot.colorbar(image, cax = colorBarAxis, ticks = [-180, -120, -60, 0, 60, 120, 180])
#                     colorBar.ax.set_ylabel("gain phase ($\degree$)")
#                     pyplot.subplots_adjust(left = .07, right = .93, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "phases2D_Ant" + str(i) + "_Pol" + polarisationName + ".pdf")
#                     pyplot.close()
#
#                     if (plot3D):
#                         # Create 3D antenna-based gain phases plot. Because 3D plotting does not correctly
#                         # work with 'numpy.nan' (color artefacts appear), we replace all 'numpy.nan' with zeros.
#                         figure = pyplot.figure(figsize = (12, 6))
#                         axes3D = figure.add_subplot(111, projection = "3d")
#                         axes3D.plot_surface(gridFrequencies, gridTimes, numpy.where(numpy.isnan(dataGrid), 0, dataGrid),
#                                             vmin = -180, vmax = 180, cmap = cm.YlOrRd, alpha = .8, rstride = 1, cstride = 1)
#                         axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                         axes3D.set_xlabel("frequency (MHz)")
#                         axes3D.set_ylabel("time (s)")
#                         axes3D.set_zlabel("gain phase ($\degree$)")
#                         axes3D.set_zlim(-180, 180)
#                         axes3D.set_zticks(numpy.array([-180, -120, -60, 0, 60, 120, 180]).astype(int))
#                         axes3D.set_title("antenna-based gain phases of uncalibrated calibrator visibilities\ndata set: "
#                                      + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                      + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                         pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                         pyplot.savefig(pathPlotDirectory + "phases3D_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                         pyplot.close()
#             else:
#                 print ("Skipping gain phases visualisation for antenna ID " + str(i) + ": all data are flagged.")
#
#
#
#     if (performStepBandpassAmplitudes):
#         '''
#         In this step, the amplitude bandpasses for each antenna are determined in two iterations.
#         Iteration 1 is subdivided into several 'subiterations' that are meant to flag outliers.
#         Iteration 2 takes the bandpass results from iteration 1 and performs median filtering and interpolation.
#         '''
#         # Create lists that store the normalised amplitude bandpasses and the bandpass normalisation factors.
#         bandpassesAmplitudePol1 = []
#         bandpassesAmplitudePol2 = []
#         bandpassNormalisationFactorsPol1 = []
#         bandpassNormalisationFactorsPol2 = []
#
#         for i in range(numberOfAntennae):
#             if (antennaeWorking[i]):
#                 print ("Starting amplitude bandpass calculation for antenna ID " + str(i) + "...")
#
#                 for subIteration in range(numberOfSubIterationsBandpassAmplitude):
#                     # Determine the amplitude bandpass (iteration 1).
#                     if (subIteration == numberOfSubIterationsBandpassAmplitude - 1): # If in the last iteration, determine the amplitude bandpass using the mean.
#                         bandpassAmplitudePol1Iter1 = numpy.nanmean(gainAmplitudesPol1[i], axis = 1)
#                         bandpassAmplitudePol2Iter1 = numpy.nanmean(gainAmplitudesPol2[i], axis = 1)
#                     else: # During earlier iterations, to avoid e.g. RFI signatures, we use the median.
#                         bandpassAmplitudePol1Iter1 = numpy.nanmedian(gainAmplitudesPol1[i], axis = 1)
#                         bandpassAmplitudePol2Iter1 = numpy.nanmedian(gainAmplitudesPol2[i], axis = 1)
#
#                     # Tile the bandpass into a grid.
#                     gridBandpassAmplitudePol1 = numpy.tile(bandpassAmplitudePol1Iter1, (numberOfTimeStamps, 1)).T
#                     gridBandpassAmplitudePol2 = numpy.tile(bandpassAmplitudePol2Iter1, (numberOfTimeStamps, 1)).T
#
#                     # Remove the bandpass from the data. Residuals will be centered around 1.
#                     gridResidualsPol1 = numpy.divide(gainAmplitudesPol1[i], gridBandpassAmplitudePol1)
#                     gridResidualsPol2 = numpy.divide(gainAmplitudesPol2[i], gridBandpassAmplitudePol2)
#
#                     # Calculate the standard deviation of the residuals.
#                     STDPol1 = numpy.nanstd(gridResidualsPol1)
#                     STDPol2 = numpy.nanstd(gridResidualsPol2)
#
#                     # Determine outliers and update the flags. This makes sure that the bandpass in the next iterations is better.
#                     gridIsOutlierPol1 = numpy.greater(numpy.absolute(gridResidualsPol1 - 1), flaggingThresholdFactorAmplitude * STDPol1)
#                     gridIsOutlierPol2 = numpy.greater(numpy.absolute(gridResidualsPol2 - 1), flaggingThresholdFactorAmplitude * STDPol2)
#                     flagsPol1[i] = numpy.logical_or(flagsPol1[i], gridIsOutlierPol1)
#                     flagsPol2[i] = numpy.logical_or(flagsPol2[i], gridIsOutlierPol2)
#
#                     # Set flagged data to 'numpy.nan'.
#                     gainAmplitudesPol1[i] = numpy.where(flagsPol1[i], numpy.nan, gainAmplitudesPol1[i])
#                     gainAmplitudesPol2[i] = numpy.where(flagsPol2[i], numpy.nan, gainAmplitudesPol2[i])
#                     gainPhasesPol1[i]     = numpy.where(flagsPol1[i], numpy.nan, gainPhasesPol1[i])
#                     gainPhasesPol2[i]     = numpy.where(flagsPol2[i], numpy.nan, gainPhasesPol2[i])
#
#                 # Median filter and interpolate the amplitude bandpass (iteration 2).
#                 bandpassAmplitudePol1Iter2 = fillGaps1D(filterMedian1D(bandpassAmplitudePol1Iter1, kernelSize = 7))
#                 bandpassAmplitudePol2Iter2 = fillGaps1D(filterMedian1D(bandpassAmplitudePol2Iter1, kernelSize = 7))
#
#                 '''
#                 # Use this piece of code to debug the second iteration of the amplitude bandpass.
#                 pyplot.plot(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter1, c = 'r', linestyle = "--")
#                 pyplot.plot(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter2, c = 'b', linestyle = ":")
#                 pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter1, s = 4, c = 'r')
#                 pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter2, s = 4, c = 'b')
#                 pyplot.show()
#                 sys.exit()
#                 '''
#
#                 # Normalise the amplitude bandpass for both polarisations.
#                 # We do not normalise by determining the factor that scales the amplitude bandpass peak to 1,
#                 # as this would result in wrong scaling behaviour in cases where the 'real' peak is absent in the
#                 # bandpass due to flagging.
#                 # We also do not normalise by determining the factor that would scale the mean of the bandpass
#                 # to 1, as this requires taking the mean of all channels except for those flagged. In such case,
#                 # when many channels at the lower and upper end of the frequency band are non-flagged, the scaling
#                 # is different than when many of those edge channels are flagged.
#                 # We conclude that the best thing to do is to normalise by determining the factor that would scale
#                 # the mean bandpass in a selected range of frequency channels to 1. These channels are chosen such
#                 # that the bandpass is approximately constant over them.
#                 # If we normalise bandpass iteration 2, then we could use 'numpy.mean' instead of 'numpy.nanmean' as well
#                 # in the determination of 'bandpassNormalisationFactorPolX'. After all, due to the linear interpolation
#                 # only the edges of the frequency band still lack a bandpass value. We still prefer using 'numpy.nanmean'
#                 # for antennae with a lot of flagged channels at the outer parts of the frequency band, as in such cases
#                 # a few 'numpy.nan's might lie into the domain '[frequencyNormalisationStart, frequencyNormalisationEnd]'.
#
#
#                 # Determine the normalisation factors.
#                 indexStart = numpy.argmax(gridFrequencies[ : , 0] > frequencyNormalisationStart)
#                 indexEnd   = numpy.argmax(gridFrequencies[ : , 0] > frequencyNormalisationEnd)
#                 bandpassNormalisationFactorPol1 = numpy.nanmean(bandpassAmplitudePol1Iter2[indexStart : indexEnd]) # 'numpy.mean' would work too, mostly
#                 bandpassNormalisationFactorPol2 = numpy.nanmean(bandpassAmplitudePol2Iter2[indexStart : indexEnd]) # 'numpy.mean' would work too, mostly
#
#                 # Divide the amplitude bandpasses by the normalisation factor.
#                 bandpassAmplitudePol1Iter1 = numpy.divide(bandpassAmplitudePol1Iter1, bandpassNormalisationFactorPol1)
#                 bandpassAmplitudePol2Iter1 = numpy.divide(bandpassAmplitudePol2Iter1, bandpassNormalisationFactorPol2)
#                 bandpassAmplitudePol1Iter2 = numpy.divide(bandpassAmplitudePol1Iter2, bandpassNormalisationFactorPol1)
#                 bandpassAmplitudePol2Iter2 = numpy.divide(bandpassAmplitudePol2Iter2, bandpassNormalisationFactorPol2)
#
#                 # Add the normalised amplitude bandpasses and normalisation factors to the appropriate lists.
#                 bandpassesAmplitudePol1.append(bandpassAmplitudePol1Iter2)
#                 bandpassesAmplitudePol2.append(bandpassAmplitudePol2Iter2)
#                 bandpassNormalisationFactorsPol1.append(bandpassNormalisationFactorPol1)
#                 bandpassNormalisationFactorsPol2.append(bandpassNormalisationFactorPol2)
#
#
#                 if (plot):
#                     print ("Starting amplitude bandpass visualisation for antenna ID " + str(i) + "...")
#
#                     # Create plot of amplitude bandpass (for both polarisations, iteration 1).
#                     pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter1, c = "navy", s = 16, lw = 0, label = "polarisation 1\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol1, 3)))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol2Iter1, c = "orangered", s = 16, lw = 0, label = "polarisation 2\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol2, 3)))
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.xlabel("frequency channel centre (MHz)")
#                     pyplot.ylabel("antenna-based gain amplitude (1)")
#                     pyplot.xlim(frequencyStart - plotFrequencyLimit, frequencyStart + frequencyRange + plotFrequencyLimit)
#                     pyplot.ylim(0, 2 + plotAmplitudeLimit)
#                     pyplot.title("normalised amplitude bandpass (iteration 1)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .07, right = .98, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "bandpassAmplitudeIter1_" + "Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     # Create plot of amplitude bandpass (for both polarisations, iteration 2).
#                     pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol1Iter2, c = "navy", s = 16, lw = 0, label = "polarisation 1\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol1, 3)))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassAmplitudePol2Iter2, c = "orangered", s = 16, lw = 0, label = "polarisation 2\nnorm. factor: " + str(numpy.round(bandpassNormalisationFactorPol2, 3)))
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.xlabel("frequency channel centre (MHz)")
#                     pyplot.ylabel("antenna-based gain amplitude (1)")
#                     pyplot.xlim(frequencyStart - plotFrequencyLimit, frequencyStart + frequencyRange + plotFrequencyLimit)
#                     pyplot.ylim(0, 2 + plotAmplitudeLimit)
#                     pyplot.title("normalised amplitude bandpass (iteration 2)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .07, right = .98, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "bandpassAmplitudeIter2_" + "Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     if (plot3D):
#                         dataGrids = [gridBandpassAmplitudePol1, gridBandpassAmplitudePol2]
#                         for dataGrid, polarisationName in zip(dataGrids, polarisationNames):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             axes3D.plot_surface(gridFrequencies, gridTimes, numpy.where(numpy.isnan(dataGrid), 0, dataGrid), cmap = cm.Spectral_r, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("amplitude ($1$)")
#                             axes3D.set_zlim(0, 1)
#                             axes3D.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1])
#                             axes3D.set_title("unnormalised amplitude bandpass (iteration 1)\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "bandpassAmplitude3DIter1_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#
#                         # Create 3D plots of antenna-based gain amplitude residuals.
#                         dataGrids = [gridResidualsPol1, gridResidualsPol2]
#                         STDs = [STDPol1, STDPol2]
#                         for dataGrid, STD, polarisationName in zip(dataGrids, STDs, polarisationNames):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             axes3D.plot_surface(gridFrequencies, gridTimes, numpy.where(numpy.isnan(dataGrid), 1, dataGrid), cmap = cm.YlOrRd, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("amplitude (1)")
#                             axes3D.set_zlim(1 - .1, 1 + .1)
#                             axes3D.set_title("antenna-based gain amplitude residuals (iteration 1) | $\sigma$ = " + str(numpy.round(STD, 3)) + "\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "residualsAmplitudeIter1_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#             else:
#                 print ("Skipping amplitude bandpass calculation for antenna ID " + str(i) + ": all data are flagged.")
#
#         # Save the normalised amplitude bandpasses and their normalisation factors, after all working antennae have been looped over.
#         numpy.save(pathOutputBPsAmplitude, [bandpassesAmplitudePol1, bandpassesAmplitudePol2, bandpassNormalisationFactorsPol1, bandpassNormalisationFactorsPol2, antennaeWorking])
#
#
#     if (performStepBandpassPhases):
#         '''
#         In this step, the phase bandpass is found in an iterative process. DTEC(t) (for the calibrator)
#         is found as a side product.
#         '''
#
#         # Create a 2D and a 1D array of inverse frequencies over the frequency band.
#         gridFrequenciesInverse = numpy.divide(1, gridFrequencies) # in MHz^-1
#         frequenciesInverse     = gridFrequenciesInverse[ : , 0] # in MHz^-1
#
#         # Create lists that will contain, for each antenna, the two phase bandpasses and the function DTEC(t).
#         bandpassesPhasePol1    = []
#         bandpassesPhasePol2    = []
#         DTECsList              = []
#
#         for i in range(numberOfAntennae):
#             if (antennaeWorking[i]):
#
#                 print ("Starting DTEC and phase bandpass calculation for antenna ID " + str(i) + "...")
#
#                 # Calculate the first derivative of gain phase to time in a way robust to phase wrapping (for both polarisations).
#                 # We do so by calculating the derivative for each time-frequency bin 2 times: one with the ordinary data, and once after shifting
#                 # - all the phases by 180 degrees to place them in the [0, 360) domain, and
#                 # - another 180 degrees to effect a shift within that domain.
#                 derivTimePol1Choice1        = computeDerivative2D(gainPhasesPol1[i], timeBinInterval, axis = 1, degree = 1, intermediateSampling = True) # in degrees per second
#                 derivTimePol1Choice2        = computeDerivative2D(numpy.mod(gainPhasesPol1[i] + 180 + 180, 360), timeBinInterval, axis = 1, degree = 1, intermediateSampling = True) # in degrees per second
#                 derivTimePol1               = numpy.where(numpy.less(numpy.absolute(derivTimePol1Choice1), numpy.absolute(derivTimePol1Choice2)),
#                                                           derivTimePol1Choice1, derivTimePol1Choice2) # in degrees per MHz^2
#
#                 derivTimePol2Choice1        = computeDerivative2D(gainPhasesPol2[i], timeBinInterval, axis = 1, degree = 1, intermediateSampling = True) # in degrees per second
#                 derivTimePol2Choice2        = computeDerivative2D(numpy.mod(gainPhasesPol2[i] + 180 + 180, 360), timeBinInterval, axis = 1, degree = 1, intermediateSampling = True) # in degrees per second
#                 derivTimePol2               = numpy.where(numpy.less(numpy.absolute(derivTimePol2Choice1), numpy.absolute(derivTimePol2Choice2)),
#                                                           derivTimePol2Choice1, derivTimePol2Choice2) # in degrees per MHz^2
#
#
#                 # Determine the DTEC time derivative for both polarisations (iteration 1).
#                 derivTimeDTECsIter1Pol1Grid = derivTimePol1 * gridFrequencies[ : , : -1] / (1210 * 400)
#                 derivTimeDTECsIter1Pol1     = numpy.nanmean(derivTimeDTECsIter1Pol1Grid, axis = 0)
#
#                 derivTimeDTECsIter1Pol2Grid = derivTimePol2 * gridFrequencies[ : , : -1] / (1210 * 400)
#                 derivTimeDTECsIter1Pol2     = numpy.nanmean(derivTimeDTECsIter1Pol2Grid, axis = 0)
#
#
#                 #'''
#                 # Create a grid of residuals. CURRENTLY NOT USED IN THE ALGORITHM!
#                 derivTimeFitPol1            = numpy.multiply(gridFrequenciesInverse[ : , : -1], numpy.tile(derivTimeDTECsIter1Pol1, (numberOfChannels, 1))) * 1210 * 400 # in degrees per second
#                 derivTimeFitPol2            = numpy.multiply(gridFrequenciesInverse[ : , : -1], numpy.tile(derivTimeDTECsIter1Pol2, (numberOfChannels, 1))) * 1210 * 400 # in degrees per second
#                 derivTimeResidualsPol1      = derivTimePol1 - derivTimeFitPol1
#                 derivTimeResidualsPol2      = derivTimePol2 - derivTimeFitPol2
#
#                 # To do: flag or mask the worst outliers, and fit 'derivTimeDTECsIter1PolX' again with the new flags.
#                 #'''
#
#
#                 # Whenever there are NaNs in 'derivTimeDTECsIter1PolX', then 'DTECsIter1PolX' is not well-determined.
#                 # If NaNs occur at the beginning or end, they are replaced by zeros.
#                 derivTimeDTECsIter1Pol1     = fillGaps1D(derivTimeDTECsIter1Pol1)
#                 derivTimeDTECsIter1Pol2     = fillGaps1D(derivTimeDTECsIter1Pol2)
#
#
#                 # Find DTECs by integrating along time for both polarisations (iteration 1).
#                 DTECsIter1Pol1              = [0]
#                 DTECsIter1Pol2              = [0]
#                 for j in range(numberOfTimeStamps - 1):
#                     DTECsIter1Pol1.append(timeBinInterval * numpy.nansum(derivTimeDTECsIter1Pol1[ : j + 1]))
#                     DTECsIter1Pol2.append(timeBinInterval * numpy.nansum(derivTimeDTECsIter1Pol2[ : j + 1]))
#                 '''
#                 # We miss a DTEC-determination for the last time slice, as the calculation of a time derivative relies
#                 # on the availability of data at both the previous and the next time slice.
#                 # For easy array manipulation, we artificially add a DTEC for the last time slice.
#                 # On the basis of time coherency, we add a copy of the last value.
#                 # In this way, the lists contain 'numberOfTimeStamps' DTECs.
#                 # Note: this makes the last time slice of the gain phases to be unreliable!
#                 DTECsIter1Pol1.append(DTECsIter1Pol1[-1])
#                 DTECsIter1Pol2.append(DTECsIter1Pol2[-1])
#                 '''
#                 # NumPy-ify the 1D arrays.
#                 DTECsIter1Pol1              = numpy.array(DTECsIter1Pol1)
#                 DTECsIter1Pol2              = numpy.array(DTECsIter1Pol2)
#
#
#                 # Flag DTEC values appropriately, or do median filtering?
#                 '''
#                 '''
#
#                 # Subtract the mean or median from the DTEC values, so that they lie around 0.
#                 DTECsIter1Pol1             -= numpy.mean(DTECsIter1Pol1) #DTECsIter1Pol1[ : -1]
#                 DTECsIter1Pol2             -= numpy.mean(DTECsIter1Pol2) #DTECsIter1Pol2[ : -1]
#
#
#                 # Calculate the fitted plasma opacity phase effect (iteration 1).
#                 gridPlasmaOpIter1Pol1       = wrapPhasesCenter0(numpy.multiply(
#                                               gridFrequenciesInverse, numpy.tile(DTECsIter1Pol1, (numberOfChannels, 1))) * 1210 * 400) # in degrees
#                 gridPlasmaOpIter1Pol2       = wrapPhasesCenter0(numpy.multiply(
#                                               gridFrequenciesInverse, numpy.tile(DTECsIter1Pol2, (numberOfChannels, 1))) * 1210 * 400) # in degrees
#
#                 # We now determine the bandpass by subtracting the plasma opacity phase effect from the original data.
#                 # As we flag outliers in the process, we repeat this procedure several times.
#                 # During the last subiteration, the bandpass is calculated by taking the mean along time, whereas in
#                 # earlier iterations we calculate the bandpass using the median.
#                 for subIteration in range(numberOfSubIterationsBandpassPhase):
#                     # Subtract the time-variable ionospheric effect from the gain phases.
#                     gridStaticIter1Pol1 = wrapPhasesCenter0(gainPhasesPol1[i] - gridPlasmaOpIter1Pol1)#numpy.ma.masked_array(wrapPhasesCenter0(gainPhasesPol1[i] - gridPlasmaOpIter1Pol1), flagsPol1[i])
#                     gridStaticIter1Pol2 = wrapPhasesCenter0(gainPhasesPol2[i] - gridPlasmaOpIter1Pol2)#numpy.ma.masked_array(wrapPhasesCenter0(gainPhasesPol2[i] - gridPlasmaOpIter1Pol2), flagsPol2[i])
#
#
#                     # Calculate the phase bandpass for both polarisations in a way robust to phase wrapping (iteration 1).
#                     # This is again done by creating two candidate means where one is found by first shifting the phases by 180 degrees.
#                     # Then the standard deviations for both mean candidates are calculated. The mean candidate is chosen that has the lowest standard deviation.
#                     if (subIteration == numberOfSubIterationsBandpassPhase - 1): # If in the last iteration, determine the phase bandpass using the mean.
#                         bandpassPhaseIter1Pol1Choice1 = numpy.nanmean(gridStaticIter1Pol1, axis = 1)
#                         bandpassPhaseIter1Pol2Choice1 = numpy.nanmean(gridStaticIter1Pol2, axis = 1)
#                         bandpassPhaseIter1Pol1Choice2 = wrapPhasesCenter0(numpy.nanmean(numpy.mod(gridStaticIter1Pol1 + 180 + 180, 360), axis = 1) - 180 - 180)
#                         bandpassPhaseIter1Pol2Choice2 = wrapPhasesCenter0(numpy.nanmean(numpy.mod(gridStaticIter1Pol2 + 180 + 180, 360), axis = 1) - 180 - 180)
#                     else: # During earlier iterations, to avoid RFI signatures, we use the median.
#                         bandpassPhaseIter1Pol1Choice1 = numpy.nanmedian(gridStaticIter1Pol1, axis = 1)
#                         bandpassPhaseIter1Pol2Choice1 = numpy.nanmedian(gridStaticIter1Pol2, axis = 1)
#                         bandpassPhaseIter1Pol1Choice2 = wrapPhasesCenter0(numpy.nanmedian(numpy.mod(gridStaticIter1Pol1 + 180 + 180, 360), axis = 1) - 180 - 180)
#                         bandpassPhaseIter1Pol2Choice2 = wrapPhasesCenter0(numpy.nanmedian(numpy.mod(gridStaticIter1Pol2 + 180 + 180, 360), axis = 1) - 180 - 180)
#
#
#                     bandpassPhaseErrorIter1Pol1Choice1 = numpy.nanstd(gridStaticIter1Pol1, axis = 1)
#                     bandpassPhaseErrorIter1Pol2Choice1 = numpy.nanstd(gridStaticIter1Pol2, axis = 1)
#                     bandpassPhaseErrorIter1Pol1Choice2 = numpy.nanstd(numpy.mod(gridStaticIter1Pol1 + 180 + 180, 360), axis = 1)
#                     bandpassPhaseErrorIter1Pol2Choice2 = numpy.nanstd(numpy.mod(gridStaticIter1Pol2 + 180 + 180, 360), axis = 1)
#
#                     # Calculate the phase bandpasses. NOTE: although these are called 'bandpassPhaseIter1PolX', we can remove 'Iter1'!
#                     # We won't update the bandpass anymore after this, so no need for the extra tag 'Iter1'.
#                     bandpassPhaseIter1Pol1      = numpy.where(numpy.less(bandpassPhaseErrorIter1Pol1Choice1, bandpassPhaseErrorIter1Pol1Choice2),
#                                                               bandpassPhaseIter1Pol1Choice1, bandpassPhaseIter1Pol1Choice2) # in degrees
#                     bandpassPhaseIter1Pol2      = numpy.where(numpy.less(bandpassPhaseErrorIter1Pol2Choice1, bandpassPhaseErrorIter1Pol2Choice2),
#                                                               bandpassPhaseIter1Pol2Choice1, bandpassPhaseIter1Pol2Choice2) # in degrees
#
#                     #print(numpy.sum(numpy.isnan(bandpassPhaseIter1Pol1)))
#
#                     # Remove the phase bandpass from the static pattern to keep residuals, which we will use for flagging.
#                     gridResidualsIter1Pol1      = wrapPhasesCenter0(gridStaticIter1Pol1 - numpy.tile(bandpassPhaseIter1Pol1, (numberOfTimeStamps, 1)).T)
#                     gridResidualsIter1Pol2      = wrapPhasesCenter0(gridStaticIter1Pol2 - numpy.tile(bandpassPhaseIter1Pol2, (numberOfTimeStamps, 1)).T)
#
#
#                     # Calculate the standard deviation of the residuals.
#                     STDIter1Pol1                = numpy.nanstd(gridResidualsIter1Pol1)
#                     STDIter1Pol2                = numpy.nanstd(gridResidualsIter1Pol2)
#
#                     # Determine outliers and update the flags. This makes sure that the bandpass in the next iterations is better.
#                     gridIsOutlierPol1           = numpy.greater(numpy.absolute(gridResidualsIter1Pol1), flaggingThresholdFactorPhase * STDIter1Pol1)
#                     gridIsOutlierPol2           = numpy.greater(numpy.absolute(gridResidualsIter1Pol2), flaggingThresholdFactorPhase * STDIter1Pol2)
#                     flagsPol1[i]                = numpy.logical_or(flagsPol1[i], gridIsOutlierPol1)
#                     flagsPol2[i]                = numpy.logical_or(flagsPol2[i], gridIsOutlierPol2)
#
#                     # Set flagged data to 'numpy.nan'.
#                     gainAmplitudesPol1[i]       = numpy.where(flagsPol1[i], numpy.nan, gainAmplitudesPol1[i])
#                     gainAmplitudesPol2[i]       = numpy.where(flagsPol2[i], numpy.nan, gainAmplitudesPol2[i])
#                     gainPhasesPol1[i]           = numpy.where(flagsPol1[i], numpy.nan, gainPhasesPol1[i])
#                     gainPhasesPol2[i]           = numpy.where(flagsPol2[i], numpy.nan, gainPhasesPol2[i])
#
#
#
#                 # Remove the phase bandpass from the original data.
#                 gridIonosOnlyPol1               = wrapPhasesCenter0(gainPhasesPol1[i] - numpy.tile(bandpassPhaseIter1Pol1, (numberOfTimeStamps, 1)).T)
#                 gridIonosOnlyPol2               = wrapPhasesCenter0(gainPhasesPol2[i] - numpy.tile(bandpassPhaseIter1Pol2, (numberOfTimeStamps, 1)).T)
#
#                 # Calculate DTECs directly by fitting to the ionospheric data only (iteration 2).
#                 DTECsIter2Pol1                  = numpy.nanmean(gridIonosOnlyPol1 * gridFrequencies / (1210 * 400), axis = 0)#numpy.linalg.lstsq(matrix, numpy.ma.masked_array(gridIonosOnlyPol1, flagsPol1[i]))[0][0] # in TECUs
#                 DTECsIter2Pol2                  = numpy.nanmean(gridIonosOnlyPol2 * gridFrequencies / (1210 * 400), axis = 0)#numpy.linalg.lstsq(matrix, numpy.ma.masked_array(gridIonosOnlyPol2, flagsPol2[i]))[0][0] # in TECUs
#
#
#                 # Calculate the fitted plasma opacity phase effect (iteration 2).
#                 gridPlasmaOpIter2Pol1           = wrapPhasesCenter0(numpy.multiply(
#                                                   gridFrequenciesInverse, numpy.tile(DTECsIter2Pol1, (numberOfChannels, 1))) * 1210 * 400) # in degrees
#                 gridPlasmaOpIter2Pol2           = wrapPhasesCenter0(numpy.multiply(
#                                                   gridFrequenciesInverse, numpy.tile(DTECsIter2Pol2, (numberOfChannels, 1))) * 1210 * 400) # in degrees
#
#
#                 # Subtract the time-variable ionospheric effect from the gain phases.
#                 gridStaticIter2Pol1             = wrapPhasesCenter0(gainPhasesPol1[i] - gridPlasmaOpIter2Pol1)#numpy.ma.masked_array(wrapPhasesCenter0(gainPhasesPol1[i] - gridPlasmaOpIter2Pol1), flagsPol1[i])
#                 gridStaticIter2Pol2             = wrapPhasesCenter0(gainPhasesPol2[i] - gridPlasmaOpIter2Pol2)#numpy.ma.masked_array(wrapPhasesCenter0(gainPhasesPol2[i] - gridPlasmaOpIter2Pol2), flagsPol2[i])
#
#                 # Remove the phase bandpass from the static pattern to keep residuals, which we will use for flagging.
#                 gridResidualsIter2Pol1          = wrapPhasesCenter0(gridStaticIter2Pol1 - numpy.tile(bandpassPhaseIter1Pol1, (numberOfTimeStamps, 1)).T)
#                 gridResidualsIter2Pol2          = wrapPhasesCenter0(gridStaticIter2Pol2 - numpy.tile(bandpassPhaseIter1Pol2, (numberOfTimeStamps, 1)).T)
#
#                 # Calculate the standard deviation of the residuals.
#                 STDIter2Pol1                    = numpy.nanstd(gridResidualsIter2Pol1)
#                 STDIter2Pol2                    = numpy.nanstd(gridResidualsIter2Pol2)
#
#
#                 # Calculate the derivative of the phase bandpass to frequency (for both polarisations).
#                 frequencyBinInterval            = gridFrequencies[1, 0] - gridFrequencies[0, 0]
#                 bandpassPhaseDerivPol1          = computeDerivative1D(bandpassPhaseIter1Pol1, stepSize = frequencyBinInterval, intermediateSampling = True)
#                 bandpassPhaseDerivPol2          = computeDerivative1D(bandpassPhaseIter1Pol2, stepSize = frequencyBinInterval, intermediateSampling = True)
#
#                 # Determine the mean, median and standard deviation of the derivative after sigma clipping (for both polarisations).
#                 meanPol1, medianPol1, sigmaPol1 = libraryRadioAstronomy.sigmaClip1D(bandpassPhaseDerivPol1, verbose = False) # in degrees per MHz
#                 meanPol2, medianPol2, sigmaPol2 = libraryRadioAstronomy.sigmaClip1D(bandpassPhaseDerivPol2, verbose = False) # in degrees per MHz
#
#                 # Determine whether the mean or median should be used for generating the bandpass fit (for both polarisations).
#                 if (sigmaPol1 < 10): # in degrees per MHz
#                     slopePol1 = meanPol1
#                     print("Phase bandpass fit | slope type: mean")
#                 else:
#                     slopePol1 = medianPol1
#                     print("Phase bandpass fit | slope type: median")
#                 if (sigmaPol2 < 10): # in degrees per MHz
#                     slopePol2 = meanPol2
#                     print("Phase bandpass fit | slope type: mean")
#                 else:
#                     slopePol2 = medianPol2
#                     print("Phase bandpass fit | slope type: median")
#
#                 # Calculate the bandpass fit (for both polarisations).
#                 BPPhaseFitPol1 = []
#                 BPPhaseFitPol2 = []
#
#                 for indexRefPol1 in range(numberOfChannels):
#                     if (not numpy.isnan(bandpassPhaseIter1Pol1[indexRefPol1])):
#                         break
#                 for indexRefPol2 in range(numberOfChannels):
#                     if (not numpy.isnan(bandpassPhaseIter1Pol2[indexRefPol2])):
#                         break
#
#                 print("indexRefPol1, indexRefPol2:", indexRefPol1, indexRefPol2)
#
#                 for j in range(numberOfChannels):
#                     BPPhaseFitPol1.append(bandpassPhaseIter1Pol1[indexRefPol1] + frequencyBinInterval * (j - indexRefPol1) * slopePol1)
#                     BPPhaseFitPol2.append(bandpassPhaseIter1Pol2[indexRefPol2] + frequencyBinInterval * (j - indexRefPol2) * slopePol2)
#
#                 # NumPy-ify the 1D arrays, and ensure phase wrapping.
#                 BPPhaseFitPol1 = wrapPhasesCenter0(numpy.array(BPPhaseFitPol1))
#                 BPPhaseFitPol2 = wrapPhasesCenter0(numpy.array(BPPhaseFitPol2))
#
#
#                 # Calculate the residuals that remain after the fit has been subtracted from the data.
#                 BPPhaseRes1Pol1 = wrapPhasesCenter0(bandpassPhaseIter1Pol1 - BPPhaseFitPol1)
#                 BPPhaseRes1Pol2 = wrapPhasesCenter0(bandpassPhaseIter1Pol2 - BPPhaseFitPol2)
#
#                 # We hope we can now apply median filtering and interpolation, as we removed the phase ramps.
#                 BPPhaseRes2Pol1 = fillGaps1D(filterMedian1D(BPPhaseRes1Pol1, kernelSize = 7))
#                 BPPhaseRes2Pol2 = fillGaps1D(filterMedian1D(BPPhaseRes1Pol2, kernelSize = 7))
#
#                 # Now add back in the phase ramps.
#                 bandpassPhaseIter2Pol1 = wrapPhasesCenter0(BPPhaseRes2Pol1 + BPPhaseFitPol1)
#                 bandpassPhaseIter2Pol2 = wrapPhasesCenter0(BPPhaseRes2Pol2 + BPPhaseFitPol2)
#
#                 '''
#                 BPPhaseFitPol1 = [bandpassPhaseIter1Pol1[0]]
#                 BPPhaseFitPol2 = [bandpassPhaseIter1Pol2[0]]
#                 for j in range(numberOfChannels - 1):
#                     BPPhaseFitPol1.append(BPPhaseFitPol1[0] + frequencyBinInterval * (j + 1) * slopePol1)
#                     BPPhaseFitPol2.append(BPPhaseFitPol2[0] + frequencyBinInterval * (j + 1) * slopePol2)
#                 '''
#                 '''
#                 pyplot.scatter(range(len(bandpassPhaseDerivPol1)), bandpassPhaseDerivPol1,s = 4)
#                 pyplot.scatter(range(len(bandpassPhaseDerivPol2)), bandpassPhaseDerivPol2, s= 4)
#                 pyplot.axhline(y = medianPol1)
#                 pyplot.axhline(y = medianPol2)
#                 pyplot.axhline(y = meanPol1, ls = "--")
#                 pyplot.axhline(y = meanPol2, ls = "--")
#                 pyplot.axhline(y = medianPol1 - sigmaPol1, ls = ":", c = 'b')
#                 pyplot.axhline(y = medianPol1 + sigmaPol1, ls = ":", c = 'b')
#                 pyplot.show()
#                 '''
#                 '''
#                 pyplot.scatter(range(len(BPPhaseRes1Pol1)), BPPhaseRes1Pol1, s= 4)
#                 pyplot.scatter(range(len(BPPhaseRes1Pol2)), BPPhaseRes1Pol2, s=4)
#                 pyplot.show()
#                 '''
#
#                 # Append the final phase bandpasses to their respective lists, which will be saved to disk.
#                 bandpassesPhasePol1.append(bandpassPhaseIter2Pol1)
#                 bandpassesPhasePol2.append(bandpassPhaseIter2Pol2)
#
#                 # Append the mean of the two DTEC time series to the final list with DTECs.
#                 # The data is averaged because both polarisations should measure the same DTEC.
#                 #DTECsList.append(numpy.add(DTECsIter2Pol1, DTECsIter2Pol2) / 2)
#                 DTECsList.append(numpy.nan_to_num(fillGaps1D(numpy.add(DTECsIter2Pol1, DTECsIter2Pol2) / 2)))
#
#
#                 if (plot):
#                     print ("Starting DeltaTEC and phase bandpass visualisation for antenna ID " + str(i) + "...")
#
#                     # Create images of the time derivative of the gain phases (for both polarisations).
#                     polarisationDataGrids = [derivTimePol1, derivTimePol2]
#                     for polarisationName, polarisationDataGrid in zip(polarisationNames, polarisationDataGrids):
#                         figure = pyplot.figure(figsize = (12, 6))
#                         image = pyplot.imshow(polarisationDataGrid, aspect = "auto", cmap = cm.coolwarm, vmin = -2, vmax = 2)
#                         pyplot.xlabel("time (s)")
#                         pyplot.ylabel("frequency (MHz)")
#                         pyplot.xticks(  numpy.linspace(0, numberOfTimeStamps - 1 - 1, num = 5, endpoint = True),
#                                         numpy.linspace(timeStart + timeBinInterval / 2, timeStart + timeRange - timeBinInterval / 2, num = 5, endpoint = True).astype(int))
#                         pyplot.yticks(  numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
#                                         numpy.linspace(frequencyStart, frequencyStart + frequencyRange, num = 5, endpoint = True).astype(int))
#                         pyplot.subplots_adjust(left = .08, right = .9)
#                         pyplot.title("time derivative of antenna-based gain phases of uncalibrated calibrator visibilities\ndata set: "
#                                      + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                      + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                         colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "3%", pad = .05)
#                         colorBar = pyplot.colorbar(image, cax = colorBarAxis)
#                         colorBar.ax.set_ylabel(r"$\frac{\partial\phi}{\partial t}$ ($\degree \cdot \mathrm{s}^{-1}$)", size = 13)
#                         pyplot.subplots_adjust(left = .06, right = .92, bottom = .08, top = .91)
#                         pyplot.savefig(pathPlotDirectory + "derivTime_Ant" + str(i) + "_Pol_" + polarisationName + ".pdf")
#                         pyplot.close()
#
#                     # Create plot of the time derivative of DeltaTEC (for both polarisations, iteration 1).
#                     figure = pyplot.figure(figsize = (12, 6))
#                     timesIntermediate = (gridTimes[0, : -1] + gridTimes[0, 1 : ]) / 2
#                     pyplot.scatter(timesIntermediate, derivTimeDTECsIter1Pol1, color = "navy", label = "data polarisation 1")
#                     pyplot.scatter(timesIntermediate, derivTimeDTECsIter1Pol2, color = "orangered", label = "data polarisation 2")
#                     pyplot.axhline(y = numpy.mean(derivTimeDTECsIter1Pol1), ls = "--", alpha = 1, color = "cornflowerblue", label = "mean polarisation 1")
#                     pyplot.axhline(y = numpy.mean(derivTimeDTECsIter1Pol2), ls = "--", alpha = 1, color = "lightsalmon", label = "mean polarisation 2")
#                     pyplot.xlabel("time (s)")
#                     pyplot.ylabel(r"$\frac{\mathrm{d}\Delta\mathrm{TEC}}{\mathrm{d}t}\ (\mathrm{TECU}\cdot\mathrm{s}^{-1})$")
#                     pyplot.ticklabel_format(style = "sci", axis = "y", scilimits = (0, 0))
#                     pyplot.xlim(timesIntermediate[0] - plotTimeLimit, timesIntermediate[-1] + plotTimeLimit)
#                     pyplot.ylim(-plotDerivTimeDTECLimit, plotDerivTimeDTECLimit)
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.title("time derivative of differential total electron content (iteration 1)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .08, right = .98, bottom = .08, top = .91)
#                     pyplot.savefig(pathPlotDirectory + "derivTimeDTECIter1_Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     # Create images of the residuals of the time derivative of the gain phases: original - fit (for both polarisations).
#                     polarisationDataGrids = [derivTimeResidualsPol1, derivTimeResidualsPol2]
#                     for polarisationName, polarisationDataGrid in zip(polarisationNames, polarisationDataGrids):
#                         figure = pyplot.figure(figsize = (12, 6))
#                         image = pyplot.imshow(polarisationDataGrid, aspect = "auto", cmap = cm.coolwarm, vmin = -2, vmax = 2)
#                         pyplot.xlabel("time (s)")
#                         pyplot.ylabel("frequency (MHz)")
#                         pyplot.xticks(  numpy.linspace(0, numberOfTimeStamps - 1 - 1, num = 5, endpoint = True),
#                                         numpy.linspace(timeStart + timeBinInterval / 2, timeStart + timeRange - timeBinInterval / 2, num = 5, endpoint = True).astype(int))
#                         pyplot.yticks(  numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
#                                         numpy.linspace(frequencyStart, frequencyStart + frequencyRange, num = 5, endpoint = True).astype(int))
#                         pyplot.subplots_adjust(left = .08, right = .9)
#                         pyplot.title("residuals of time derivative of antenna-based gain phases of uncalibrated calibrator visibilities\ndata set: "
#                                      + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                      + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                         colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "3%", pad = .05)
#                         colorBar = pyplot.colorbar(image, cax = colorBarAxis)
#                         colorBar.ax.set_ylabel(r"$\frac{\partial\phi}{\partial t}\ (\degree\cdot\mathrm{s}^{-1})$", size = 13)
#                         pyplot.subplots_adjust(left = .06, right = .92, bottom = .08, top = .91)
#                         pyplot.savefig(pathPlotDirectory + "derivTimeResiduals_Ant" + str(i) + "_Pol_" + polarisationName + ".pdf")
#                         pyplot.close()
#
#
#                     # Create plot of DeltaTEC (for both polarisations, iteration 1).
#                     figure = pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridTimes[0, : ], DTECsIter1Pol1, color = "navy", label = "data polarisation 1")
#                     pyplot.scatter(gridTimes[0, : ], DTECsIter1Pol2, color = "orangered", label = "data polarisation 2")
#                     pyplot.axhline(y = numpy.mean(DTECsIter1Pol1), ls = "--", alpha = 1, color = "cornflowerblue", label = "mean polarisation 1")
#                     pyplot.axhline(y = numpy.mean(DTECsIter1Pol2), ls = "--", alpha = 1, color = "lightsalmon", label = "mean polarisation 2")
#                     pyplot.xlabel("time (s)")
#                     pyplot.ylabel(r"$\Delta\mathrm{TEC}\ (\mathrm{TECU})$")
#                     pyplot.xlim(timeStart - plotTimeLimit, timeStart + timeRange + plotTimeLimit)
#                     pyplot.ylim(-plotDTECLimit, plotDTECLimit)
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.title("differential total electron content (iteration 1)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .08, right = .98, bottom = 0.08, top = .91)
#                     pyplot.savefig(pathPlotDirectory + "DTECIter1_Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     # Create plot of DeltaTEC (for both polarisations, iteration 2).
#                     figure = pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridTimes[0, : ], DTECsIter2Pol1, color = "navy", label = "data polarisation 1")
#                     pyplot.scatter(gridTimes[0, : ], DTECsIter2Pol2, color = "orangered", label = "data polarisation 2")
#                     pyplot.axhline(y = numpy.nanmean(DTECsIter2Pol1), ls = "--", alpha = 1, color = "cornflowerblue", label = "mean polarisation 1")
#                     pyplot.axhline(y = numpy.nanmean(DTECsIter2Pol2), ls = "--", alpha = 1, color = "lightsalmon", label = "mean polarisation 2")
#                     pyplot.xlabel("time (s)")
#                     pyplot.ylabel(r"$\Delta\mathrm{TEC}\ (\mathrm{TECU})$")
#                     pyplot.xlim(timeStart - plotTimeLimit, timeStart + timeRange + plotTimeLimit)
#                     pyplot.ylim(-plotDTECLimit, plotDTECLimit)
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.title("differential total electron content (iteration 2)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = .08, right = .98, bottom = .08, top = .91)
#                     pyplot.savefig(pathPlotDirectory + "DTECIter2_Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#
#                     # Create plot of phase bandpass (for both polarisations, iteration 1).
#                     pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassPhaseIter1Pol1, c = "navy", lw = 0, label = "polarisation 1", s = 16)
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassPhaseIter1Pol2, c = "orangered", lw = 0, label = "polarisation 2", s = 16)
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.xlabel("frequency channel centre (MHz)")
#                     pyplot.ylabel("gain phase ($\degree$)")
#                     pyplot.xlim(frequencyStart - plotFrequencyLimit, frequencyStart + frequencyRange + plotFrequencyLimit)
#                     pyplot.ylim(-180 - plotPhaseLimit, 180 + plotPhaseLimit)
#                     pyplot.yticks(numpy.linspace(-180, 180, num = 9, endpoint = True))
#                     pyplot.title("phase bandpass (iteration 1)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "bandpassPhaseIter1_" + "Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     # Create plot of phase bandpass (for both polarisations, iteration 2).
#                     pyplot.figure(figsize = (12, 6))
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassPhaseIter2Pol1, c = "navy", lw = 0, label = "polarisation 1", s = 16)
#                     pyplot.scatter(gridFrequencies[ : , 0], bandpassPhaseIter2Pol2, c = "orangered", lw = 0, label = "polarisation 2", s = 16)
#                     pyplot.grid(linestyle = "--")
#                     pyplot.legend()
#                     pyplot.xlabel("frequency channel centre (MHz)")
#                     pyplot.ylabel("gain phase ($\degree$)")
#                     pyplot.xlim(frequencyStart - plotFrequencyLimit, frequencyStart + frequencyRange + plotFrequencyLimit)
#                     pyplot.ylim(-180 - plotPhaseLimit, 180 + plotPhaseLimit)
#                     pyplot.yticks(numpy.linspace(-180, 180, num = 9, endpoint = True))
#                     pyplot.title("phase bandpass (iteration 2)\ndata set: "
#                                  + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: $\mathbf{" + str(i) + "}$ | calibrator: "
#                                  + nameField) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                     pyplot.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.08, top = 0.91)
#                     pyplot.savefig(pathPlotDirectory + "bandpassPhaseIter2_" + "Ant" + str(i) + ".pdf")
#                     pyplot.close()
#
#                     if (plot3D):
#                         # Create 3D plots of bandpass-subtracted phases only (iteration 1).
#                         polarisationDataGrids = [gridIonosOnlyPol1, gridIonosOnlyPol2]
#                         polarisationFlagGrids = [flagsPol1[i], flagsPol2[i]]
#                         for polarisationName, polarisationDataGrid, polarisationFlagGrid in zip(polarisationNames, polarisationDataGrids, polarisationFlagGrids):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             gridDataMasked = polarisationDataGrid * (1 - polarisationFlagGrid)
#                             axes3D.plot_surface(gridFrequencies, gridTimes, gridDataMasked, cmap = cm.YlGnBu, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("phase ($\degree$)")
#                             axes3D.set_zlim(-60, 60)
#                             axes3D.set_title("bandpass-subtracted antenna-based gain phases (iteration 1)\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "ionosOnly3DIter1_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#
#                         # Create 3D plots of antenna-based gain phase residuals (iteration 1).
#                         polarisationDataGrids = [gridResidualsIter1Pol1, gridResidualsIter1Pol2]
#                         polarisationSTDs = [STDIter1Pol1, STDIter1Pol2]
#                         polarisationFlagGrids = [flagsPol1[i], flagsPol2[i]]
#                         for polarisationName, polarisationDataGrid, polarisationFlagGrid, polarisationSTD in zip(polarisationNames, polarisationDataGrids, polarisationFlagGrids, polarisationSTDs):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             axes3D.plot_surface(gridFrequencies, gridTimes, polarisationDataGrid.data * (1 - polarisationFlagGrid), cmap = cm.YlOrRd, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("phase ($\degree$)")
#                             axes3D.set_zlim(-9, 9)
#                             axes3D.set_zticks(numpy.array([-9, -6, -3, 0, 3, 6, 9]).astype(int))
#                             axes3D.set_title("antenna-based gain phase residuals (iteration 1) | $\sigma = " + str(numpy.round(polarisationSTD, 2)) + "\degree$\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "residualsIter1_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#
#                         # Create 3D plots of antenna-based gain phase residuals (iteration 2).
#                         polarisationDataGrids = [gridResidualsIter2Pol1, gridResidualsIter2Pol2]
#                         polarisationSTDs = [STDIter2Pol1, STDIter2Pol2]
#                         polarisationFlagGrids = [flagsPol1[i], flagsPol2[i]]
#                         for polarisationName, polarisationDataGrid, polarisationFlagGrid, polarisationSTD in zip(polarisationNames, polarisationDataGrids, polarisationFlagGrids, polarisationSTDs):
#                             figure = pyplot.figure(figsize = (12, 6))
#                             axes3D = figure.add_subplot(111, projection = "3d")
#                             axes3D.plot_surface(gridFrequencies, gridTimes, polarisationDataGrid.data * (1 - polarisationFlagGrid), cmap = cm.YlOrRd, alpha = .8, rstride = 1, cstride = 1)
#                             axes3D.view_init(elev = angleElevation, azim = angleAzimuthal)
#                             axes3D.set_xlabel("frequency (MHz)")
#                             axes3D.set_ylabel("time (s)")
#                             axes3D.set_zlabel("phase ($\degree$)")
#                             axes3D.set_zlim(-9, 9)
#                             axes3D.set_zticks(numpy.array([-9, -6, -3, 0, 3, 6, 9]).astype(int))
#                             axes3D.set_title("antenna-based gain phase residuals (iteration 2) | $\sigma = " + str(numpy.round(polarisationSTD, 2)) + "\degree$\ndata set: "
#                                          + nameMS[ : -3] + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
#                                          + nameField + " | polarisation: " + polarisationName) # 'nameMS[ : -3]' ensures that we remove '.MS' from the name.
#                             pyplot.subplots_adjust(left = .08, right = .9, bottom = .05, top = .95)
#                             pyplot.savefig(pathPlotDirectory + "residualsIter2_Ant" + str(i) + "_Pol" + polarisationName + ".png")
#                             pyplot.close()
#             else:
#                 print ("Skipping DeltaTEC and phase bandpass calculation for antenna ID " + str(i) + ": all data are flagged.")
#
#         # Save the DTEC values found to a NumPy array, alongside the IDs of the antennae.
#         numpy.save(pathOutputBPsPhase, [bandpassesPhasePol1, bandpassesPhasePol2])
#         numpy.save(pathOutputDTECs, [DTECsList, antennaeWorking])
#
#
# if (__name__ = "__main__"):
#
#     # If the program is run from the command line, parse arguments.
#     parser                      = argparse.ArgumentParser(description = "Pipeline step 3: Generation of bandpasses.")
#
#     # Temporary!
#     pathH5Parm                  = "./scanID1.H5"
#
#     dedicated_uGMRT_bandpass(pathH5Parm)