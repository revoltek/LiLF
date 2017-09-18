#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse, logging, sys

from losoto import h5parm
from matplotlib import cm
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy

import lib_ms, lib_util



def filterMedian1D(array, kernelSize):
    """
    Median filter a 1D array (except for the edges).
    'kernelSize' must be an odd integer, bigger than or equal to 3.
    """
    kernelRadius               = (kernelSize - 1) / 2

    arrayNew                   = numpy.zeros_like(array)
    arrayNew[ : kernelRadius]  = array[ : kernelRadius]
    arrayNew[-kernelRadius : ] = array[-kernelRadius : ]

    for index in range(kernelRadius, len(array) - kernelRadius):
        arrayNew[index] = numpy.nanmedian(array[index - kernelRadius : index + kernelRadius + 1])

    return arrayNew



def fillGaps1D(array):
    """
    Replace 'numpy.nan' values wherever they lie between parts of 'array' that have measured values
    via linear interpolation.
    Occurences of 'numpy.nan' at the beginning or end of the array are not replaced, as it is unclear
    how this should be done via interpolation.
    """
    indexDatumLast = None # index of the last encountered valid value (not 'numpy.nan')
    arrayNew       = numpy.copy(array)

    for index in range(len(array)):
        if (not numpy.isnan(array[index])):

            if (not (indexDatumLast == None) and indexDatumLast < index - 1):
                interpolation                        = numpy.linspace(array[indexDatumLast], array[index], num = index - indexDatumLast + 1, endpoint = True)[1 : -1]
                arrayNew[indexDatumLast + 1 : index] = interpolation
            indexDatumLast = index

    return arrayNew



def sigmaClip1D(array, numberOfSigmas = 3.0, iterationNumberMax = 100, verbose = True):
    '''
    This function takes a 1D array and performs sigma clipping on it.
    The mean, median and standard deviation are returned.

    'array' may contain NaNs. 'arraySliced' does not.
    '''

    isNaN       = numpy.isnan(array)
    arraySliced = array[numpy.logical_not(isNaN)]

    if (verbose):
        print ("sigma clipping has started...")
        print ("# array elements (initially): " + str(len(array)) + " | # NaNs: " + str(numpy.sum(isNaN)) + " | # array elements (remaining): " + str(len(arraySliced)))


    for iterationCurrent in range(1, iterationNumberMax + 1):
        if (verbose):
            print ("- iteration " + str(iterationCurrent) + " of " + str(iterationNumberMax) + " -")

        median           = numpy.median(arraySliced)
        sigma            = numpy.std(arraySliced)
        outlierBooleans  = numpy.greater(numpy.absolute(arraySliced - median), numberOfSigmas * sigma)
        numberOfOutliers = numpy.sum(outlierBooleans)


        if (numberOfOutliers == 0): # If convergence is reached, the loop can end.
            if (verbose):
                print("# new outliers: 0 | # array elements (remaining): " + str(len(arraySliced)))
            break
        else:
            arraySliced = arraySliced[numpy.logical_not(outlierBooleans)]
            if (verbose):
                print ("# new outliers: " + str(numberOfOutliers) + " | # array elements (remaining): " + str(len(arraySliced)))

    if (verbose):
        print ("sigma clipping has ended...")

    # Computes the 1D mean, median and standard deviation arrays.
    mean        = numpy.mean(arraySliced)
    median      = numpy.median(arraySliced)
    sigma       = numpy.std(arraySliced)

    return (mean, median, sigma)



def wrapPhasesZeroCentred(phases, unitDegree = True):
        """
        This method assumes that 'phases' is expressed on a -180 to 180 degree scale (or -pi to pi if 'unitDegree' is False),
        with some of the phases lying out of bounds - which is the motivation to call this function.
        It correctly phase wraps these phases, and returns them expressed on the same scale.
        """
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



def computeDerivative1D(array, stepSize, degree = 1, intermediateSampling = False):
    """
    This function computes a numerical approximation to the first or second derivative, determined via 'degree'.

    'intermediateSampling = False': if coordinates are t_0, t_1, ..., t_{N-2}, t_{N-1},
    then the derivatives are given at matching points t_1, t_2, ..., t_{N-3}, t_{N-2}.

    'intermediateSampling = True': if coordinats are t_0, t_1, ..., t_{N_2}, t_{N-1},
    then the derivatives are given at intermediate points (t_0 + t_1) / 2, ..., (t_{N-2} + t_{N-1}) / 2.
    """
    if (degree == 0):
        arrayResult = array
    elif (degree == 1):
        if (intermediateSampling):
            arrayResult = (array[1 : ] - array[ : -1]) / stepSize
        else:
            arrayResult = (array[2 : ] - array[ : -2]) / (2 * stepSize)
    elif (degree == 2):
        if (intermediateSampling):
            print ("Error: intermediate sampling is not yet featured for the second derivative.")
            sys.exit()
        else:
            arrayResult = (array[4 : ] + array[ : -4] - 2 * array[2 : -2]) / (4 * stepSize ** 2)
    else:
        print ("Error: 'computeDerivative1D' does not support calculating derivatives of degree " + str(degree) + ".")
        sys.exit()

    return arrayResult



def computeDerivative2D(grid, stepSize, axis = 0, degree = 1, intermediateSampling = False):
    """
    This function computes a numerical approximation to the first or second derivative.
    Derivatives can either be taken along the rows ('axis = 0'), or along the columns ('axis = 1').

    If the coordinates of the columns of 'grid' are t_0, t_1, t_2, ..., t_{N-1}, then
    the derivatives are given at t_1, t_2, ..., t_{N-2}.

    Sometimes, e.g. for numerical integration, it is desirable that derivatives are given at
    the intermediate points: (t_0 + t_1) / 2, (t_1 + t_2) / 2, ..., (t_{N-2} + t_{N-1}) / 2.
    In such cases, set 'intermediateSampling' to True. Beware that output array dimensions change!
    """
    if (axis == 0):
        if (degree == 0):
            gridResult = grid
        elif (degree == 1):
            if (intermediateSampling):
                gridResult = (grid[1 : , : ] - grid[ : -1, : ]) / stepSize
            else:
                gridResult = (grid[2 : , : ] - grid[ : -2, : ]) / (2 * stepSize)
        elif (degree == 2):
            if (intermediateSampling):
                print ("Error: intermediate sampling is not yet featured for the second derivative.")
                sys.exit()
            else:
                gridResult = (grid[4 : , : ] + grid[ : -4, : ] - 2 * grid[2 : -2, : ]) / (4 * stepSize ** 2)
        else:
            print ("Error: 'computeDerivative2D' does not support calculating derivatives of degree " + str(degree) + ".")
            sys.exit()
    else:
        if (degree == 0):
            gridResult = grid
        elif (degree == 1):
            if (intermediateSampling):
                gridResult = (grid[ : , 1 : ] - grid[ : , : -1]) / stepSize
            else:
                gridResult = (grid[ : , 2 : ] - grid[ : , : -2]) / (2 * stepSize)
        elif (degree == 2):
            if (intermediateSampling):
                print ("Error: intermediate sampling is not yet featured for the second derivative.")
                sys.exit()
            else:
                gridResult = (grid[ : , 4 : ] + grid[ : , : -4] - 2 * grid[ : , 2 : -2]) / (4 * stepSize ** 2)
        else:
            print ("Error: 'computeDerivative2D' does not support calculating derivatives of degree " + str(degree) + ".")
            sys.exit()

    return gridResult



def calculateBandpassAmplitude(gainAmplitudes, flags, numberOfSubiterations = 3, flaggingThresholdFactor = 3, medianFilteringKernelSize = 7):
    """
    Generate an amplitude bandpass from a time-frequency grid of gain amplitudes, for a certain polarisation and antenna.
    The process will output 2 amplitude bandpasses, one for each of 2 iterations.
    Iteration 1 consists out of subiterations, during which outliers are flagged.
    Iteration 2 takes the results from iteration 1 and performs median filtering and interpolation.

    gainAmplitudes: two-dimensional grid of gain amplitudes
    flags:          two-dimensional grid of flags
    """

    _, numberOfTimeStamps = gainAmplitudes.shape

    #
    # Determine amplitude bandpass (iteration 1).
    #
    for subiteration in range(numberOfSubiterations):

        # If in the last subiteration, determine the amplitude bandpass using the mean.
        # During earlier iterations, to avoid e.g. RFI signatures, we use the median.
        if (subiteration == numberOfSubiterations - 1):
            bandpassAmplitudeIter1 = numpy.nanmean(gainAmplitudes, axis = 1)
        else:
            bandpassAmplitudeIter1 = numpy.nanmedian(gainAmplitudes, axis = 1)

            # Tile the bandpass into a grid.
            gridBandpassAmplitude = numpy.tile(bandpassAmplitudeIter1, (numberOfTimeStamps, 1)).T

            # Divide the data by the provisional bandpass. Residuals will be centered around 1.
            gridResiduals         = numpy.divide(gainAmplitudes, gridBandpassAmplitude)

            # Calculate the standard deviation of the residuals.
            STDResiduals          = numpy.nanstd(gridResiduals)

            # Determine outliers and update the flags. This makes sure that the bandpass in the next iterations is better.
            gridIsOutlier         = numpy.greater(numpy.absolute(gridResiduals - 1), flaggingThresholdFactor * STDResiduals)
            flags                 = numpy.logical_or(flags, gridIsOutlier)

            # Set flagged amplitudes to 'numpy.nan'.
            gainAmplitudes        = numpy.where(flags, numpy.nan, gainAmplitudes)

    #
    # Determine amplitude bandpass (iteration 2).
    #

    # Median filter and interpolate the amplitude bandpass.
    bandpassAmplitudeIter2 = fillGaps1D(filterMedian1D(bandpassAmplitudeIter1, kernelSize = medianFilteringKernelSize))

    return bandpassAmplitudeIter1, bandpassAmplitudeIter2



def calculateBandpassPhaseNoIonosphere(gainPhases, useMedian = True):
    """
    Generate a phase bandpass from a time-frequency grid of gain phases, for a certain polarisation and antenna.
    The method is robust to phase wrapping.
    Do this under the assumption that the input gain phases are approximately time-independent; that is,
    ionospheric imprints are already taken out.

    gainPhases: two-dimensional grid of gain phases (zero-centred, in degrees)
    useMedian:  if True, use median to determine bandpass; if False, use mean
    """
    if (useMedian):
        bandpassPhaseChoice1 = numpy.nanmedian(gainPhases, axis = 1)                                                                # in degrees
        bandpassPhaseChoice2 = wrapPhasesZeroCentred(numpy.nanmedian(numpy.mod(gainPhases + 180 + 180, 360), axis = 1) - 180 - 180) # in degrees
    else:
        bandpassPhaseChoice1 = numpy.nanmean(gainPhases, axis = 1)                                                                  # in degrees
        bandpassPhaseChoice2 = wrapPhasesZeroCentred(numpy.nanmean(numpy.mod(gainPhases + 180 + 180, 360), axis = 1) - 180 - 180)   # in degrees

    # Calculate, for each frequency channel, the variability in the gain phases. Variability can arise from phase wrapping.
    STDChoice1    = numpy.nanstd(gainPhases, axis = 1)                                                                              # in degrees
    STDChoice2    = numpy.nanstd(numpy.mod(gainPhases + 180 + 180, 360), axis = 1)                                                  # in degrees

    # Construct the bandpass by choosing, for each frequency channel, the bandpass option that is associated with least variability.
    bandpassPhase = numpy.where(numpy.less(STDChoice1, STDChoice2), bandpassPhaseChoice1, bandpassPhaseChoice2)                     # in degrees

    return bandpassPhase



def calculateBandpassPhase(gainPhases, flags, frequencies, times, numberOfSubiterations = 3, flaggingThresholdFactor = 3, medianFilteringKernelSize = 7):
    """
    Generate phase bandpasses from a time-frequency grid of gain phases, for a certain polarisation and antenna.
    The process will output 2 phase bandpasses, one for each of 2 iterations.
    Iteration 1 consists out of subiterations, during which outliers are flagged.
    Iteration 2 takes the results from iteration 1 and performs median filtering and interpolation.

    gainPhases:  two-dimensional array of gain phases
    flags:       two-dimensional array of flags
    frequencies: one-dimensional array of frequencies
    times:       one-dimensional array of times
    """

    numberOfChannels, numberOfTimeStamps = gainPhases.shape

    timeStampDuration             = times[1] - times[0]             # in s (we assume that the time stamp lengths are uniform)
    frequencyChannelWidth         = frequencies[1] - frequencies[0] # in MHz

    # Create a 2D and a 1D array of inverse frequencies over the frequency band.
    gridIndexRow, gridIndexColumn = numpy.mgrid[0 : numberOfChannels, 0 : numberOfTimeStamps]
    gridFrequencies               = numpy.divide(gridIndexRow, numberOfChannels - 1) * (frequencies[-1] - frequencies[0]) + frequencies[0] # in MHz
    gridFrequenciesInverse        = numpy.divide(1, gridFrequencies)                                                                       # in MHz^-1
    frequenciesInverse            = gridFrequenciesInverse[ : , 0]                                                                         # in MHz^-1



    # Calculate the first derivative of gain phase to time in a way robust to phase wrapping (for both polarisations).
    # We do so by calculating the derivative for each time-frequency bin 2 times: one with the ordinary data, and once after shifting
    # - all the phases by 180 degrees to place them in the [0, 360) domain, and
    # - another 180 degrees to effect a maximum shift within that domain.

    derivTimeChoice1        = computeDerivative2D(gainPhases,                             timeStampDuration, axis = 1, degree = 1, intermediateSampling = True) # in degrees per second
    derivTimeChoice2        = computeDerivative2D(numpy.mod(gainPhases + 180 + 180, 360), timeStampDuration, axis = 1, degree = 1, intermediateSampling = True) # in degrees per second
    derivTime               = numpy.where(numpy.less(numpy.absolute(derivTimeChoice1), numpy.absolute(derivTimeChoice2)), derivTimeChoice1, derivTimeChoice2)   # in degrees per MHz^2


    # Determine the DTEC time derivative for both polarisations (iteration 1).
    derivTimeDTECsIter1Grid = derivTime * gridFrequencies[ : , : -1] / (1210 * 400)
    derivTimeDTECsIter1     = numpy.nanmean(derivTimeDTECsIter1Grid, axis = 0)

    '''
    # Create a grid of residuals. CURRENTLY NOT USED IN THE ALGORITHM!
    derivTimeFitPol1            = numpy.multiply(gridFrequenciesInverse[ : , : -1], numpy.tile(derivTimeDTECsIter1Pol1, (numberOfChannels, 1))) * 1210 * 400 # in degrees per second
    derivTimeResidualsPol1      = derivTimePol1 - derivTimeFitPol1

    # To do: flag or mask the worst outliers, and fit 'derivTimeDTECsIter1PolX' again with the new flags.
    '''

    # Whenever there are NaNs in 'derivTimeDTECsIter1', then 'DTECsIter1' is not well-determined.
    # We assume continuity of the time derivative of DTEC, and thus interpolate.
    derivTimeDTECsIter1     = fillGaps1D(derivTimeDTECsIter1)

    # Find DTECs by integrating along time (iteration 1).
    DTECsIter1              = [0]
    for j in range(numberOfTimeStamps - 1):
        DTECsIter1.append(timeStampDuration * numpy.nansum(derivTimeDTECsIter1[ : j + 1]))

    # NumPy-ify the 1D array.
    DTECsIter1              = numpy.array(DTECsIter1)

    # Flag DTEC values appropriately, or do median filtering?
    '''
    '''

    # Subtract the mean or median from the DTEC values, so that they lie around 0.
    DTECsIter1             -= numpy.mean(DTECsIter1)

    # Calculate the fitted plasma opacity phase effect (iteration 1).
    gridPlasmaOpIter1       = wrapPhasesZeroCentred(numpy.multiply(gridFrequenciesInverse, numpy.tile(DTECsIter1, (numberOfChannels, 1))) * 1210 * 400) # in degrees


    # We now determine the bandpass by subtracting the plasma opacity phase effect from the original data.
    # As we flag outliers in the process, we repeat this procedure several times.
    # During the last subiteration, the bandpass is calculated by taking the mean along time, whereas in
    # earlier iterations we calculate the bandpass using the median.
    for subiteration in range(numberOfSubiterations):

        # Subtract the time-variable ionospheric effect from the gain phases.
        gridStaticIter1    = wrapPhasesZeroCentred(gainPhases - gridPlasmaOpIter1)

        # Determine bandpass. Use median for all but the last subiteration.
        bandpassPhaseIter1 = calculateBandpassPhaseNoIonosphere(gridStaticIter1, useMedian = (subiteration < numberOfSubiterations - 1))

        # Remove the phase bandpass from the static pattern to keep residuals, which we will use for flagging.
        gridResidualsIter1 = wrapPhasesZeroCentred(gridStaticIter1 - numpy.tile(bandpassPhaseIter1, (numberOfTimeStamps, 1)).T)

        # Calculate the standard deviation of the residuals.
        STDIter1           = numpy.nanstd(gridResidualsIter1)

        # Determine outliers and update the flags. This makes sure that the bandpass in the next iterations is better.
        gridIsOutlier      = numpy.greater(numpy.absolute(gridResidualsIter1), flaggingThresholdFactor * STDIter1)
        flags              = numpy.logical_or(flags, gridIsOutlier)

        # Set flagged phases to 'numpy.nan'.
        gainPhases         = numpy.where(flags, numpy.nan, gainPhases)

        # Temporary debug output!
        print ("# NaNs bandpassPhaseIter1:", numpy.sum(numpy.isnan(bandpassPhaseIter1)))
        print ("# NaNs gainPhases:", numpy.sum(numpy.isnan(gainPhases)))


    # Remove the phase bandpass from the original data.
    gridIonosOnly               = wrapPhasesZeroCentred(gainPhases - numpy.tile(bandpassPhaseIter1, (numberOfTimeStamps, 1)).T)

    # Calculate DTECs directly by fitting to the ionospheric data only (iteration 2).
    DTECsIter2                  = numpy.nanmean(gridIonosOnly * gridFrequencies / (1210 * 400), axis = 0) # in TECUs

    # Calculate the fitted plasma opacity phase effect (iteration 2).
    gridPlasmaOpIter2           = wrapPhasesZeroCentred(numpy.multiply(gridFrequenciesInverse, numpy.tile(DTECsIter2, (numberOfChannels, 1))) * 1210 * 400) # in degrees



    # Subtract the time-variable ionospheric effect from the gain phases.
    gridStaticIter2             = wrapPhasesZeroCentred(gainPhases - gridPlasmaOpIter2)

    # Determine bandpass.
    bandpassPhaseIter2          = calculateBandpassPhaseNoIonosphere(gridStaticIter2, useMedian = False)

    '''
    # Remove the phase bandpass from the static pattern to keep residuals, which we will use for flagging.
    gridResidualsIter2          = wrapPhasesZeroCentred(gridStaticIter2 - numpy.tile(bandpassPhaseIter1, (numberOfTimeStamps, 1)).T)

    # Calculate the standard deviation of the residuals.
    STDIter2                    = numpy.nanstd(gridResidualsIter2)
    '''

    # Calculate the derivative of the phase bandpass to frequency.
    bandpassPhaseDeriv          = computeDerivative1D(bandpassPhaseIter2, stepSize = frequencyChannelWidth, intermediateSampling = True)

    # Determine the mean, median and standard deviation of the derivative after sigma clipping.
    mean, median, sigma         = sigmaClip1D(bandpassPhaseDeriv, verbose = False) # in degrees per MHz

    # Determine whether the mean or median should be used for generating the bandpass fit.
    if (sigma < 10): # in degrees per MHz
        slope = mean
        print("Phase bandpass fit | slope type: mean")
    else:
        slope = median
        print("Phase bandpass fit | slope type: median")

    # Calculate the bandpass fit.
    BPPhaseFit                  = []

    for indexRef in range(numberOfChannels):
        if (not numpy.isnan(bandpassPhaseIter2[indexRef])):
            break

    for j in range(numberOfChannels):
        BPPhaseFit.append(bandpassPhaseIter2[indexRef] + frequencyChannelWidth * (j - indexRef) * slope)

    # NumPy-ify the 1D array, and ensure phase wrapping.
    BPPhaseFit                  = wrapPhasesZeroCentred(numpy.array(BPPhaseFit))

    # Calculate the residuals that remain after the fit has been subtracted from the data.
    BPPhaseRes1                 = wrapPhasesZeroCentred(bandpassPhaseIter2 - BPPhaseFit)

    # We hope we can now apply median filtering and interpolation, as we removed the phase ramps.
    BPPhaseRes2                 = fillGaps1D(filterMedian1D(BPPhaseRes1, kernelSize = medianFilteringKernelSize))

    # Now add back in the phase ramps.
    bandpassPhaseIter3          = wrapPhasesZeroCentred(BPPhaseRes2 + BPPhaseFit)

    '''
    BPPhaseFit = [bandpassPhaseIter1[0]]
    BPPhaseFitPol2 = [bandpassPhaseIter1Pol2[0]]
    for j in range(numberOfChannels - 1):
        BPPhaseFit.append(BPPhaseFit[0] + frequencyBinInterval * (j + 1) * slope)
        BPPhaseFitPol2.append(BPPhaseFitPol2[0] + frequencyBinInterval * (j + 1) * slopePol2)

    pyplot.scatter(range(len(bandpassPhaseDeriv)), bandpassPhaseDeriv,s = 4)
    pyplot.scatter(range(len(bandpassPhaseDerivPol2)), bandpassPhaseDerivPol2, s= 4)
    pyplot.axhline(y = median)
    pyplot.axhline(y = medianPol2)
    pyplot.axhline(y = mean, ls = "--")
    pyplot.axhline(y = meanPol2, ls = "--")
    pyplot.axhline(y = median - sigma, ls = ":", c = 'b')
    pyplot.axhline(y = median + sigma, ls = ":", c = 'b')
    pyplot.show()

    pyplot.scatter(range(len(BPPhaseRes1)), BPPhaseRes1, s= 4)
    pyplot.scatter(range(len(BPPhaseRes1Pol2)), BPPhaseRes1Pol2, s=4)
    pyplot.show()
    '''

    return bandpassPhaseIter2, bandpassPhaseIter3, DTECsIter2



def plotAmplitudes2D(amplitudes, times, frequencies, antennaeWorking, pathDirectoryPlots,
                     namePolarisation = "?", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate time-frequency plots of antenna-based gain amplitudes, of one specific polarisation.
    """

    # Load dimensions.
    numberOfAntennae, numberOfChannels, numberOfTimeStamps = amplitudes.shape

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            logging.info("Starting gain amplitudes visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + "...")

            # Create 2D antenna-based gain amplitudes plot.
            figure = pyplot.figure(figsize = (12, 6))
            image  = pyplot.imshow(amplitudes[i], aspect = "auto", interpolation = "nearest", cmap = cm.viridis, vmin = 0, vmax = 1)

            pyplot.xlabel("time (s)")
            pyplot.ylabel("frequency (MHz)")
            pyplot.xticks(numpy.linspace(0, numberOfTimeStamps - 1, num = 5, endpoint = True),
                          numpy.linspace(times[0], times[-1], num = 5, endpoint = True).astype(int))
            pyplot.yticks(numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
                          numpy.linspace(frequencies[0], frequencies[-1], num = 5, endpoint = True).astype(int))
            pyplot.title("$\mathbf{antenna-based\ gain\ amplitudes\ of\ uncalibrated\ calibrator\ visibilities}$\ndata set: "
                         + nameDataSet + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
                         + nameField + " | polarisation: " + namePolarisation, fontsize = 9)

            colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
            colorBar     = pyplot.colorbar(image, cax = colorBarAxis, ticks = [0, 0.2, 0.4, 0.6, 0.8, 1])
            colorBar.ax.set_ylabel("antenna-based gain amplitude (1)")

            pyplot.subplots_adjust(left = .05, right = .93, bottom = 0.08, top = 0.91)
            pyplot.savefig(pathDirectoryPlots + "/amplitudes2D_ant" + str(i) + "_pol" + namePolarisation + ".pdf")
            pyplot.close()
        else:
            logging.info("Skipping gain amplitudes visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + ": all data are flagged.")



def plotPhases2D(phases, times, frequencies, antennaeWorking, pathDirectoryPlots,
                 namePolarisation = "?", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate time-frequency plots of antenna-based gain phases, of one specific polarisation.
    """

    # Load dimensions.
    numberOfAntennae, numberOfChannels, numberOfTimeStamps = phases.shape

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            logging.info("Starting gain phases visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + "...")

            # Create 2D antenna-based gain phases plot.
            figure = pyplot.figure(figsize = (12, 6))
            image  = pyplot.imshow(phases[i], aspect = "auto", interpolation = "none", cmap = cm.hsv, vmin = -180, vmax = 180)

            pyplot.xlabel("time (s)")
            pyplot.ylabel("frequency (MHz)")
            pyplot.xticks(numpy.linspace(0, numberOfTimeStamps - 1, num = 5, endpoint = True),
                          numpy.linspace(times[0], times[-1], num = 5, endpoint = True).astype(int))
            pyplot.yticks(numpy.linspace(0, numberOfChannels - 1, num = 5, endpoint = True),
                          numpy.linspace(frequencies[0], frequencies[-1], num = 5, endpoint = True).astype(int))
            pyplot.title("$\mathbf{antenna-based\ gain\ phases\ of\ uncalibrated\ calibrator\ visibilities}$\ndata set: "
                         + nameDataSet + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
                         + nameField + " | polarisation: " + namePolarisation, fontsize = 9)

            colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
            colorBar     = pyplot.colorbar(image, cax = colorBarAxis, ticks = [-180, -120, -60, 0, 60, 120, 180])
            colorBar.ax.set_ylabel("antenna-based gain phase $(\degree)$")

            pyplot.subplots_adjust(left = .05, right = .93, bottom = 0.08, top = 0.91)
            pyplot.savefig(pathDirectoryPlots + "/phases2D_ant" + str(i) + "_pol" + namePolarisation + ".pdf")
            pyplot.close()
        else:
            logging.info("Skipping gain phases visualisation for antenna ID " + str(i) + " and polarisation " + namePolarisation + ": all data are flagged.")



def plotBandpassesAmplitude(bandpassesAmplitudePol1, bandpassesAmplitudePol2, frequencies, antennaeWorking, pathDirectoryPlots,
                            namesPolarisation = ["?", "?"], nameIteration = "?", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate plots of amplitude bandpasses, for two polarisations.
    """
    numberOfAntennae   = len(antennaeWorking)
    plotFrequencyLimit = 1 # in MHz

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            logging.info("Starting amplitude bandpass visualisation for antenna ID " + str(i) + "...")

            # Create plot of amplitude bandpass (for both polarisations).
            pyplot.figure(figsize = (12, 6))

            pyplot.scatter(frequencies, bandpassesAmplitudePol1[i], c = "navy",      s = 12, lw = 0, label = namesPolarisation[0])
            pyplot.scatter(frequencies, bandpassesAmplitudePol2[i], c = "orangered", s = 12, lw = 0, label = namesPolarisation[1])

            pyplot.grid(linestyle = "--", alpha = 0.5)
            pyplot.legend(title = "polarisation")
            pyplot.xlabel("frequency channel centre (MHz)")
            pyplot.ylabel("antenna-based gain amplitude (1)")
            pyplot.xlim(frequencies[0] - plotFrequencyLimit, frequencies[-1] + plotFrequencyLimit)
            pyplot.ylim(0, 1.05)
            pyplot.title("$\mathbf{amplitude\ bandpass\ (iteration\ " + nameIteration + ")}$\ndata set: "
                         + nameDataSet + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
                         + nameField, fontsize = 9)

            pyplot.subplots_adjust(left = .07, right = .98, bottom = 0.08, top = 0.91)
            pyplot.savefig(pathDirectoryPlots + "/bandpassAmplitude_ant" + str(i) + "_iter" + nameIteration + ".pdf")
            pyplot.close()
        else:
            logging.info("Skipping amplitude bandpass visualisation for antenna ID " + str(i) + ": all data are flagged.")



def plotBandpassesPhase(bandpassesPhasePol1, bandpassesPhasePol2, frequencies, antennaeWorking, pathDirectoryPlots,
                            namesPolarisation = ["?", "?"], nameIteration = "?", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate plots of phase bandpasses, for two polarisations.
    """
    numberOfAntennae   = len(antennaeWorking)
    plotFrequencyLimit = 1 # in MHz

    for i in range(numberOfAntennae):
        if (antennaeWorking[i]):
            logging.info("Starting phase bandpass visualisation for antenna ID " + str(i) + "...")

            # Create plot of phase bandpass (for both polarisations).
            pyplot.figure(figsize = (12, 6))

            pyplot.scatter(frequencies, bandpassesPhasePol1[i], c = "navy",      s = 12, lw = 0, label = namesPolarisation[0])
            pyplot.scatter(frequencies, bandpassesPhasePol2[i], c = "orangered", s = 12, lw = 0, label = namesPolarisation[1])

            pyplot.grid(linestyle = "--", alpha = 0.5)
            pyplot.legend(title = "polarisation")
            pyplot.xlabel("frequency channel centre (MHz)")
            pyplot.ylabel("antenna-based gain phase $(\degree)$")
            pyplot.xlim(frequencies[0] - plotFrequencyLimit, frequencies[-1] + plotFrequencyLimit)
            pyplot.ylim(-180, 180)
            pyplot.title("$\mathbf{phase\ bandpass\ (iteration\ " + nameIteration + ")}$\ndata set: "
                         + nameDataSet + " | telescope: " + nameTelescope + " | antenna ID: " + str(i) + " | calibrator: "
                         + nameField, fontsize = 9)

            pyplot.subplots_adjust(left = .07, right = .98, bottom = 0.08, top = 0.91)
            pyplot.savefig(pathDirectoryPlots + "/bandpassPhase_ant" + str(i) + "_iter" + nameIteration + ".pdf")
            pyplot.close()
        else:
            logging.info("Skipping phase bandpass visualisation for antenna ID " + str(i) + ": all data are flagged.")



def plotBandpassesAmplitude2D(bandpassesAmplitude, pathDirectoryPlots, namePolarisation = "?", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate frequency-antenna plots of amplitude bandpasses, of one specific polarisation.
    """
    logging.info("Starting amplitude bandpass visualisation for polarisation " + namePolarisation + "...")


    pyplot.figure(figsize = (12, 6))

    image = pyplot.imshow(numpy.array(bandpassesAmplitude), aspect = "auto", interpolation = "none")

    pyplot.xlabel("frequency channel index")
    pyplot.ylabel("antenna index")
    pyplot.title("$\mathbf{amplitude\ bandpasses}$\ndata set: " + nameDataSet
                 + " | telescope: " + nameTelescope + " | calibrator: " + nameField
                 + " | polarisation: " + namePolarisation, fontsize = 9)

    colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
    colorBar     = pyplot.colorbar(image, cax = colorBarAxis, ticks = [0, 0.2, 0.4, 0.6, 0.8, 1])
    colorBar.ax.set_ylabel("antenna-based gain amplitude (1)")

    pyplot.subplots_adjust(left = .05, right = .93, bottom = 0.08, top = 0.91)
    pyplot.savefig(pathDirectoryPlots + "/bandpassAmplitudeAll2D_pol" + namePolarisation + ".pdf")
    pyplot.close()



def plotBandpassesPhase2D(bandpassesPhase, pathDirectoryPlots, namePolarisation = "?", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate frequency-antenna plots of phase bandpasses, of one specific polarisation.
    """
    logging.info("Starting phase bandpass visualisation for polarisation " + namePolarisation + "...")


    pyplot.figure(figsize = (12, 6))

    image = pyplot.imshow(numpy.array(bandpassesPhase), aspect = "auto", interpolation = "none", cmap = cm.hsv, vmin = -180, vmax = 180)

    pyplot.xlabel("frequency channel index")
    pyplot.ylabel("antenna index")
    pyplot.title("$\mathbf{phase\ bandpasses}$\ndata set: " + nameDataSet
                 + " | telescope: " + nameTelescope + " | calibrator: " + nameField
                 + " | polarisation: " + namePolarisation, fontsize = 9)

    colorBarAxis = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
    colorBar     = pyplot.colorbar(image, cax = colorBarAxis, ticks = [-180, -120, -60, 0, 60, 120, 180])
    colorBar.ax.set_ylabel("antenna-based gain phase $(\degree)$")

    pyplot.subplots_adjust(left = .05, right = .93, bottom = 0.08, top = 0.91)
    pyplot.savefig(pathDirectoryPlots + "/bandpassPhaseAll2D_pol" + namePolarisation + ".pdf")
    pyplot.close()


def plotFunctionsDTEC2D(functionsDTEC, pathDirectoryPlots, suffixFilename = "", comment = "-", nameField = "?", nameDataSet = "?", nameTelescope = "uGMRT"):
    """
    Generate time-antenna plots of DTECs.
    """
    logging.info("Starting DTEC function visualisation with comment '" + comment + "'...")

    colorScaleMaxDTEC = numpy.round(numpy.nanmax(numpy.absolute(functionsDTEC)), 3)

    pyplot.figure(figsize = (12, 6))

    image             = pyplot.imshow(numpy.array(functionsDTEC), aspect = "auto", interpolation = "none", cmap = cm.coolwarm,
                                      vmin = -colorScaleMaxDTEC, vmax = colorScaleMaxDTEC)

    pyplot.xlabel("time stamp index")
    pyplot.ylabel("antenna index")
    pyplot.title("$\mathbf{\Delta TEC}$\ndata set: " + nameDataSet
                 + " | telescope: " + nameTelescope + " | calibrator: " + nameField
                 + " | comment: " + comment, fontsize = 9)

    colorBarAxis      = make_axes_locatable(pyplot.gca()).append_axes("right", size = "2%", pad = .05)
    colorBar          = pyplot.colorbar(image, cax = colorBarAxis,
                                        ticks = numpy.round(numpy.array([-1, -2 / 3., -1 / 3., 0, 1 / 3., 2 / 3., 1]) * colorScaleMaxDTEC, 3))
    colorBar.ax.set_ylabel("$\Delta\mathrm{TEC}\ (\mathrm{TECU})$")

    pyplot.subplots_adjust(left = .05, right = .93, bottom = 0.08, top = 0.91)
    pyplot.savefig(pathDirectoryPlots + "/functionsDTECAll2D" + suffixFilename + ".pdf")
    pyplot.close()


def dedicated_uGMRT_bandpass(pathDirectoryMS, referenceAntennaID = 0, verbose = False):

    # Initialise logistics.
    nameMS                     = pathDirectoryMS[pathDirectoryMS.rfind('/') + 1 : ]
    pathMS                     = pathDirectoryMS + '/' + nameMS + ".MS"
    pathH5ParmInput            = pathDirectoryMS + "/solutions/gainsRaw.h5"
    pathH5ParmOutput           = pathDirectoryMS + "/solutions/bandpassesTECs.h5"
    pathDirectoryPlots         = pathDirectoryMS + "/plots"


    # Initialise H5Parm file objects.
    objectH5Parm               = h5parm.h5parm(pathH5ParmInput, readonly = True)
    objectSolSet               = objectH5Parm.getSolset("sol000")
    objectSolTabGainAmplitudes = objectSolSet.getSoltab("amplitude000")
    objectSolTabGainPhases     = objectSolSet.getSoltab("phase000")

    # Load antenna-based gains and weights.
    # 'objectSolTabGainAmplitudes.getValues(retAxesVals = False).shape' is e.g. (2, 1, 30, 2048, 75):
    # 2 polarisations, 1 direction, 30 antennae, 2048 frequency channels, 75 time stamps.
    _, axes                    = objectSolTabGainAmplitudes.getValues(retAxesVals = True)

    gainAmplitudes             = objectSolTabGainAmplitudes.getValues(retAxesVals = False, weight = False)[ : , 0, : , : , : ]
    gainPhases                 = objectSolTabGainPhases.getValues(    retAxesVals = False, weight = False)[ : , 0, : , : , : ]

    weightsForAmplitudes       = objectSolTabGainAmplitudes.getValues(retAxesVals = False, weight = True) [ : , 0, : , : , : ]
    weightsForPhases           = objectSolTabGainPhases.getValues(    retAxesVals = False, weight = True) [ : , 0, : , : , : ]

    # Close H5Parm file, freeing up memory.
    objectH5Parm.close()


    # Calculate flags from weights, which are (inverted) generalised flags, in a sense. This program uses only flags.
    flagsForAmplitudes         = numpy.logical_not(weightsForAmplitudes)
    flagsForPhases             = numpy.logical_not(weightsForPhases)


    # Load dimensions.
    numberOfPolarisations, numberOfAntennae, numberOfChannels, numberOfTimeStamps = gainAmplitudes.shape

    # Stop program if deviant data shape found.
    if (numberOfPolarisations != 2):
        logging.error("2 polarisations expected, but " + str(numberOfPolarisations) + " received. Aborting...")
        sys.exit()


    # Initialise axes arrays.
    namesPolarisation          = axes["pol"]
    namesAntenna               = axes["ant"]
    frequencies                = axes["freq"] / 1e6             # in MHz
    times                      = axes["time"] - axes["time"][0] # in s


    # Make gain phases relative to reference antenna, clip and convert from radians to degrees.
    gainPhases                -= numpy.tile(gainPhases[ : , referenceAntennaID : referenceAntennaID + 1, : , : ], (1, numberOfAntennae, 1, 1))
    gainPhases                 = wrapPhasesZeroCentred(gainPhases, unitDegree = False)
    gainPhases                 = numpy.degrees(gainPhases)


    # Split-up gains by polarisation.
    gainAmplitudesPol1         = gainAmplitudes[0]
    gainAmplitudesPol2         = gainAmplitudes[1]
    gainPhasesPol1             = gainPhases[0]
    gainPhasesPol2             = gainPhases[1]

    # Split-up flags by polarisation.
    flagsForAmplitudesPol1   = flagsForAmplitudes[0]
    flagsForAmplitudesPol2   = flagsForAmplitudes[1]
    flagsForPhasesPol1       = flagsForPhases[0]
    flagsForPhasesPol2       = flagsForPhases[1]


    # Flagged data should not be used in calculations.
    # Masking using 'numpy.ma.masked_array(...)' is not always practical - the mask is lost during some NumPy operations.
    # We choose to set flagged amplitudes and phases to 'numpy.nan'. (This leads to undesired colormap behaviour in 3D plotting, however.)
    gainAmplitudesPol1         = numpy.where(flagsForAmplitudesPol1, numpy.nan, gainAmplitudesPol1)
    gainAmplitudesPol2         = numpy.where(flagsForAmplitudesPol2, numpy.nan, gainAmplitudesPol2)
    gainPhasesPol1             = numpy.where(flagsForPhasesPol1,     numpy.nan, gainPhasesPol1)
    gainPhasesPol2             = numpy.where(flagsForPhasesPol2,     numpy.nan, gainPhasesPol2)


    # Load the field name.
    objectMS                   = lib_ms.MS(pathMS)
    nameField                  = objectMS.getNameField()


    # Determine which antennae are working.
    # Should we differentiate between polarisation 1 and 2, because the data streams are handled independently?
    antennaeWorking            = []
    for i in range(numberOfAntennae):
        if (numpy.all(flagsForAmplitudes[ : , i, : , : ]) and numpy.all(flagsForPhases[ : , i, : , : ])):
            antennaeWorking.append(False)
        else:
            antennaeWorking.append(True)


    # As long as the wrong polarisation names are in the H5Parm file, we replace these ourselves. Temporary!
    namesPolarisation          = ["RR", "LL"]


    # Plot gain amplitudes.
    #plotAmplitudes2D(gainAmplitudesPol1, times, frequencies, antennaeWorking, pathDirectoryPlots, namePolarisation = namesPolarisation[0], nameField = nameField, nameDataSet = pathH5ParmInput)
    #plotAmplitudes2D(gainAmplitudesPol2, times, frequencies, antennaeWorking, pathDirectoryPlots, namePolarisation = namesPolarisation[1], nameField = nameField, nameDataSet = pathH5ParmInput)

    # Plot gain phases.
    #plotPhases2D(    gainPhasesPol1,     times, frequencies, antennaeWorking, pathDirectoryPlots, namePolarisation = namesPolarisation[0], nameField = nameField, nameDataSet = pathH5ParmInput)
    #plotPhases2D(    gainPhasesPol2,     times, frequencies, antennaeWorking, pathDirectoryPlots, namePolarisation = namesPolarisation[1], nameField = nameField, nameDataSet = pathH5ParmInput)


    # Initialise bandpass determination settings.
    bandpassAmplitudeNumberOfSubiterations     = 3
    bandpassAmplitudeFlaggingThresholdFactor   = 3
    bandpassAmplitudeMedianFilteringKernelSize = 7

    bandpassPhaseNumberOfSubiterations         = 3
    bandpassPhaseFlaggingThresholdFactor       = 3
    bandpassPhaseMedianFilteringKernelSize     = 7

    bandpassNaNs                               = numpy.ones_like(frequencies) * numpy.nan
    functionDTECNaNs                           = numpy.ones_like(times) * numpy.nan



    #
    # Generate amplitude bandpasses (in an iterative way).
    #

    # Create lists that store, for each antenna, 4 amplitude bandpasses (2 polarisations, 2 iterations).
    bandpassesAmplitudePol1Iter1 = []
    bandpassesAmplitudePol1Iter2 = []
    bandpassesAmplitudePol2Iter1 = []
    bandpassesAmplitudePol2Iter2 = []

    for i in range(numberOfAntennae):

        # Produce amplitude bandpasses.
        if (antennaeWorking[i]):
            print ("Starting amplitude bandpass calculation for antenna ID " + str(i) + "...")
            bandpassAmplitudePol1Iter1, bandpassAmplitudePol1Iter2 = calculateBandpassAmplitude(gainAmplitudesPol1[i], flagsForAmplitudesPol1[i],
                                                                                                numberOfSubiterations = bandpassAmplitudeNumberOfSubiterations, flaggingThresholdFactor = bandpassAmplitudeFlaggingThresholdFactor, medianFilteringKernelSize = bandpassAmplitudeMedianFilteringKernelSize)
            bandpassAmplitudePol2Iter1, bandpassAmplitudePol2Iter2 = calculateBandpassAmplitude(gainAmplitudesPol2[i], flagsForAmplitudesPol2[i],
                                                                                                numberOfSubiterations = bandpassAmplitudeNumberOfSubiterations, flaggingThresholdFactor = bandpassAmplitudeFlaggingThresholdFactor, medianFilteringKernelSize = bandpassAmplitudeMedianFilteringKernelSize)
        else:
            bandpassAmplitudePol1Iter1, bandpassAmplitudePol1Iter2, bandpassAmplitudePol2Iter1, bandpassAmplitudePol2Iter2 = [bandpassNaNs] * 4

        # Save amplitude bandpasses.
        bandpassesAmplitudePol1Iter1.append(bandpassAmplitudePol1Iter1)
        bandpassesAmplitudePol1Iter2.append(bandpassAmplitudePol1Iter2)
        bandpassesAmplitudePol2Iter1.append(bandpassAmplitudePol2Iter1)
        bandpassesAmplitudePol2Iter2.append(bandpassAmplitudePol2Iter2)


    #
    # Generate phase bandpasses (in an iterative way). Calibrator DTECs are found as a side product.
    #

    # Create lists that store, for each antenna, 4 phase bandpasses (2 polarisations, 2 iterations) and 2 functions DTEC(t) (2 polarisations; should be identical).
    bandpassesPhasePol1Iter1     = []
    bandpassesPhasePol1Iter2     = []
    bandpassesPhasePol2Iter1     = []
    bandpassesPhasePol2Iter2     = []
    functionsDTECPol1            = []
    functionsDTECPol2            = []

    for i in range(numberOfAntennae):

        # Produce phase bandpasses and DTECs.
        if (antennaeWorking[i]):
            print ("Starting phase bandpass (and DTEC) calculation for antenna ID " + str(i) + "...")

            bandpassPhasePol1Iter1, bandpassPhasePol1Iter2, functionDTECPol1 = calculateBandpassPhase(gainPhasesPol1[i], flagsForPhasesPol1[i], frequencies, times,
                                                                                                      numberOfSubiterations = bandpassPhaseNumberOfSubiterations, flaggingThresholdFactor = bandpassPhaseFlaggingThresholdFactor, medianFilteringKernelSize = bandpassPhaseMedianFilteringKernelSize)
            bandpassPhasePol2Iter1, bandpassPhasePol2Iter2, functionDTECPol2 = calculateBandpassPhase(gainPhasesPol2[i], flagsForPhasesPol2[i], frequencies, times,
                                                                                                      numberOfSubiterations = bandpassPhaseNumberOfSubiterations, flaggingThresholdFactor = bandpassPhaseFlaggingThresholdFactor, medianFilteringKernelSize = bandpassPhaseMedianFilteringKernelSize)
        else:
            bandpassPhasePol1Iter1, bandpassPhasePol1Iter2, bandpassPhasePol2Iter1, bandpassPhasePol2Iter2 = [bandpassNaNs] * 4
            functionDTECPol1, functionDTECPol2                                                             = [functionDTECNaNs] * 2

        # Save phase bandpasses and DTECs.
        bandpassesPhasePol1Iter1.append(bandpassPhasePol1Iter1)
        bandpassesPhasePol1Iter2.append(bandpassPhasePol1Iter2)
        bandpassesPhasePol2Iter1.append(bandpassPhasePol2Iter1)
        bandpassesPhasePol2Iter2.append(bandpassPhasePol2Iter2)
        functionsDTECPol1.append(functionDTECPol1)
        functionsDTECPol2.append(functionDTECPol2)



    # Plot amplitude bandpasses.
    plotBandpassesAmplitude(bandpassesAmplitudePol1Iter1, bandpassesAmplitudePol2Iter1, frequencies, antennaeWorking, pathDirectoryPlots, namesPolarisation = namesPolarisation, nameIteration = "1", nameField = nameField, nameDataSet = pathH5ParmInput)
    plotBandpassesAmplitude(bandpassesAmplitudePol1Iter2, bandpassesAmplitudePol2Iter2, frequencies, antennaeWorking, pathDirectoryPlots, namesPolarisation = namesPolarisation, nameIteration = "2", nameField = nameField, nameDataSet = pathH5ParmInput)

    # Plot amplitude bandpass overviews.
    plotBandpassesAmplitude2D(bandpassesAmplitudePol1Iter2, pathDirectoryPlots, namePolarisation = namesPolarisation[0], nameField = nameField, nameDataSet = pathH5ParmInput)
    plotBandpassesAmplitude2D(bandpassesAmplitudePol2Iter2, pathDirectoryPlots, namePolarisation = namesPolarisation[1], nameField = nameField, nameDataSet = pathH5ParmInput)


    # Plot phase bandpasses.
    plotBandpassesPhase(bandpassesPhasePol1Iter1, bandpassesPhasePol2Iter1, frequencies, antennaeWorking, pathDirectoryPlots, namesPolarisation = namesPolarisation, nameIteration = "1", nameField = nameField, nameDataSet = pathH5ParmInput)
    plotBandpassesPhase(bandpassesPhasePol1Iter2, bandpassesPhasePol2Iter2, frequencies, antennaeWorking, pathDirectoryPlots, namesPolarisation = namesPolarisation, nameIteration = "2", nameField = nameField, nameDataSet = pathH5ParmInput)

    # Plot phase bandpass overviews.
    plotBandpassesPhase2D(bandpassesPhasePol1Iter2, pathDirectoryPlots, namePolarisation = namesPolarisation[0], nameField = nameField, nameDataSet = pathH5ParmInput)
    plotBandpassesPhase2D(bandpassesPhasePol2Iter2, pathDirectoryPlots, namePolarisation = namesPolarisation[1], nameField = nameField, nameDataSet = pathH5ParmInput)


    # Plot TECs.


    # Plot TEC overviews.
    plotFunctionsDTEC2D(functionsDTECPol1, pathDirectoryPlots, suffixFilename = "_pol" + namesPolarisation[0], comment = "derived from polarisation " + namesPolarisation[0], nameField = nameField, nameDataSet = pathH5ParmInput)
    plotFunctionsDTEC2D(functionsDTECPol2, pathDirectoryPlots, suffixFilename = "_pol" + namesPolarisation[1], comment = "derived from polarisation " + namesPolarisation[1], nameField = nameField, nameDataSet = pathH5ParmInput)



    # Create final data products.
    cubeBandpassAmplitudeValues  = numpy.array([bandpassesAmplitudePol1Iter2, bandpassesAmplitudePol2Iter2])
    cubeBandpassAmplitudeWeights = numpy.logical_not(numpy.isnan(cubeBandpassAmplitudeValues))

    cubeBandpassPhaseValues      = numpy.array([bandpassesPhasePol1Iter2, bandpassesPhasePol2Iter2])
    cubeBandpassPhaseWeights     = numpy.logical_not(numpy.isnan(cubeBandpassPhaseValues))

    # Calculate final DTECs by averaging, per antenna, the DTECs found for the 2 polarisations: both polarisations should measure the same DTEC.
    gridTECValues                = []
    for functionDTECPol1, functionDTECPol2 in zip(functionsDTECPol1, functionsDTECPol2):
        gridTECValues.append(numpy.nan_to_num(fillGaps1D(numpy.add(functionDTECPol1, functionDTECPol2) / 2)))
    gridTECValues                = numpy.array(gridTECValues)
    gridTECWeights               = numpy.ones_like(gridTECValues)



    # Save the amplitude and phase bandpasses.
    lib_util.check_rm(pathH5ParmOutput)
    objectH5Parm = h5parm.h5parm(pathH5ParmOutput, readonly = False)
    objectSolSet = objectH5Parm.makeSolset(solsetName = "sol000", addTables = False)
    objectSolSet.makeSoltab(soltype = "amplitude", soltabName = "bandpassAmplitude", axesNames = ["pol", "ant", "freq"], axesVals = [namesPolarisation, namesAntenna, frequencies], vals = cubeBandpassAmplitudeValues, weights = cubeBandpassAmplitudeWeights)
    objectSolSet.makeSoltab(soltype = "phase",     soltabName = "bandpassPhase",     axesNames = ["pol", "ant", "freq"], axesVals = [namesPolarisation, namesAntenna, frequencies], vals = cubeBandpassPhaseValues,     weights = cubeBandpassPhaseWeights)
    objectSolSet.makeSoltab(soltype = "tec",       soltabName = "TEC",               axesNames = ["ant", "time"],        axesVals = [namesAntenna, times],                          vals = gridTECValues,               weights = gridTECWeights)
    objectH5Parm.close()



if (__name__ == "__main__"):
    # If the program is run from the command line, parse arguments.
    parser                      = argparse.ArgumentParser(description = "Pipeline step 3: Generation of bandpasses.")
    parser.add_argument("pathDirectoryMS", help = "Path of the directory that contains the calibrator scan MS (directory).")
    arguments                   = parser.parse_args()

    # Temporary!
    arguments.pathDirectoryMS = "/disks/strw3/oei/uGMRTCosmosCut-PiLF/fieldsCalibrator/scanID1"

    lib_util.printLineBold("Parameters to use:")
    print (arguments)

    dedicated_uGMRT_bandpass(arguments.pathDirectoryMS)



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
#     nameReferenceAntenna        = "C13" # user choice, e.g. "C00" or "C13"
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
#erivTimePol1, derivTimePol2]
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