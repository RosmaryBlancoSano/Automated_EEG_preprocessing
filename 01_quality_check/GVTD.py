import numpy as np
from enum import Enum
import matplotlib.pyplot as plt


def gvtd_processing(dataMatrix):
    # Step 1: Find the matrix of the temporal derivatives
    # dataMatrix = dataMatrix.T
    dataDiff = dataMatrix - np.roll(dataMatrix, -1, axis=1)

    # Step 2: Find the RMS across the channels for each time-point of dataDiff
    # gvtdTimeTrace = np.sqrt(np.mean(np.square(dataDiff[:, :-1]), axis=1))
    # gvtdTimeTrace = np.sqrt(np.mean(np.diff(dataDiff[:, :-1], axis=1) ** 2, axis=0))
    gvtdTimeTrace = np.sqrt(np.mean(dataDiff[:, :-1] ** 2, axis=0))

    # Step 3: Add a zero in the beginning for GVTD to have the same number of
    # time-points as your original dataMatrix
    gvtdTimeTrace = np.concatenate(([0], gvtdTimeTrace))

    return gvtdTimeTrace

class StatType(Enum):
    Default = 0
    Histogram_Mode = 1
    Kdensity_Mode = 2
    Parabolic_Mode = 3
    Mean = 4
    Median = 5
    MAD = 6

def find_gvtd_thresh(gvtdTimeTrace, statType = None, nStd = None):
    if statType is None:
        statType = 'Default'

    if statType != 'Default' and nStd is None:
        nStd = 10

    if gvtdTimeTrace.ndim == 2 and gvtdTimeTrace.shape[1] == 1:
        gvtdTimeTrace = gvtdTimeTrace.flatten()

    if statType == 'Default':
        minCountsEachBin = 5
        nBins = int(np.round(len(gvtdTimeTrace) / minCountsEachBin))
        del_ = np.max(gvtdTimeTrace) / nBins
        N, edges = np.histogram(gvtdTimeTrace, bins=nBins, range=(0, np.max(gvtdTimeTrace)))
        argmax = np.argmax(N)
        runMode = edges[argmax] + del_ / 2
        minValToModeDist = runMode - np.min(gvtdTimeTrace)
        thresh = runMode + minValToModeDist

    elif statType == 'Histogram_Mode':
        minCountsEachBin = 5
        nBins = int(np.round(len(gvtdTimeTrace) / minCountsEachBin))
        del_ = np.max(gvtdTimeTrace) / nBins
        N, edges = np.histogram(gvtdTimeTrace, bins=nBins, range=(0, np.max(gvtdTimeTrace)))
        argmax = np.argmax(N)
        runMode = edges[argmax] + del_ / 2
        minValToModeDist = runMode - np.min(gvtdTimeTrace)
        thresh = runMode + nStd * minValToModeDist

    elif statType == 'Mean':
        runMean = np.mean(gvtdTimeTrace)

        numPointsBelowMode = np.sum(gvtdTimeTrace < runMean)
        rmsPointBelowMode = np.sum((gvtdTimeTrace[gvtdTimeTrace < runMean] - runMean) ** 2)
        leftStdRun = np.sqrt(rmsPointBelowMode / numPointsBelowMode)
        thresh = runMean + nStd * leftStdRun

    elif statType == 'Median':
        runMedian = np.median(gvtdTimeTrace)

        numPointsBelowMedian = np.sum(gvtdTimeTrace < runMedian)
        rmsPointBelowMedian = np.sum((gvtdTimeTrace[gvtdTimeTrace < runMedian] - runMedian) ** 2)
        leftStdRun = np.sqrt(rmsPointBelowMedian / numPointsBelowMedian)

        thresh = runMedian + nStd * leftStdRun;

    elif statType == 'MAD':
        runMAD = np.median(np.abs(gvtdTimeTrace - np.median(gvtdTimeTrace)))
        thresh = np.median(gvtdTimeTrace) + nStd * runMAD

    elif statType == 'Left_Std':
        runMode = np.histogram_bin_edges(gvtdTimeTrace, bins='fd')[0]
        minValToModeDist = runMode - np.min(gvtdTimeTrace)
        runStd = np.std(gvtdTimeTrace[gvtdTimeTrace < runMode])
        thresh = runMode + nStd * runStd

    return thresh


def make_gvtd_hist(gvtdTimeTrace, plotThresh=None, statType=None, nStd=None, binSize=None):
    if plotThresh is None:
        plotThresh = 1

    if statType is None:
        statType = "Histogram_Mode"

    if binSize is None:
        minCountsEachBin = 5
        nBins = int(np.round(gvtdTimeTrace.shape[0] / minCountsEachBin))
        binSize = max(gvtdTimeTrace) / nBins

    plt.hist(gvtdTimeTrace, bins=np.arange(0, max(gvtdTimeTrace) + binSize, binSize))

    if plotThresh:
        if nStd is None:
            nStd = 4
        thresh = find_gvtd_thresh(gvtdTimeTrace, statType, nStd)
        hist, bins = np.histogram(gvtdTimeTrace, bins=np.arange(0, max(gvtdTimeTrace) + binSize, binSize))
        plt.plot([thresh, thresh], [0, max(hist)])

    plt.title("GVTD Histogram")
    plt.xlabel("GVTD")
    plt.ylabel("Counts")
    plt.show()