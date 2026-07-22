# Package index

## Call density estimation

The main entry points – the functions most users will call directly.

- [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  : Call Density Estimate using auxiliary information from sonar
  equation
- [`cdeFromParamFile()`](https://briansmiller.github.io/callDensity/reference/cdeFromParamFile.md)
  : Call density estimate from a parameter file (deprecated)
- [`callDensity`](https://briansmiller.github.io/callDensity/reference/callDensity-package.md)
  [`callDensity-package`](https://briansmiller.github.io/callDensity/reference/callDensity-package.md)
  : callDensity: Analysis for bioacoustic call density estimation
  (ABCDE)

## Uncertainty and coefficients of variation

Propagating and combining error terms into CV.Nc, CV.c, CV.pa, CV.Dc.

- [`Dc_CV()`](https://briansmiller.github.io/callDensity/reference/Dc_CV.md)
  : Coefficient of Variation (CV) for call density (Dc)
- [`Nc_CV()`](https://briansmiller.github.io/callDensity/reference/Nc_CV.md)
  : Nc_CV: Coefficient of Variation (CV) of the number of detected calls
  (Nc).
- [`c_CV()`](https://briansmiller.github.io/callDensity/reference/c_CV.md)
  : c_CV: Coefficient of Variation (CV) of false discovery rate (c).
- [`pa_CV()`](https://briansmiller.github.io/callDensity/reference/pa_CV.md)
  : pa_CV: Coefficient of Variation (CV) of the probability of detection
  (pa).
- [`ciFromCV()`](https://briansmiller.github.io/callDensity/reference/ciFromCV.md)
  : Confidence interval from density and it's CV
- [`falseDiscoveryCV()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryCV.md)
  : WARNING: This function is not supported, and probably does not do
  what whatever you were hoping it might do. Perhaps consider function
  c_CV instead. CV of false discovery rate using Cochran approximation
- [`var.wtd.mean.cochran()`](https://briansmiller.github.io/callDensity/reference/var.wtd.mean.cochran.md)
  : Variance of a weighted mean following Chochran 1977 definition

## Detection functions

Fitting, predicting, and displaying SNR-based detection curves –
glm/gam/scam for observer-ground designs, vglm for adjudicated
capture-recapture.

- [`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md)
  : Fit an SNR detection function
- [`predictDetFun()`](https://briansmiller.github.io/callDensity/reference/predictDetFun.md)
  : Predict detection function from fitted SNR models
- [`predictDetFunList()`](https://briansmiller.github.io/callDensity/reference/predictDetFunList.md)
  : Predict detection functions for multiple models
- [`predict(`*`<detFun>`*`)`](https://briansmiller.github.io/callDensity/reference/predict.detFun.md)
  : predict method for detFun objects
- [`showDetFun()`](https://briansmiller.github.io/callDensity/reference/showDetFun.md)
  : Plot SNR detection function with observed SNR distributions
- [`fitSNRdetectionFunc()`](https://briansmiller.github.io/callDensity/reference/fitSNRdetectionFunc.md)
  : Fit an SNR-detection function
- [`fitSNRvgam()`](https://briansmiller.github.io/callDensity/reference/fitSNRvgam.md)
  : Fit an SNR-detection function with a closed population
  capture-recapture gam
- [`fitSNRvglm()`](https://briansmiller.github.io/callDensity/reference/fitSNRvglm.md)
  : Fit an SNR-detection function with a closed population
  capture-recapture glm
- [`fitSNRbySeason()`](https://briansmiller.github.io/callDensity/reference/fitSNRbySeason.md)
  : Model SNR-detection function with season as a factor (deprecated?)
- [`predVglmPDet()`](https://briansmiller.github.io/callDensity/reference/predVglmPDet.md)
  : Predict probability of detection from a VGLM SNR-detection function
- [`showSNRdetectionFunc()`](https://briansmiller.github.io/callDensity/reference/showSNRdetectionFunc.md)
  : Plot an SNR detection function with mirrored SNR distributions

## Probability of detection and the sonar equation

- [`pDetInArea()`](https://briansmiller.github.io/callDensity/reference/pDetInArea.md)
  : Monte-Carlo simulation to predict average probability of detection
  in area
- [`pDetGivenNL()`](https://briansmiller.github.io/callDensity/reference/pDetGivenNL.md)
  : Probability of detection in the study area at a fixed noise level
- [`studyArea()`](https://briansmiller.github.io/callDensity/reference/studyArea.md)
  : Calculate study area for call density estimate.

## Noise level estimation

- [`nlFromDetections()`](https://briansmiller.github.io/callDensity/reference/nlFromDetections.md)
  : Estimate the noise level distribution from noise measured at
  detections
- [`nlFromSnrInfo()`](https://briansmiller.github.io/callDensity/reference/nlFromSnrInfo.md)
  : Estimate noise level distribution for MC simulation from the noise
  measurements included in the snrInfo file, plus the SNR detection
  function.
- [`noiseLevelDistribution()`](https://briansmiller.github.io/callDensity/reference/noiseLevelDistribution.md)
  : Mean and standard deviation of noise levels by month, season, or
  year.
- [`predictSampledNL()`](https://briansmiller.github.io/callDensity/reference/predictSampledNL.md)
  : Predict the mean noise level you would measure at detections
- [`kerguelen2015TL`](https://briansmiller.github.io/callDensity/reference/kerguelen2015TL.md)
  : pyRAM transmission loss, Kerguelen2015 summer, 25 Hz
- [`kerguelen2015TLmeta`](https://briansmiller.github.io/callDensity/reference/kerguelen2015TLmeta.md)
  : Transect metadata for the Kerguelen2015 pyRAM transmission loss
  table

## Capture history construction

Building and reshaping the tables cde() and the detection-function
functions expect, including the adjudicated-verdict convention CR
analyses require.

- [`capHist2snrInfo()`](https://briansmiller.github.io/callDensity/reference/capHist2snrInfo.md)
  : Convert a capture history table into an SNRinfo data.frame
- [`capHistTosnrInfo()`](https://briansmiller.github.io/callDensity/reference/capHistTosnrInfo.md)
  : Convert capture history DATA.FRAME into the 'SNRinfo' format used by
  the callDensity package capHistTosnrInfo
- [`readCapHist()`](https://briansmiller.github.io/callDensity/reference/readCapHist.md)
  : Read a capture history csv file (e.g. created in Matlab)
- [`simsTocaptureHistoryTable()`](https://briansmiller.github.io/callDensity/reference/simsTocaptureHistoryTable.md)
  : Create a capture history table from the two simulated detection
  tables
- [`mchToCR()`](https://briansmiller.github.io/callDensity/reference/mchToCR.md)
  : Convert a multi-observer capture history table to a two-observer CR
  table
- [`resolveColumns()`](https://briansmiller.github.io/callDensity/reference/resolveColumns.md)
  : Resolve column names from a prefix pattern or explicit vector
- [`pivotSNR()`](https://briansmiller.github.io/callDensity/reference/pivotSNR.md)
  : Pivot SNR observer columns to long format
- [`addObserverMeans()`](https://briansmiller.github.io/callDensity/reference/addObserverMeans.md)
  : Compute row-wise mean across observer columns
- [`multiObserverDetectionCount()`](https://briansmiller.github.io/callDensity/reference/multiObserverDetectionCount.md)
  : Multi-Observer Detection Count

## False discovery rate

- [`falseDiscoveryRate()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryRate.md)
  : Estimate \$c\$ (False discovery rate)
- [`falseDiscoveryRateFromNth()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryRateFromNth.md)
  : WARNING: This function is not supported, and probably does not do
  what whatever you were hoping it might do. Use functions
  falseDiscoverRate and c_CV instead. False discovery rate from
  inspection of every Nth detection

## Simulation

Generating synthetic calls, detectors, and detections – used throughout
the vignettes and test suite.

- [`simCallLocation()`](https://briansmiller.github.io/callDensity/reference/simCallLocation.md)
  : Simulate animals calls using a uniform distribution time and space
- [`simCallAcoustics()`](https://briansmiller.github.io/callDensity/reference/simCallAcoustics.md)
  : Simulate acoustic properties of calls
- [`simulateDetector()`](https://briansmiller.github.io/callDensity/reference/simulateDetector.md)
  : Simulate a detector for a callDensity simulation.
- [`simTLradials_20logR()`](https://briansmiller.github.io/callDensity/reference/simTLradials_20logR.md)
  : Simulate transmission loss (TL) radials following geometric
  (spherical) spreading law.
- [`subsampleSimInTime()`](https://briansmiller.github.io/callDensity/reference/subsampleSimInTime.md)
  : Subsample from a simulation at evenly spaced time intervals.

## Time and season utilities

- [`capHistTimeSeason()`](https://briansmiller.github.io/callDensity/reference/capHistTimeSeason.md)
  : Add R datetime, season, and month codes to a capture history table
- [`capHistTimeSeason2()`](https://briansmiller.github.io/callDensity/reference/capHistTimeSeason2.md)
  : Add R datetime and season to a capture history table (version 2)
- [`mat2Rdate()`](https://briansmiller.github.io/callDensity/reference/mat2Rdate.md)
  : Convert matlab datenum to R POSIXct
- [`time2monthCode()`](https://briansmiller.github.io/callDensity/reference/time2monthCode.md)
  : Lookup month codes from a POSIXct
- [`time2season()`](https://briansmiller.github.io/callDensity/reference/time2season.md)
  : Lookup the season for a given POSIXct, x.
- [`subsetByTimeCode()`](https://briansmiller.github.io/callDensity/reference/subsetByTimeCode.md)
  : Subset a data.frame by by months or season
- [`deploymentDuration()`](https://briansmiller.github.io/callDensity/reference/deploymentDuration.md)
  : Calculate deployment duration from recording start times and
  durations
- [`deploymentDurationFromsoundFolderCsv()`](https://briansmiller.github.io/callDensity/reference/deploymentDurationFromsoundFolderCsv.md)
  : Calculate deployment duration from a Matlab soundFolder csv file.

## Workflow and file helpers

- [`defaultOutputFileNames()`](https://briansmiller.github.io/callDensity/reference/defaultOutputFileNames.md)
  : Default names of output files for call density estimation
- [`listTLFiles()`](https://briansmiller.github.io/callDensity/reference/listTLFiles.md)
  : Return a list of existing TL files associated with the call density
- [`updateDetectionFolder()`](https://briansmiller.github.io/callDensity/reference/updateDetectionFolder.md)
  : Change path to detection csv files
- [`countDetections()`](https://briansmiller.github.io/callDensity/reference/countDetections.md)
  : Title
- [`densityInputTable()`](https://briansmiller.github.io/callDensity/reference/densityInputTable.md)
  : Show table of inputs into call density estimate
- [`densityResultsTable()`](https://briansmiller.github.io/callDensity/reference/densityResultsTable.md)
  : Table of call densities including CVs (ggplot2)
- [`collectDensityResults()`](https://briansmiller.github.io/callDensity/reference/collectDensityResults.md)
  : Load specified density results text files in a directory

## Plotting

- [`densityPlot()`](https://briansmiller.github.io/callDensity/reference/densityPlot.md)
  : Bar-plot of call density by timeCodes (ggplot2)
- [`detectionRateCorrectedPlot()`](https://briansmiller.github.io/callDensity/reference/detectionRateCorrectedPlot.md)
  : Bar-plot of detection rates by timeCodes and scaled (corrected) by
  precision (ggplot2)
- [`detectionRatePlot()`](https://briansmiller.github.io/callDensity/reference/detectionRatePlot.md)
  : Bar-plot of detection rates by timeCodes (ggplot2)
- [`plotDetectionDistribution()`](https://briansmiller.github.io/callDensity/reference/plotDetectionDistribution.md)
  : SNR histogram of true & false positives, and false negatives for
  callDensity simulation
- [`plotSNRHistogram()`](https://briansmiller.github.io/callDensity/reference/plotSNRHistogram.md)
  : Plot SNR histogram by observer
- [`plotSNRTimeSeries()`](https://briansmiller.github.io/callDensity/reference/plotSNRTimeSeries.md)
  : Plot SNR time series for the adjudicated subset
- [`plotSpatialDetections()`](https://briansmiller.github.io/callDensity/reference/plotSpatialDetections.md)
  : Plot the spatial detection density for a callDensity simulation
- [`computeSNRLims()`](https://briansmiller.github.io/callDensity/reference/computeSNRLims.md)
  : Compute pretty SNR limits across observer columns

## Internal utilities

- [`tidyeval`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  [`enquo`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  [`enquos`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  [`.data`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  [`:=`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  [`as_name`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  [`as_label`](https://briansmiller.github.io/callDensity/reference/tidyeval.md)
  : Tidy eval helpers
