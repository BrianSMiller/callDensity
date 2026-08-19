# library("stats",quietly = T, warn.conflicts = F)
# library("dplyr",quietly = T, warn.conflicts = F)
# library("ggplot2",quietly=T,warn.conflicts = F)
# library("kableExtra", quietly=T, warn.conflicts = F)
# library("lubridate", quietly = T, warn.conflicts = F)
# library("Hmisc",quietly = T, warn.conflicts = F)

#' Convert matlab datenum to R POSIXct
#'
#' @param x List of Matlab datenums
#' @returns posixCT corresponding to the Matlab datenum that was input
#' @export
mat2Rdate <- function (x){ as.POSIXct((x - 719529)*86400, origin = "1970-01-01",
                                     tz = "UTC")}

#' Lookup the season for a given POSIXct, x.
#'
#' @description
#' Seasons are defined by 3 month
#' periods from Summer=(Dec,Jan,Feb); Autumn=(Mar,Apr,May), Winter=(Jun,Jul,Aug),
#' Spring=(Sep,Oct,Nov)
#'
#' @param x POSIXct
#' @returns Factor containing seasons summer, autumn, winter, or spring
#'
#' @export
time2season <- function(x) {
  mon<-lubridate::month(x)
  seasons <- rep(c("winter","spring","summer","autumn"), each = 3)
  names(seasons) <- month.abb[c(6:12,1:5)]
  season <- factor(seasons[month.abb[mon]],
                   levels = c("summer","autumn","winter","spring"))
}

#' Lookup month codes from a POSIXct
#'
#' @description
#' Lookup the month codes for a given POSIXct, x. MonthCodes are a factor from
#' '01'-'12' respectively for months from Jan-Dec
#'
#' @param x POSIXct
#' @returns Factor containing month codes '01'-'12' for months from Jan-Dec.
#'
#' @export
time2monthCode <- function(x) {
  mon<-lubridate::month(x)
  monthCodes <- sprintf("%02d",mon)
  monthCodes <- factor(monthCodes,
                   levels = c('01','02','03','04','05','06',
                              '07','08','09','10','11','12'))
}

#' Add R datetime, season, and month codes to a capture history table
#'
#' @param x Capture history table with column t0_tableN of Matlab datenum's
#' @returns Capture history table with column t, season, and month containing
#'     POSIXct, and timeCode factors for season and month.
#' @export
capHistTimeSeason <- function(x){
  x$t <- mat2Rdate(x$t0_table1)
  x$t[is.na(x$t)] <- mat2Rdate(x$t0_table2[is.na(x$t)])
  x$season <- time2season(x$t)
  x$month <- time2monthCode(x$t)
  # }
  return(x)
}

#' Add R datetime and season to a capture history table (version 2)
#'
#' @param x A capture history table with column t0 containing a matlab datenum
#'
#' @returns A capture history table with original columns plus:
#'     t: a posixCT corresponding to the date and time of t0
#'     season: a factor containing summer, autumn, winter, spring, or year
#'       corresponding to the season of t0
#'     month: a factor containing '01'-'12' corresponding to the month of t0
#' @export
capHistTimeSeason2 <- function(x){
  x$t <- mat2Rdate(x$t0)
  x$season <- time2season(x$t)
  x$month <- time2monthCode(x$t)
  return(x)
}

#' Subset a data.frame by by months or season
#'
#' @description
#' Subset a data.frame by timeCode, where timeCode can be either a 'season'
#' (summer, autumn, winter, spring), a month characters ('01'-'12'), or 'year',
#' with 'year' the same as including all data
#'
#' @param df A data.frame
#' @param dt A posixCT datetime
#' @param timeCode A timecode indicating the season, month, or 'year' for all data
#'
#' @export
subsetByTimeCode<- function(df, dt, timeCode){

  if (timeCode !='year'){
    if (timeCode == 'summer' || timeCode=='autumn' ||
        timeCode=='winter' || timeCode=='spring'){
      times <- time2season(dt)
    } else {
      times <- time2monthCode(dt)
    }
    df <- df[times==timeCode,]
  }
  return(df)
}

#' Default names of output files for call density estimation
#'
#' @description
#' Given data.frame of call density parameters, p, create default parameters and
#' output file-names for Monte-Carlo modelling and estimating call-density.
#' Append all of these parameters to the data.frame and return it.
#'
#' @param p Data.frame containing parameters for call density estimation
#' @param season TimeCode to filter month, season, or full year
#' @param outputFolder Location/path/folder where output files will be saved
#'   default is '.' (the current working directory).
#'
#' @returns Parameter file with file names corresponding to the correct timeCode
#'
#' @export
defaultOutputFileNames <- function(p,season, outputFolder='.'){

  ### Inputs
  ## Density equation
  p$k <- 1; # Number of sensors
  p$w <- 1000; # Defined radius of simulation area (km)

  # Set to TRUE if false-positive automated detections have been independently
  # verified from manual inspection of every Nth detection
  p$useSeparateFPdata <- FALSE;
  p$manualFPfileName <- 'falsePositiveProportion_50annotations_Kerguelen2014_FRC.xlsx';

  ## P_a & Monte-Carlo model

  # Source level distribution
  # p$SLmean <- 189;
  # p$SLsd <- 8.0;
  # p$SLsamplesize <- 350;

  # Whether to use GLM, GAM, SCAM, or VGLM (not yet implemented)
  # p$useGLM <- TRUE;
  # Whether to constrained gam from scam or plain mgcv::gam
  # p$useSCAM <- TRUE;
  p$modelType <- 'scam';
    # Number of knots used in GAM
  p$numKnots <- 3

  #number of iterations for the outer loop
  p$outerloop<-1000

  # range resolution in m for detection probabilities along  each transect
  p$output.resolution.m <- 100

  # Seasonal TL file name
  p$tlFile <- paste0('TL/TL_', siteCode, '_woa2018_meanSSP_',
                     p$output.resolution.m, 'm_res_', season, '.csv')
  p$tlFile <- ifelse(file.exists(p$tlFile), p$tlFile,
                     paste0('TL/TL_', siteCode, '_woa2018_meanSSP_', season, '.csv') )

  ### Outputs
  ## Output file names of Monte-Carlo model of Pa
  p$paFile <- file.path(outputFolder,
                        paste0("Pa_", p$siteCode, "_", season, ".txt") )
  p$transectFile <- file.path(outputFolder,
                        paste0("pDet_",p$siteCode,'_',season,'_transects.csv') )
  p$simResultsFile <- file.path(outputFolder,
                        paste0('simResult_', p$siteCode, '_',season,'.txt') )

  ## Overall call density output file name
  p$densityResultsFile <- file.path(outputFolder,
                        paste0('density_',p$siteCode, '_', season,'.csv') )

  p$outputFolder <- outputFolder
  if (!dir.exists(outputFolder)){
    dir.create(outputFolder)
  }

  return(p)
}

#' Return a list of existing TL files associated with the call density
#'
#' @param p Data.frame containing parameters for call density estimation
#' @param season TimeCode to filter month, season, or full year
#'
#' @returns a list of all transmission loss files
#' @export
listTLFiles <- function(p,season=''){
  # Seasonal TL file name
  tlFiles <- dir(p$tlParams$tl_root, pattern=paste0(season,'*.csv'))
  tlFiles <- paste0(p$tlParams$tl_root,tlFiles)

  return(tlFiles)
}

#' Change path to detection csv files
#'
#' @param p Call density detector parameter data.frame
#' @param newFolder relative or absolute path to where the detection files are
#'   located
#'
#' @returns new detection parameters with filenames that have updated folders
#' @export
updateDetectionFolder <- function(p, newFolder){
  d <- p$detectorParams
  d$folder <- newFolder
  d$fullYearDetectionCsv <- file.path(newFolder,basename(d$fullYearDetectionCsv))
  d$subsetDetCsv <- file.path(newFolder,basename(d$subsetDetCsv))
  d$fullYearEffortFile <- file.path(newFolder, basename(d$fullYearEffortFile))
  p$annotationParams$annotationFile <- file.path(
    newFolder, basename(p$annotationParams$annotationFile) )
  p$detectorParams <- d
  return(p)
}



#' @title Mean and standard deviation of noise levels by month, season, or year.
#'
#' @description
#' Read a table of noise levels containing a columns t and NL for time and
#' level, respectively.
#' @returns The mean and standard deviation from the table of noise levels.
#'
#' @param nlFile File name of csv containing noise level times and values
#'
#' @param season TimeCode for season, month, or 'year'
#'
#' @export
noiseLevelDistribution <- function(nlFile,season='year'){
  nl <- utils::read.csv(nlFile)
  nl$t <- lubridate::dmy_hms(nl$t)
  nl$season <- time2season(nl$t)
  time2season(nl$t)
  NLseason <- nl$nl
  if (season != 'year'){
    NLseason <- nl$nl[nl$season==season]
  }
  sampleSize<-sum(!is.na(NLseason)) #(needed for step 4b)
  mean<- mean(NLseason,na.rm=TRUE) # in dB
  sd<- sd(NLseason,na.rm = TRUE)
  NL <- data.frame(mean,sd,sampleSize)
  return(NL)
}


#' Estimate noise level distribution for MC simulation from the noise
#' measurements included in the snrInfo file, plus the SNR detection function.
#'
#' @param snrInfo - Table of SNR information containing column named NoiseRL
#'   with noise level measurements in dB
#' @param snrDetFun - SNR detection function (e.g. from fitSNRdetectionFunc or
#'   fitSNRvglm)
#'
#' @returns Data.frame with 1 row and 3 columns containing parameterised
#'   distribution of noise levels. Column names are mean, sd, and sampleSize.
#'   Distribution assumed to be normal.
#' @export
nlFromSnrInfo <- function(snrInfo, snrDetFun){

  NL = snrInfo %>% dplyr::summarise(
    mean=mean(NoiseRL,na.rm = TRUE),
    sd=sd(NoiseRL,na.rm = TRUE),
    sampleSize=dplyr::n()-sum(is.na(NoiseRL))
  )

  # Use the SNR intercept estimate from the VGLM to add back into the NL
  # distribution. Have resulted to a clunky kludge since I can't figure out how to
  # get the intercept on the response scale. Instead, just predict the
  # probabilities on the response scale for a grid of SNR values and lookup the
  # median value.
  snrGrid   <- data.frame(SNR = seq(from = -20, to = 20, by = 0.1))
  probs     <- predictDetFun(snrDetFun, newdata = snrGrid)$fit
  snrOffset <- snrGrid[which.min(abs(probs - 0.5)), ]
  NL$mean   <- NL$mean + snrOffset
  return(NL)
}


#' Convert a capture history table to SNRinfo format
#'
#' @description
#' Canonical converter from a capture history table to the 'SNRinfo' format
#' used throughout callDensity for detection-function fitting
#' (\code{\link{fitDetFun}}) and \code{\link{cde}}. Replaces
#' \code{capHist2snrInfo}/\code{capHistTosnrInfo} (now deprecated thin
#' wrappers around this function -- see \code{\link{mchToCR}}, also
#' deprecated, for the older 2-observer-reduction route these superseded).
#'
#' Supports both:
#' \itemize{
#'   \item \strong{OG} (observer-ground-truth): a single \code{observers}
#'     column, treated as the detector under evaluation against a
#'     (possibly imperfect) \code{groundTruth} column. This is the
#'     original, degenerate-format shape (Danielle Harris' original
#'     script): every retained row is treated as a real call, and
#'     \code{Detected} records whether one automated detector also caught
#'     it.
#'   \item \strong{CR} (capture-recapture): two or more \code{observers}
#'     columns, fit jointly via
#'     \code{fitDetFun(modelType = 'vglm', yColNames = c(groundTruthCol, observerCols))}.
#'     Any number of observers is supported -- \code{VGAM::posbernoulli.t}
#'     is not limited to two occasions, despite most examples in this
#'     package using two.
#' }
#'
#' Row selection is a union filter across \code{groundTruth} and all
#' \code{observers} columns (flagged by at least one), followed by
#' restriction to \code{groundTruth == 1} (i.e. adjudicated/assumed-true
#' calls only). This is the sample construction \code{posbernoulli.t}
#' (the CR/vglm case) assumes generated the data -- it conditions on
#' detection by at least one occasion internally, so handing it a sample
#' that was itself assembled that way is correct, not circular. The same
#' construction is also correct, and simpler, for the OG/glm/gam/scam
#' case. What is NOT safe is filtering to \code{groundTruth == 1} and then
#' fitting an ordinary \code{glm} on a \code{vglm}-style union-filtered
#' sample -- \code{glm} has no equivalent internal correction for that
#' selection mechanism, and doing so silently biases the fit (see
#' \code{notes/cde_nl_and_truncation_handover.md}'s ablation section for a
#' worked example of exactly this mistake).
#'
#' @param cht Capture history table. Matchbox-native column names
#'   (\code{t0}, \code{signalRMSdB}, \code{noiseRMSdB}, \code{verdict},
#'   \code{detect_observerN}) are assumed by default. Legacy
#'   \code{detect_table1}/\code{detect_table2} naming also works directly
#'   without any renaming step -- e.g. \code{groundTruth = 'table1'},
#'   \code{observers = 'table2'} reproduces the old
#'   \code{capHist2snrInfo}/\code{capHistTosnrInfo} row selection exactly
#'   (a single-element \code{observers} degenerates algebraically to
#'   their \code{groundTruth == 1} filter, since the union filter with a
#'   single observer is redundant with the ground-truth restriction that
#'   follows it).
#' @param groundTruth Column name (or bare suffix, resolved as
#'   \code{detect_<groundTruth>} if not found directly) to treat as
#'   ground truth. Default \code{'verdict'}.
#' @param observers Character vector of one or more column names (or bare
#'   suffixes, same resolution as \code{groundTruth}) for the detector(s)
#'   under evaluation. A single element reproduces the old OG/SNRinfo
#'   shape; two or more is a genuine multi-observer CR sample.
#' @param signalCol,noiseCol Column(s) supplying CallRL/NoiseRL. Default
#'   \code{'signalRMSdB'}/\code{'noiseRMSdB'} (matchbox-native: one shared
#'   value per matched event, SNR measured after matching). A length-2
#'   vector instead averages a per-observer pair (legacy shape, SNR
#'   measured before matching -- see \code{notes/cde_pipeline_handover.md}'s
#'   "SNR join timing" note).
#' @param timeCol Column to derive \code{t}/\code{month}/\code{season}
#'   from. Default \code{'t0'} (matchbox-native, a Matlab datenum,
#'   converted via \code{\link{mat2Rdate}}). Any other column is used
#'   as-is (assumed already POSIXct-like) -- e.g. \code{timeCol = 't'}
#'   for the legacy pre-computed-time shape.
#' @param season A timeCode corresponding to months, seasons, or 'year'
#'   (the default, which returns all rows).
#'
#' @returns SNRinfo data.frame: the resolved \code{groundTruth} and all
#'   \code{observers} columns (kept intact, unrenamed -- ready to pass as
#'   \code{yColNames} to \code{fitDetFun(modelType = 'vglm')}), plus
#'   \code{Detected} (the last-named observer, for the single-observer
#'   OG case and backward compatibility), \code{CallRL}, \code{NoiseRL},
#'   \code{SNR}, \code{t}, \code{month}, and \code{season}.
#' @export
chtToSNRinfo <- function(cht,
                          groundTruth = 'verdict',
                          observers,
                          signalCol = 'signalRMSdB',
                          noiseCol  = 'noiseRMSdB',
                          timeCol   = 't0',
                          season    = 'year') {

  resolveCol <- function(name, df) {
    if (name %in% names(df)) return(name)
    alt <- paste0('detect_', name)
    if (alt %in% names(df)) return(alt)
    stop("chtToSNRinfo: column '", name, "' not found (also tried '", alt,
         "'). Pass the exact column name if neither matches.")
  }

  gtCol   <- resolveCol(groundTruth, cht)
  obsCols <- vapply(observers, resolveCol, character(1), df = cht)

  # Union filter (flagged by ground truth or any observer under
  # evaluation), then restrict to ground-truth-positive rows. See the
  # function-level @description for why this order, and why it's correct
  # for both OG and CR rather than a compromise between them.
  keep <- Reduce(`|`, as.list(cht[c(gtCol, obsCols)]))
  cht  <- cht[keep, ]
  cht  <- cht[cht[[gtCol]] == 1, ]

  CallRL  <- if (length(signalCol) == 1) as.numeric(cht[[signalCol]])
             else as.numeric(rowMeans(cht[, signalCol], na.rm = TRUE))
  NoiseRL <- if (length(noiseCol) == 1) as.numeric(cht[[noiseCol]])
             else as.numeric(rowMeans(cht[, noiseCol], na.rm = TRUE))
  SNR <- CallRL - NoiseRL

  if (!timeCol %in% names(cht)) {
    stop("chtToSNRinfo: time column '", timeCol, "' not found.")
  }
  t <- if (timeCol == 't0') mat2Rdate(cht[[timeCol]]) else cht[[timeCol]]

  SNRinfo <- data.frame(cht[, c(gtCol, obsCols), drop = FALSE],
                        Detected = cht[[utils::tail(obsCols, 1)]],
                        CallRL = CallRL, NoiseRL = NoiseRL, SNR = SNR,
                        t = t,
                        month  = time2monthCode(t),
                        season = time2season(t))

  infRows <- which(is.infinite(SNRinfo$SNR))
  if (length(infRows) > 0) {
    warning(sprintf('%g of the SNR were infinite:\n', length(infRows)))
    SNRinfo <- SNRinfo[-infRows, ]
  }

  if (!(season %in% c('year', '00'))) {
    if (season %in% c('summer','autumn','winter','spring')) {
      SNRinfo <- SNRinfo[SNRinfo$season == season, ]
    } else {
      SNRinfo <- SNRinfo[SNRinfo$month == season, ]
    }
  }

  SNRinfo
}


#' Convert a capture history table into an SNRinfo data.frame
#'
#' @description
#' \strong{Deprecated.} A thin wrapper around \code{\link{chtToSNRinfo}},
#' kept so the published Common Ground and Beyond Counting Calls analyses
#' keep reproducing exactly against current \code{main}. New code should
#' call \code{chtToSNRinfo} directly.
#'
#' @param snr A capture history data.frame containing detections and SNR info.
#'   Must have columns \code{detect_table1}, \code{detect_table2},
#'   \code{signalRMSdB}, \code{noiseRMSdB}, \code{t}, \code{month}, and
#'   \code{season}. Passing a file path here is deprecated and will error.
#'
#'   \code{detect_table1} is always read as ground truth: rows are filtered to
#'   \code{detect_table1 == 1} before anything else happens. For an
#'   observer-ground (OG) analysis this is the trusted observer's own column.
#'   For an adjudicated capture-recapture (CR) analysis, where neither raw
#'   observer is ground truth, \code{detect_table1} must be overwritten with
#'   the adjudicator's verdict before calling this function -- see
#'   \code{\link{cde}}'s \code{capHistTab} documentation for the same
#'   requirement and a pointer to a worked example.
#' @param season A timeCode corresponding to months, seasons, or 'year'
#'   (the default, which returns all rows).
#' @param groundTruthCol Column name treated as ground truth. Default
#'   \code{'detect_table1'}, matching original behaviour.
#' @param detectedCol Column name for the detector under evaluation.
#'   Default \code{'detect_table2'}, matching original behaviour.
#'
#' @returns SNRinfo data.frame with columns \code{Detected}, \code{CallRL},
#'   \code{NoiseRL}, \code{SNR}, \code{t}, \code{month}, and \code{season},
#'   filtered to the requested timeCode.
#' @export
capHist2snrInfo <- function(snr, season = 'year',
                            groundTruthCol = 'detect_table1',
                            detectedCol    = 'detect_table2'){
  .Deprecated('chtToSNRinfo')
  if (class(snr)=='character'){
  stop(paste0("Calling capHistToSnrInfo with a file name is deprecated.\n",
                 "Instead call capHist2SnrInfo with the capture ",
                 "history table in a data.frame.")
       )
  }
  timeCol <- if ('t0' %in% names(snr)) 't0' else 't'
  out <- chtToSNRinfo(snr,
                      groundTruth = sub('^detect_', '', groundTruthCol),
                      observers   = sub('^detect_', '', detectedCol),
                      signalCol = 'signalRMSdB', noiseCol = 'noiseRMSdB',
                      timeCol = timeCol, season = season)
  out[, c('Detected', 'CallRL', 'NoiseRL', 'SNR', 't', 'month', 'season')]
}


#' Convert capture history DATA.FRAME into the 'SNRinfo' format used
#' by the callDensity package
#'
#' @description
#' \strong{Deprecated.} A thin wrapper around \code{\link{chtToSNRinfo}},
#' kept so the published Common Ground and Beyond Counting Calls analyses
#' keep reproducing exactly against current \code{main}. New code should
#' call \code{chtToSNRinfo} directly.
#'
#' @param capHistTab - Capture history table of detections. Must have
#'   \code{detect_table1}/\code{detect_table2}, a \code{signalRMSdB1}/
#'   \code{signalRMSdB2} and \code{noiseRMSdB1}/\code{noiseRMSdB2} pair
#'   (averaged, matching original behaviour), and \code{datetime}.
#' @export
capHistTosnrInfo <- function(capHistTab){
  .Deprecated('chtToSNRinfo')
  out <- chtToSNRinfo(capHistTab,
                      groundTruth = 'table1',
                      observers   = 'table2',
                      signalCol = c('signalRMSdB1', 'signalRMSdB2'),
                      noiseCol  = c('noiseRMSdB1', 'noiseRMSdB2'),
                      timeCol   = 'datetime', season = 'year')
  out[, c('Detected', 'CallRL', 'NoiseRL', 'SNR', 't', 'month', 'season')]
}



#' Variance of a weighted mean following Chochran 1977 definition
#'
#' @param x list of data
#' @param w weights for data
#'
#' @export
var.wtd.mean.cochran <- function(x,w){
  ix <- !is.na(x)
  x <- x[ix]
  w <- w[ix]
  # Computes the variance of a weighted mean following Cochran 1977 definition
  n = length(w)
  xWbar = Hmisc::wtd.mean(x,w, na.rm = TRUE)
  wbar = mean(w, na.rm=TRUE)
  out = n/( (n-1)*sum(w)^2)*
    (sum((w*x-wbar*xWbar)^2)-2*xWbar*
       sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2)
    )
  return(out)
}

#' Fit an SNR-detection function
#'
#' @description
#' SNR detection functions are binomial models of the form detected ~ SNR. These
#' can be a GAM (default), GLM, or SCAM.
#'
#' @param SNRinfo A data.frame containing detection information in columns
#'     Detected and SNR
#'
#' @param modelType one of either 'GAM', 'GLM', or 'SCAM'. ModelType determines
#' which type of model will be used to estimate the detection-SNR function.
#'
#' @param numKnots Number of knots in the GAM or SCAM
#'
#' @importFrom stats binomial coef predict quasibinomial rnorm sd vcov
#'
#' @export
fitSNRdetectionFunc <- function(SNRinfo, modelType = "gam", numKnots = 3) {
  .Deprecated("fitDetFun")

  fitDetFun(
    SNRinfo = SNRinfo,
    modelType = modelType,
    numKnots = numKnots
  )
}

#' Fit an SNR-detection function with a closed population capture-recapture glm
#'
#' @description
#' VGLM SNR-detection functions are positive bernoulli models of the form
#' (y1,y2,...yN) ~ SNR. Here yi is a column of detections from the ith observer
#' (with a 1 for detected and 0 for not detected by that observer). The
#' assumptions for these models is that detection probability depends on
#' heterogeneity from SNR, observers are independent, and there are
#' no false positive detections included in the data for either observer.
#' Models are fitted using the VGAM package (Yee, Stoklosa & Huggins 2015).
#'
#' Yee, Thomas W., Jakub Stoklosa, and Richard M. Huggins. “The VGAM Package for
#'  Capture-Recapture Data Using the Conditional Likelihood.” Journal of
#'  Statistical Software 65 (June 1, 2015): 1–33.
#'  https://doi.org/10.18637/jss.v065.i05.
#'
#' @param SNRinfo A data.frame containing detection information in columns
#'     Detected and SNR
#'
#' @param yColNames a list of strings containing the column names for each
#' set of observer detections (default: 'detect_observer1','detect_observer2')
#'
#' @param whichObserver column name of the observer to use for predictions
#'
#'
#' @export
fitSNRvglm <- function(SNRinfo,
                       yColNames = c("detect_observer1","detect_observer2"),
                       whichObserver = "detect_observer2") {

  .Deprecated("fitDetFun")

  fitDetFun(
    SNRinfo = SNRinfo,
    modelType = "vglm",
    yColNames = yColNames,
    whichObserver = whichObserver
  )
}

#' Fit an SNR-detection function with a closed population capture-recapture gam
#'
#' @description
#' VGAM SNR-detection functions are positive bernoulli models of the form
#' (y1,y2,...yN) ~ SNR. Here yi is a column of detections from the ith observer
#' (with a 1 for detected and 0 for not detected by that observer). The
#' assumptions for these models is that detection probability depends on
#' heterogeneity from SNR, observers are independent, and there are
#' no false positive detections included in the data for either observer.
#' Models are fitted using the VGAM package (Yee, Stoklosa & Huggins 2015).
#'
#' Yee, Thomas W., Jakub Stoklosa, and Richard M. Huggins. “The VGAM Package for
#'  Capture-Recapture Data Using the Conditional Likelihood.” Journal of
#'  Statistical Software 65 (June 1, 2015): 1–33.
#'  https://doi.org/10.18637/jss.v065.i05.
#'
#' @param SNRinfo A data.frame containing detection information in columns
#'     Detected and SNR
#'
#' @param yColNames a list of strings containing the column names for each
#' set of observer detections (default: 'detect_observer1','detect_observer2')
#'
#' @param whichObserver column name of the observer to use for predictions
#'
#'
#' @export
fitSNRvgam <- function(...) {
  .Deprecated("fitDetFun (vglm)")
  stop("vgam support is experimental and not currently maintained.")
}


#' Model SNR-detection function with season as a factor (deprecated?)
#'
#' @param SNRinfo Data.frame containing columns SNR, Detected, and season
#'
#' @param season Factor containing a timeCode
#' @param useGLM Boolean, if true fit a GLM instead of GAM
#' @param numKnots Number of knots to use in the GAM
#' @importFrom grid unit
#'
#' @export
fitSNRbySeason <- function(SNRinfo, season=year, useGLM=TRUE, numKnots=3){
  # Initial glm:

  if (useGLM==TRUE){
    res.2=stats::glm(Detected~(SNR*season), data = SNRinfo,
              family = binomial(link="logit")) # or: family quasibinomial
  } else {
    res.2<-mgcv::gam(Detected~s(SNR,by=season, k=numKnots), data=SNRinfo,
               family=quasibinomial) # Test to include k=3
  }

  preds <- ggeffects::ggpredict(res.2,terms=c("SNR[all]","season[all]"),
                                pretty=FALSE, ci.lvl = 0.95)
  p2<-plot(preds, rawdata=TRUE, dot.alpha=0.01) + #facet_wrap(~1,nrow = 1)+
    ggplot2::geom_hline(yintercept = .50)+
    ggplot2::labs(x = c("SNR (dB)"), y = "Probability of detecting a call" ,
         title = NULL, color=NULL)+
    ggplot2::theme(legend.position="bottom",legend.key.height = unit(0.25,"cm"),
          legend.key.width = unit(0.25,"cm"),
          legend.text = element_text(size = 10))
  return(p2)
}

# Helper function to convert predicted response from vglms into a detection
# functions that we can use here
#' Predict probability of detection from a VGLM SNR-detection function
#'
#' Helper that combines the per-observer detection probability and the
#' "any-observer detected" probability from a fitted VGLM into a single
#' marginal probability of detection at each supplied SNR.
#'
#' @param snrDetFun.vglm A fitted VGLM SNR-detection function (e.g. from
#'   \code{\link{fitSNRvglm}}) carrying \code{@extra$whichObserver}.
#' @param snrs A data.frame with at least one column named \code{SNR}, giving
#'   the SNR values at which to predict the probability of detection.
#' @param whichObserver Character. Name of the observer column to use for
#'   predictions. Defaults to \code{snrDetFun.vglm@extra$whichObserver} (the
#'   observer baked into the fitted model by \code{\link{fitSNRvglm}}).
#'
#' @returns Numeric vector of probabilities of detection (one per row in
#'   \code{snrs}) for the chosen observer.
#' @export
predVglmPDet <- function(...){
    .Deprecated("predictDetFun")
    predictDetFun(...)$fit

}



#' Load specified density results text files in a directory
#'
#' @returns Data.frame of all call density results in the specified path
#'
#' @param path Location (path) of text files containing call density results
#'
#' @param densityString Prefix of _density files to be loaded
#'
#' @export
collectDensityResults <- function (path,
                                   densityString=paste0('^density_',siteCode)){
  files <- dir(path=path, pattern = densityString, full.names = TRUE)
  d <- utils::read.csv(files[1])

  for (i in 2:length(files))   {
    d <- rbind(d,read.csv(files[i]))
  }
  d$Dc <- d$Dc*1e3
  d$season <- factor(d$season, levels=c('1','2','3','4','5','6',
                                        '7','8','9','10','11','12',
                                        'summer','autumn','winter','spring','year') )
  d<- with(d, d[order(season),])

  # Trying to get 95% confidence intervals here, but not very confident
  d$CI.low <- (d$Dc-(1.96*d$Dc*d$CV.Dc))
  d$CI.high <- (d$Dc+(1.96*d$Dc*d$CV.Dc))
  # d <- format(d,digits=3,nsmall=0, scientific=T)
  return(d)
}

#' Show table of inputs into call density estimate
#'
#' @param d Data.frame containing call density results
#'
#' @param separateCommon Boolean - if true input parameters k, w, SLmean, SLsd
#'     will be put into a separate table with only a single row.
#'
#' @export
densityInputTable <-function(d, separateCommon=TRUE){

  if (separateCommon){
    commonTable <- subset(d[1,], select=c('k','A','SLmean','SLsd'))
    rownames(commonTable) <- 'All'
    common <- kableExtra::kbl(commonTable, row.names=TRUE) %>%
      kableExtra::kable_classic_2(full_width=FALSE)

    inputTable <- subset(d,select=c('season','Nc','c','T','pa',
                                    'NLmean','NLsd') )
    inputs <- kableExtra::kbl( inputTable, row.names = FALSE,
                   digits=c(0,0,3,1,4,1,1) ) %>%
      kableExtra::kable_classic_2(full_width=FALSE)

    return(list(inputTable,inputs,commonTable,common))
  } else {
    inputTable <- subset(d,select=c('season','Nc','c','k','T','w','pa',
                                    'NLmean','NLsd','SLmean','SLsd') )

    inputs <- kableExtra::kbl( inputTable, row.names = FALSE,
                   digits=c(0,0,3,0,1,0,4,1,1,1,1) ) %>%
      kableExtra::kable_classic_2(full_width=FALSE)
  }
  return(list(inputTable,inputs))
}

#' Table of call densities including CVs (ggplot2)
#' @param d Data.frame of call density results#'
#' @export
densityResultsTable<- function (d){

  resultsTable <- subset(d, select=c('season','Dc',
                                     'CV.Nc','CV.c','CV.pa','CV.Dc') )

  result<- kableExtra::kbl( resultsTable, row.names = FALSE,
                            digits=c(0,3,4,3,3,3) ) %>%
    kableExtra::kable_classic_2(full_width = FALSE)
  return (list(resultsTable,result))
}


#' Calculate study area for call density estimate.
#'
#' Estimate study area for call density given a radius from the hydrophone. If
#' only a radius is provided, then the area=pi*radius^2 is returned. Optionally,
#' if a matrix of truncation distances are provided, then the study area is
#' calculated as the sum of the areas of each radial transect (assuming N
#' transects will have the same 360/N degrees of angular coverage around the
#' circle).
#'
#' @param w  - Radius of study area
#' @param truncationDistance - Matrix of truncation distances \code{[1 x Nr]}. Nr is
#'   the number of radial transects in the simulation (i.e. the same as the
#'   number of TL profiles).
#'
#' @returns Area of the simulation
#' @export
studyArea <- function(w, truncationDistance=w){
  A = pi*w^2
  if (any(truncationDistance < w)     ){
    if (length(truncationDistance) == 1) {
      A = pi*truncationDistance^2
    } else {
      A = sum(pi*truncationDistance^2/length(truncationDistance))
    }
  }
  return(A)
}


#' Confidence interval from density and it's CV
#'
#' Quick function to convert an abundance with a coefficient of variation (CV)
#' into a 95% confidence interval.
#'
#' @param a abundance (or density)
#' @param cv coefficient of variation, in decimal format, i.e., a CV of 10%
#'   would be 0.1.
#'
#' @returns list with lower and upper 95% confidence intervals
#' @export
ciFromCV <- function(a, cv){
  var.a<- (a*cv)^2

  var.log.a<- log( 1+ (var.a/(a^2)))

  C<- exp( 1.96 * sqrt(var.log.a))

  lower.95CI<- a/C
  upper.95CI<- a*C

  ans<- list( lower.CI95 = lower.95CI, upper.CI95 = upper.95CI )
  return(ans)
}

#' Read a capture history csv file (e.g. created in Matlab)
#'
#' Required columns for each observer are detect_, snr_, signalRMSdB_,
#' noiseRMSdB_, and the underscore is followed by input parameter observerNames.
#'    -detect_ columns are binary i.e. 1 a detection and 0 for absence
#'    -snr_, signalRMSdB_, and noiseRMSdB_ are all numeric and in dB
#'    -groundTruth is the binary column to use as ground-truth detections. It
#'     can be a detect_observerName or another column.
#' OG capture history tables have exactly two observers with suffixes,
#'   observerNames=c('table1','table2'). These are the type of tables used in
#'   the 'Beyond Counting Calls' manuscript.
#' CR capture history tables can have N observers and can have any number of
#'   suffixes of any form. MultiCaptureHistoryTable from the common ground
#'   manuscript produces suffixes of the form:
#'   ('observer1','observerN-1','observerN')
#'
#' @param capHistFile - name of the csv or txt file that contains capture
#'   histories. TODO: Document the required columns/file format better
#' @param observerNames - list of strings that correspond to the suffixes used
#'   for each observer. Default: c('table1','table2')
#' @param groundTruth - Name of column with binary data used as ground-truth.
#' @param whichObserver - Name of observer for which detection function will be
#'   modelled, and call density estimated
#' @returns capture history table as data.frame
#' @export
readCapHist <- function(capHistFile,
                        observerNames = c('table1','table2'),
                        groundTruth = 'table1',
                        whichObserver = 'detect_table2'){
  ch <- read.csv(file = capHistFile,na.strings = c('NaN','NA'))

  # TODO: Check for required input columns

  # TODO: Convert columns into format that can be used downstream

  # Add columns for times and seasons
  ch <- capHistTimeSeason(ch)
  return(ch)
}

# Convert multi-observer capture history table into an CR capture history table
#' Convert a multi-observer capture history table to a two-observer CR table
#'
#' Extracts a two-observer (table1, table2) capture-recapture capture history
#' table from a multi-observer source table by selecting only rows where at
#' least one of the two named observers detected the event, and renaming the
#' relevant columns to the \code{_table1}/\code{_table2} suffixes expected by
#' the rest of the callDensity package. The \code{verdict} column from the
#' source is used as the ground-truth detection for table1.
#'
#' @param d A multi-observer capture history data.frame containing per-event
#'   columns \code{t0}, \code{tEnd}, \code{fLow}, \code{fHigh}, \code{SNR},
#'   \code{signalRMSdB}, \code{noiseRMSdB}, \code{noiseDev}, \code{verdict},
#'   plus per-observer detection columns named
#'   \code{detect_<observerSuffix>}.
#' @param table1suffix Character. Suffix of the observer columns to remap to
#'   the \code{_table1} role (treated as ground truth via \code{verdict}).
#' @param table2suffix Character. Suffix of the observer columns to remap to
#'   the \code{_table2} role (the detector being evaluated).
#'
#' @returns A capture history data.frame in the two-observer (table1, table2)
#'   format expected by \code{\link{cde}}, with time/season columns added by
#'   \code{\link{capHistTimeSeason}}.
#'
#' @section Deprecated:
#' No current analysis needs this reduction. \code{cde()}'s only use of
#' \code{capHistTab} is \code{falseDiscoveryRate()}, which already accepts
#' \code{gtColName}/\code{testColName} directly against native
#' \code{detect_observerN} columns (once \code{cde()} forwards them --
#' see its own argument list). \code{fitDetFun(modelType = 'vglm')}
#' accepts any number of native-named observer columns directly via
#' \code{yColNames}, with no reduction to two required at all --
#' \code{VGAM::posbernoulli.t} is not limited to two occasions. This
#' function's body is unchanged and will keep working, kept only so the
#' published Common Ground and Beyond Counting Calls analysis scripts
#' keep reproducing exactly. New work should skip the 2-observer
#' reduction entirely and use \code{\link{chtToSNRinfo}} with a
#' \code{observers} vector of whatever length is needed.
#' @export
mchToCR <- function (d, table1suffix, table2suffix){
  .Deprecated(msg = paste(
    "mchToCR() is deprecated: cde()/fitDetFun(modelType='vglm') no longer",
    "need the 2-observer table1/table2 reduction. Use chtToSNRinfo()",
    "directly on the native multi-observer table instead -- see",
    "?mchToCR's Deprecated section."))
  # Column names
  d1 <- paste0('detect_',table1suffix)      # detection from first observer
  d2 <- paste0('detect_',table2suffix)      # detection from detector of interest
  ch <- subset(d, subset = (d[,d1] | d[,d2]),
               select = c('t0','tEnd','fLow','fHigh','SNR','signalRMSdB',
                          'noiseRMSdB','noiseDev','verdict',d1,d2))
  ch$t <- ch$t0
  names(ch)[1:4]<- paste0(names(ch)[1:4],'_table1')
  names(ch)[names(ch)=='verdict']<-'detect_table1'
  # names(ch) <- gsub(table1suffix,'table1',names(ch))
  names(ch) <- gsub(table2suffix,'table2',names(ch))
  ch<-capHistTimeSeason(ch)
  return(ch)
}



#' Resolve column names from a prefix pattern or explicit vector
#'
#' Helper to standardise column selection across functions. If \code{prefix} is
#' a single string it is used as a \code{grepl} pattern against \code{df_names};
#' if it is a character vector of length > 1 it is used directly as column names.
#'
#' @param df_names Character vector of column names (i.e. \code{names(df)}).
#' @param prefix Either a single character string pattern or a character vector
#'   of explicit column names.
#'
#' @return A character vector of resolved column names.
#'
#' @keywords internal
resolveColumns <- function(df_names, prefix) {
  if (length(prefix) > 1) {
    prefix
  } else {
    df_names[grepl(prefix, df_names)]
  }
}


#' Pivot SNR observer columns to long format
#'
#' Helper used by the SNR plotting functions. Selects columns matching
#' \code{snr_prefix} and pivots them to long format, attaching observer labels.
#'
#' @param df A data frame containing SNR observer columns.
#' @param snr_prefix Either a single character string pattern identifying SNR
#'   columns (e.g. \code{"snr_observer"}), or a character vector of explicit
#'   column names. Defaults to \code{"snr_observer"}.
#' @param obs_labels Named character vector mapping numeric observer indices
#'   (as character) to display labels, e.g. \code{c("1"="Alice","2"="Bob")}.
#' @param extra_cols Integer or character vector of additional column indices or
#'   names to retain (e.g. \code{"i"} for a time index). Defaults to
#'   \code{NULL}.
#'
#' @return A long-format data frame with columns \code{observer} (factor) and
#'   \code{snr}, plus any columns specified in \code{extra_cols}.
#'
#' @importFrom tidyr pivot_longer
#' @keywords internal
pivotSNR <- function(df, snr_prefix = "snr_observer", obs_labels,
                     extra_cols = NULL) {
  snr_cols <- resolveColumns(names(df), snr_prefix)
  nObs     <- length(snr_cols)
  col_ix   <- c(snr_cols, extra_cols)

  # strip the common prefix for names_prefix only when pattern-based
  names_pfx <- if (length(snr_prefix) == 1) snr_prefix else ""

  out <- tidyr::pivot_longer(df[, col_ix, drop = FALSE],
                             cols           = tidyr::all_of(snr_cols),
                             names_prefix   = names_pfx,
                             names_to       = "observer",
                             values_to      = "snr",
                             values_drop_na = TRUE)
  out$observer <- obs_labels[out$observer]
  out$observer <- factor(out$observer, levels = obs_labels)
  out
}

#' Compute pretty SNR limits across observer columns
#'
#' Finds the minimum and maximum SNR values across all matching observer
#' columns and returns pretty axis limits.
#'
#' @param df A data frame containing SNR observer columns.
#' @param snr_prefix Either a single character string pattern identifying SNR
#'   columns, or a character vector of explicit column names. Defaults to
#'   \code{"snr_observer"}.
#'
#' @return A numeric vector of length 2 giving \code{c(min, max)} pretty limits.
#'
#' @examples
#' \dontrun{
#' computeSNRLims(ap)
#' computeSNRLims(d, snr_prefix = c("snr_observer1", "snr_observer2"))
#' }
#'
#' @export
computeSNRLims <- function(df, snr_prefix = "snr_observer") {
  snr_cols <- resolveColumns(names(df), snr_prefix)
  all_vals <- unlist(df[, snr_cols, drop = FALSE], use.names = FALSE)
  lims     <- pretty(range(all_vals, na.rm = TRUE))
  c(min(lims), max(lims))
}

#' Compute row-wise mean across observer columns
#'
#' For each entry in \code{prefixes}, finds matching columns and appends a new
#' column to \code{d} containing the row-wise mean, ignoring \code{NA}s.
#'
#' @param d A data frame containing observer metric columns.
#' @param prefixes A named list or named character vector where names become the
#'   new column names and values are either a single string pattern or a
#'   character vector of explicit column names. Defaults to
#'   \code{list(SNR = "snr_observer", noiseRMSdB = "noiseRMSdB_observer",
#'   signalRMSdB = "signalRMSdB_observer")}.
#'
#' @return The input data frame \code{d} with one additional column per entry
#'   in \code{prefixes}.
#'
#' @examples
#' \dontrun{
#' # Default: pattern-based
#' d <- addObserverMeans(d)
#'
#' # Mixed: pattern-based and explicit
#' d <- addObserverMeans(d, prefixes = list(
#'   SNR        = "snr_observer",
#'   noiseRMSdB = c("noiseRMSdB_observer1", "noiseRMSdB_observer2")
#' ))
#' }
#'
#' @export
addObserverMeans <- function(d,
                             prefixes = list(
                               SNR         = "snr_observer",
                               noiseRMSdB  = "noiseRMSdB_observer",
                               signalRMSdB = "signalRMSdB_observer"
                             )) {
  for (col_name in names(prefixes)) {
    matched    <- resolveColumns(names(d), prefixes[[col_name]])
    d[[col_name]] <- rowMeans(d[, matched, drop = FALSE], na.rm = TRUE)
  }
  d
}

#' Plot an SNR detection function with mirrored SNR distributions
#'
#' Shows the fitted detection function together with the distribution of
#' detected and missed true positives. Detected calls are plotted above the
#' detection curve and missed calls below.
#'
#' @param model A fitted detection model returned by fitSNRdetectionFunc(),
#'   fitSNRvglm(), or fitSNRvgam().
#' @param SNRinfo Original data.frame used to fit the model.
#' @param whichObserver Observer column for vglm/vgam models. Defaults to the
#'   observer stored in model@extra$whichObserver.
#' @param distribution One of "density", "histogram", or "none".
#' @param mirror.height Height of the mirrored distributions.
#' @param npoints Number of SNR values for prediction.
#'
#' @return A ggplot object.
#'
#' @export
showSNRdetectionFunc <- function(...) {
  .Deprecated("showDetFun")
  showDetFun(...)
}



#' Bar-plot of detection rates by timeCodes (ggplot2)
#' @param d Data.frame of call density results
#'
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar ylab xlab labs aes
#'   theme_minimal theme element_text
#' @export
detectionRatePlot <- function(d){

  gRate <- ggplot2::ggplot(data=d, aes(x=season, y=Nc/T, fill=season))+
    ggplot2::geom_bar(stat="identity",show.legend=F)+
    # scale_y_continuous(limits=c(0,15))+
    # geom_errorbar(aes(ymin=CI.low, ymax=CI.high),
    #      width=0.5, position=position_dodge(0.9))+
    ggplot2::ylab(latex2exp::TeX(
      "Mean detection rate ($ calls \\cdot h^{-1})" ) )+
    ggplot2::xlab('')+
    ggplot2::ggtitle(siteCode, )+
    ggplot2::labs(tag="(A)") +
    ggplot2::theme_minimal()+
    ggplot2::theme(plot.tag = element_text(),
                   plot.tag.position = c(0.2, 0.975),
                   plot.title = element_text(hjust = 0.5))
}

#' Bar-plot of detection rates by timeCodes and scaled (corrected) by precision (ggplot2)
#' @param d Data.frame of call density results
#'
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar ylab xlab labs aes
#'   theme_minimal theme element_text
#' @export
detectionRateCorrectedPlot <- function(d){

  gRateCorrected <- ggplot2::ggplot(data=d, aes(x=season, y=(Nc*c)/T, fill=season))+
    ggplot2::geom_bar(stat="identity",show.legend=F)+
    # scale_y_continuous(limits=c(0,15))+
    # geom_errorbar(aes(ymin=CI.low, ymax=CI.high),
    #               width=0.5, position=position_dodge(0.9))+
    ggplot2::ylab(latex2exp::TeX(
      "Mean detection rate ($ calls \\cdot h^{-1})" ) )+
    ggplot2::xlab('')+
    ggplot2::ggtitle(siteCode, )+
    ggplot2::labs(tag="(A)") +
    ggplot2::theme_minimal()+
    ggplot2::theme(plot.tag = element_text(),
                   plot.tag.position = c(0.2, 0.975),
                   plot.title = element_text(hjust = 0.5))
}

#' Bar-plot of call density by timeCodes (ggplot2)
#' @param d Data.frame of call density results
#'
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar ylab xlab labs aes
#'   theme_minimal theme element_text
#' @importFrom latex2exp TeX
#' @export
densityPlot <- function(d){

  gDen <- ggplot2::ggplot(data=d, aes(x=season, y=Dc, fill=season))+
    ggplot2::geom_bar(stat="identity",show.legend=F)+
    ggplot2::geom_errorbar(aes(ymin=CI.low, ymax=CI.high), width=0.5,)+
    # scale_y_continuous(limits=c(0,0.25))+
    ggplot2::ylab(latex2exp::TeX(
      "Call density ($ calls \\cdot h^{-1} \\cdot 1000 \\cdot km^{-2}$)" ) )+
    ggplot2::xlab('')+
    ggplot2::labs(tag="(B)") +
    ggplot2::theme_minimal()+
    ggplot2::theme(plot.tag = element_text(),
                   plot.tag.position = c(0.225, 0.975))
}



#' Multi-Observer Detection Count
#'
#' Counts and visualises detection agreement patterns across multiple observers.
#' Each row is encoded as a binary string representing each observer's detection
#' decision, and the frequency of each unique pattern is returned as a data frame
#' and displayed as a bar chart.
#'
#' @param ch A data frame containing observer detection columns. Each detection
#'   column should contain binary values (0/1).
#' @param detect_prefix Either a single character string used as a pattern to
#'   identify detection columns by name (e.g. \code{"detect_observer"} will
#'   match \code{detect_observer1}, \code{detect_observer2}, etc.), or a
#'   character vector of explicit column names (e.g.
#'   \code{c("detect_obs1", "detect_obs2")}). Defaults to
#'   \code{"detect_observer"}.
#' @param ... Additional ggplot2 layers (e.g. scales, themes, annotations)
#'   passed to the plot via \code{+}.
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{category}{A factor giving the binary string pattern of observer
#'       detections (e.g. \code{"101"} means observers 1 and 3 detected,
#'       observer 2 did not).}
#'     \item{count}{An integer giving the number of rows matching that pattern.}
#'   }
#'
#' @details
#' When \code{detect_prefix} is a single string, columns are identified by
#' \code{grepl(detect_prefix, names(ch))}. When \code{detect_prefix} is a
#' character vector of length > 1, those column names are used directly.
#'
#' @examples
#' ch <- data.frame(
#'   detect_observer1 = c(1, 0, 1, 1),
#'   detect_observer2 = c(1, 1, 0, 1),
#'   detect_observer3 = c(0, 0, 1, 1)
#' )
#'
#' # Default
#' result <- multiObserverDetectionCount(ch)
#'
#' # With extra ggplot2 layers
#' result <- multiObserverDetectionCount(ch,
#'   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)))
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_bw theme element_line
#' @export
multiObserverDetectionCount <- function(ch,
                                        detect_prefix = "detect_observer",
                                        ...) {
  detectColumns <- resolveColumns(names(ch), detect_prefix)
  nObservers    <- length(detectColumns)

  detBin <- apply(ch[, detectColumns, drop = FALSE], 1, function(row) {
    paste(row, collapse = "")
  })

  ven       <- factor(detBin)
  counts_df <- as.data.frame(table(ven))
  colnames(counts_df) <- c("category", "count")

  p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = category, y = count)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      x = sprintf("Observer binary mask (%s)", paste(1:nObservers, collapse = "")),
      y = "Number of detections"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey80"),
                   panel.grid.minor = ggplot2::element_line(colour = "grey90")) +
    Reduce(`+`, list(...))

  print(p)
  return(counts_df)
}


#' Plot SNR time series for the adjudicated subset
#'
#' Pivots SNR observer columns to long format and produces a scatter plot of
#' SNR against a time index column, coloured by observer.
#'
#' @param ap A data frame containing SNR observer columns and a time index column.
#' @param obs_labels Named character vector mapping numeric observer indices
#'   (as character) to display labels.
#' @param snr_prefix Either a single character string pattern identifying SNR
#'   columns, or a character vector of explicit column names. Defaults to
#'   \code{"snr_observer"}.
#' @param time_col Name of the time index column. Defaults to \code{"i"}.
#' @param snr_lims Numeric vector of length 2 giving y-axis limits. If
#'   \code{NULL} (default), limits are computed from the data using
#'   \code{\link{computeSNRLims}}.
#' @param title Plot title. Defaults to
#'   \code{"Adjudicated subset: SNR time series by detector"}.
#' @param ... Additional ggplot2 layers (e.g. scales, themes, annotations)
#'   passed to the plot via \code{+}.
#'
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' \dontrun{
#' plotSNRTimeSeries(ap, obs_labels)
#'
#' # Shared limits with extra layer
#' lims <- computeSNRLims(d)
#' plotSNRTimeSeries(ap, obs_labels, snr_lims = lims,
#'   ggplot2::theme(legend.position = "bottom"))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point ylim labs theme_bw
#' @export
plotSNRTimeSeries <- function(ap,
                              obs_labels,
                              snr_prefix = "snr_observer",
                              time_col   = "i",
                              snr_lims   = NULL,
                              title      = "Adjudicated subset: SNR time series by detector",
                              ...) {
  snr_cols <- resolveColumns(names(ap), snr_prefix)
  snrs     <- pivotSNR(ap, snr_prefix, obs_labels, extra_cols = time_col)
  snr_lims <- snr_lims %||% computeSNRLims(ap, snr_prefix)

  p <- ggplot2::ggplot(snrs, ggplot2::aes(x = .data[[time_col]],
                                          y = snr,
                                          colour = observer)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::ylim(snr_lims) +
    ggplot2::labs(title = title, x = "Date", y = "SNR (dB)") +
    ggplot2::theme_bw() +
    Reduce(`+`, list(...))

  print(p)
  invisible(p)
}


#' Plot SNR histogram by observer
#'
#' Pivots SNR observer columns to long format and produces a faceted histogram
#' of SNR distributions, one panel per observer.
#'
#' @param df A data frame containing SNR observer columns.
#' @param obs_labels Named character vector mapping numeric observer indices
#'   (as character) to display labels.
#' @param snr_prefix Either a single character string pattern identifying SNR
#'   columns, or a character vector of explicit column names. Defaults to
#'   \code{"snr_observer"}.
#' @param snr_lims Numeric vector of length 2 giving x-axis limits. If
#'   \code{NULL} (default), limits are computed from the data using
#'   \code{\link{computeSNRLims}}.
#' @param binwidth Histogram bin width in dB. Defaults to \code{0.5}.
#' @param title Plot title. Defaults to \code{"SNR distribution by detector"}.
#' @param ... Additional ggplot2 layers (e.g. scales, themes, annotations)
#'   passed to the plot via \code{+}.
#'
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' \dontrun{
#' plotSNRHistogram(ap, obs_labels,
#'   title = "Adjudicated subset: SNR distribution by detector")
#'
#' # Shared limits with extra layer
#' lims <- computeSNRLims(d)
#' plotSNRHistogram(ap, obs_labels, snr_lims = lims,
#'   title = "Adjudicated subset: SNR distribution by detector",
#'   ggplot2::scale_fill_brewer(palette = "Dark2"))
#' plotSNRHistogram(d, obs_labels, snr_lims = lims,
#'   title = "Full dataset: SNR distribution by detector",
#'   ggplot2::scale_fill_brewer(palette = "Dark2"))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid xlim labs theme_bw
#' @export
plotSNRHistogram <- function(df,
                             obs_labels,
                             snr_prefix = "snr_observer",
                             snr_lims   = NULL,
                             binwidth   = 0.5,
                             title      = "SNR distribution by detector",
                             ...) {
  snr_cols <- resolveColumns(names(df), snr_prefix)
  snrLng   <- pivotSNR(df, snr_prefix, obs_labels)
  snr_lims <- snr_lims %||% computeSNRLims(df, snr_prefix)

  p <- ggplot2::ggplot(snrLng, ggplot2::aes(snr, fill = observer)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.4,
                            binwidth = binwidth) +
    ggplot2::facet_grid(observer ~ .) +
    ggplot2::xlim(snr_lims) +
    ggplot2::labs(title = title, x = "SNR (dB)", y = "Count") +
    ggplot2::theme_bw() +
    Reduce(`+`, list(...))

  print(p)
  invisible(p)
}

#' Plot probability of detection as a directional range footprint
#'
#' @description
#' Renders p(det) as a function of range and azimuth as a polar heatmap
#' centred on the recorder -- ports Brian's MATLAB \code{plotPDetRadials.m}
#' to a native ggplot2 \code{coord_polar()} plot, rather than manually
#' unrolling to Cartesian coordinates (\code{pcolor(x = r*cos(theta),
#' y = r*sin(theta), z)}) the way MATLAB does. ggplot2 has a polar
#' coordinate system built in, so that workaround isn't needed here.
#'
#' @param pDetResults Either the full list returned by
#'   \code{\link{pDetInArea}} (in which case \code{$allDetFunctions} is
#'   used), or that data.frame directly -- \code{range_m} in the first
#'   column, one column per transect named \code{tl<angle>} (matching
#'   \code{\link{simTLradials_20logR}}'s convention) giving p(det) at
#'   every range along that transect.
#' @param maxRange Optional upper limit (m) for the range axis. Default
#'   \code{NULL} uses the full range present in the data.
#'
#' @details
#' Assumes a constant range step (true of \code{simTLradials_20logR}'s
#' output and of \code{pDetInArea}'s own \code{output.resolution.m}
#' spacing) so each tile can use a fixed height. Irregular range bins
#' (e.g. the log-spaced bins in \code{plotSnrRadials.m}) would need
#' explicit per-row bin edges via \code{geom_rect} instead -- not
#' implemented here, since nothing in the package currently produces that
#' shape.
#'
#' @returns A ggplot object. Azimuth is plotted compass-style (0 = North,
#'   increasing clockwise), matching \code{simTLradials_20logR}'s
#'   convention -- \code{coord_polar}'s \code{start} is defined as an
#'   offset from 12 o'clock, and \code{direction = 1} is clockwise, so
#'   \code{start = 0, direction = 1} maps azimuth directly onto compass
#'   bearing with no rotation needed. Verified by rendering a synthetic
#'   case (known high p(det) at due north, medium at due east) and
#'   visually confirming the bright/grey wedges land where expected,
#'   rather than assumed from the parameter names alone.
#' @importFrom ggplot2 ggplot aes geom_tile coord_polar scale_x_continuous
#'   scale_fill_gradient labs theme_minimal
#' @export
plotPDetRadials <- function(pDetResults, maxRange = NULL) {

  pdet <- if (!is.null(pDetResults$allDetFunctions)) {
    pDetResults$allDetFunctions
  } else {
    pDetResults
  }

  transectCols <- setdiff(names(pdet), "range_m")
  if (length(transectCols) == 0) {
    stop("plotPDetRadials: no transect columns found (expected tl<angle> ",
         "columns alongside range_m).")
  }

  az <- suppressWarnings(as.numeric(sub("^tl", "", transectCols)))
  if (any(is.na(az))) {
    stop("plotPDetRadials: could not parse azimuth from column names ",
         "(expected 'tl<angle>', e.g. 'tl0', 'tl90'). Got: ",
         paste(transectCols, collapse = ", "))
  }

  long <- data.frame(
    range_km = rep(pdet$range_m / 1e3, times = length(transectCols)),
    azimuth  = rep(az, each = nrow(pdet)),
    pDet     = unlist(pdet[transectCols], use.names = FALSE)
  )

  # Wraparound: duplicate the smallest-azimuth transect's data at
  # az + 360, so the polar tile grid closes seamlessly rather than
  # leaving a wedge-shaped gap between the last and first transect --
  # the same fix MATLAB's manual az=360 duplicate column applies.
  az0  <- min(az)
  wrap <- long[long$azimuth == az0, ]
  wrap$azimuth <- az0 + 360
  long <- rbind(long, wrap)

  if (!is.null(maxRange)) {
    long <- long[long$range_km <= maxRange / 1e3, ]
  }

  angleStep <- if (length(az) > 1) min(diff(sort(unique(az)))) else 360
  rangeStep <- if (nrow(pdet) > 1) {
    min(diff(sort(unique(pdet$range_m)))) / 1e3
  } else {
    max(pdet$range_m) / 1e3
  }

  ggplot2::ggplot(long, ggplot2::aes(x = azimuth, y = range_km, fill = pDet)) +
    ggplot2::geom_tile(width = angleStep, height = rangeStep) +
    ggplot2::coord_polar(theta = "x", start = 0, direction = 1) +
    ggplot2::scale_x_continuous(breaks = seq(0, 270, 90),
                                labels = c("N", "E", "S", "W")) +
    ggplot2::scale_fill_gradient(low = "black", high = "white",
                                 limits = c(0, 1), name = "P(det)") +
    ggplot2::labs(x = NULL, y = "Range (km)") +
    ggplot2::theme_minimal()
}
