#' Call Density Estimate - Main function that accepts parameter file and season#' title: "Call Density"
#' author: "Brian Miller"
#' date: "2022-10-12"
#'
#' Density of blue whale D-calls (BmD) in long-term moored recording dataset.
#' Apply methods of Castro, Harris, et al (in prep) to obtain call density.
#'
#' $D_c = \\frac{N_c*(1-c)}{kTP_a{\\pi}w^2}$
#
#' where:
#' $D_c$ is call density
#' $N_c$ is number of calls
#' $c$ is false discovery rate
#' $k$ is number of sensors (here always 1)
#' $T$ is the duration of data analysed
#' $P_a$ is the probability of detection in the study area
#' ${\\pi}w^2$ is the study area (in km\\^2)
#' @param p A data.frame containing parameters for call density estimation.
#'    Usually created by calling function defaultOutputFileNames
#' @param season TimeCode specifying month, season, or year for outputs
#' @param snrDetFun OPTIONAL linear-model like structure (GLM,GAM,SCAM,etc)
#'    specifying the SNR-detection function to use. If this is not included,
#'    then the SNR-detection function will be derived from the capture history
#'    table.
#' @param truncationDistance scalar or matrix of truncation distances.
#'    If a matrix is provided, then the dimensions should be 1xN with N being
#'    the same as the number of transects
#' @param snrTruncationThreshold scalar SNR in dB below which the probability of
#'    detection will be forcibly set to zero.
#' @param NL data.frame containing distribution of noise level parameters. This
#'   data.frame must contain the rows mean, sd, and sampleSize (similar to SL).
#'
#' @returns data.frame containing call density inputs, results, and CVs
#'
#' @importFrom stats rnorm sd vcov
#' @importFrom utils read.csv write.table
#' @importFrom magrittr %>%
#' @export
#'
cdeFromParamFile <- function (p,season, snrDetFun=NULL, truncationDistance=Inf,
                 snrTruncationThreshold=-Inf, NL=NULL){
  # Check inputs and create outputs
  # Store all parameters for call-density estimation in data frame called 'p'
  #
  # InputCheck------------------------------------------------------------------

  # Require that parameters have already been defined and are in data.frame, p
  if (!exists("p")){
    warning(c(
      'Required parameter variable, p, not found.
  Please define all call density parameters in data.frame, p, before running
  this script.\n'))
  }

  if (!exists("season")){
    warning(c(
      'Required parameter variable, season, not found.
  Please define a variable called season, and assign it a value of
  summer, autumn, winter, spring, or year.\n\n
  Now assuming season<-year'))
    season <-'year'
  }

  # Number of calls
  # Nc--------------------------------------------------------------------------
  Nc <- countDetections(p$detectorParams$fullYearDetectionCsv,season,
                        snrTruncationThreshold=snrTruncationThreshold)

  # Duration of time monitored
  #T----------------------------------------------------------------------------
  T <- deploymentDurationFromsoundFolderCsv(p$detectorParams$fullYearEffortFile, season)

  #  ### Calculate study area
  # A---------------------------------------------------------------------------
  A = studyArea(p$w,truncationDistance)

  # c, CV_c---------------------------------------------------------------------
  ch <- readCapHist(p$capHistFile)
  fdr <- falseDiscoveryRate(ch, season, snrTruncationThreshold)

  # CV.c and c are packed in a list, so unpack
  CV.c <- fdr$cv.c
  c <- fdr$c


  #   ### $p_a$ (Overall probability of detection)
  # This is the most complicated part and takes the longest time
  # Pa--------------------------------------------------------------------------

  SNRinfo <- capHist2snrInfo(ch,season)

  # Read SL distribution parameters from overall parameters.
  SL <- data.frame(mean=p$slParams$slMean,
                   sd=p$slParams$slStd,
                   sampleSize=p$slParams$slSampleSize)

  TL<-utils::read.csv(p$tlParams$tlFile)

  # Check whether user supplied an snr detection function or estimate from data
  # NB: This needs to occur AFTER any SNR truncation.
  if (is.null(snrDetFun)){
    # Estimate snrDetFun from SNRinfo
    snrDetFun <- fitSNRdetectionFunc(
      subset(SNRinfo,SNR>=snrTruncationThreshold),
      modelType=p$modelType, p$numKnots)
  }

  # If user has not specified NL distribution (data.frame with columns mean, sd,
  # samplesize), then extract this information it from SNRinfo.
  # NB: This needs to occur with untruncated SNRInfo (i.e. prior to truncation).
  if (is.null(NL)){
    NL <- nlFromSnrInfo(SNRinfo, snrDetFun)
  }

  pDetResults <- pDetInArea(snrDetFun, SL, TL,  NL, # Sonar equation inputs
                            output.resolution.m=p$output.resolution.m,
                            outerloop=p$outerloop,
                            truncationDistance=truncationDistance,
                            snrTruncationThreshold= snrTruncationThreshold,
                            p$transectFile,  # file output names
                            p$simResultsFile,
                            p$paFile)

  # The above function, pDetInArea writes results to a bunch of files
  # Load the file that we need that has the p_a in them, and ignore the others
  # pa.all.transects <- read.csv(p$paFile, header = TRUE, sep = ' ');
  pa.all.transects <- pDetResults$perTransectMeanSD
  no.transects <- dim(pa.all.transects)[1]-1;
  pa <- pa.all.transects[no.transects+1,1]

  #CV_Pa------------------------------------------------------------------------
  # The line below will read the per-transect pa from file.
  # pa.all.transects <- read.csv(p$paFile, header = TRUE, sep = ' ');

  # Per-transect mean and SD of pa now returned in the list of pDetResults, so
  # no need to read or write it to a file
  pa.all.transects <- pDetResults$perTransectMeanSD
  CV.pa <- pa_CV(pa.all.transects)

  #CV_Nc------------------------------------------------------------------------
  CV.Nc <- Nc_CV(Nc,pa,c)

  #  ### Finally, estimate Dc (Density of calls) [calls $h^{-1} km^{-2}$]
  # Dc --------------------------------------------------------------------------
  Dc <- (Nc * (1-c) )/( p$k * A * pa * T )
  Dc     # per h per km^2
  Dc*1e3 # per h per 1000 km^2?

  # CV_Total--------------------------------------------------------------------
  CV.Dc <- Dc_CV(CV.Nc,CV.pa,CV.c)

  #  ## Collate into data frame and write to csv file
  # Format output---------------------------------------------------------------

  result <- data.frame(season, p$siteCode,Nc,c,p$k,T,p$w,pa,
                       SL$mean,SL$sd,NL$mean,NL$sd,p$modelType,
                       CV.Nc,CV.c,CV.pa,Dc,CV.Dc)
  names(result) <- c('season','siteCode','Nc','c','k','T','w','pa',
                     'SLmean','SLsd','NLmean','NLsd','modelType',
                     'CV.Nc','CV.c','CV.pa','Dc','CV.Dc')

  utils::write.csv(result, file = p$densityResultsFile, row.names=F)
  return(result)
}


### Supporting functions are all below
#
# Supporting functions----------------------------------------------------------

# ### $N_c$ (number of detections).
#
# $T$ and $N_c$ were calculated based, respectively, in the raw data duration
# and number of automatic detections correspondent to the full dataset.
#
## CountDetections--------------------------------------------------------------
#' Title
#'
#' @param detectionFile - File containing detections for the full dataset.
#'   File format must be tab separated with header, and must contain the
#'   following columns:
#'    t0: <double> Matlab datenum
#'    snr: <double> Signal to noise ratio of that detection (required only if
#'      snrTruncationThreshold is other than default value of -Inf)
#' @param season - Time period over which to subset detections. Can follow World
#'   Ocean Atlas numeric time codes (0-16), or be name of season
#'   ('summer','autumn','winter','spring'), month name/abbreviation, or 'year'
#'   (default).
#' @param snrTruncationThreshold <double> detections with SNR less than
#'   snrTruncation threshold (in dB) will be excluded from count
#' @param snrColName - name of column containing SNR values (only used if SNR
#'   truncation>-Inf)
#'
#' @returns <numeric> number of detections in the detectionFile within specified
#'   season and >= snrTruncationThreshold
#' @export
#'
countDetections <- function(detectionFile,
                            season='year',
                            snrTruncationThreshold=-Inf,
                            snrColName = 'snr') {
  # detectionFile must be
  det <- read.delim(detectionFile, sep = '\t', header = TRUE)

  if (snrTruncationThreshold != -Inf){ # Only if truncation requested
    if ( snrColName %in% colnames(det) ){
      det <- subset(det,det$snr>=snrTruncationThreshold)
    }  else {
      stop("Truncation by SNR requires column `snr` in fullYearDetectionCsv.\n\nSNR Truncation not applied\n\n")
    }
  }

  det$t <- mat2Rdate(det$t0)

  det$season <- time2season(det$t)
  det$month <- time2monthCode(det$t)

  det<- subsetByTimeCode(det,det$t,season)

  Nc <- dim(det)[1]
  return(Nc)
}

#  ### $T$ (time of deployment)
#
#  Sum the duration of the audio that has been analysed
#
## deploymentDuration-----------------------------------------------------------
#' Calculate deployment duration from recording start times and durations
#'
#' @param duration - list of recording durations (in s)
#' @param startDate - list of start dates for each duration
#' @param season - Subset the startDate by this season/month/year (default='year')
#'
#' @returns - duration of effort (as numeric)
#' @export
#'
deploymentDuration <- function(duration, startDate, season ='year'){

  wavInfo <- data.frame(startDate=startDate, duration=duration)
  wavInfo <- subsetByTimeCode(wavInfo,wavInfo$startDate,season)

  T <- sum(wavInfo$duration)/3600 # in hours
  return(T)
}

#' Calculate deployment duration from a Matlab soundFolder csv file.
#' @param fullYearEffortFile - csv file containing a row for each time period
#'   that was included in the study and columns named startDate and duration.
#'   startDate should be a matlab datenum, and duration should be numeric
#'   and in seconds.
#'
#' @param season - default='year'
#'
#' @export
deploymentDurationFromsoundFolderCsv <- function(fullYearEffortFile,
                                                season ='year'){
    # File generated from Matlab wavFolderInfo
    wavInfo <- read.csv(fullYearEffortFile)
    startDate <- mat2Rdate(wavInfo$startDate)
    duration <- wavInfo$duration
    deploymentDuration(duration, startDate, season)
}

#' #  ### Estimate $c$ (False discovery rate)
#'
#' False discovery rate, c, is defined as: FP/(FP+TP), where FP is the number of
#' false-positives and TP is the number of true positives
#' (https://en.wikipedia.org/wiki/Precision_and_recall; TODO: find a primary
#' literature citation to use instead of citing wikipedia).
#'
#' Here we estimate false-positive rate of the automated detector from the
#' manually annotated dataset. The data-file for this is the capture-history
#' table that contains reconciled annotated, and automated detections.
#'
#' First, we count the number of FP that the detector produced from the capture
#' history table. These are the sum of the rows where detect_table2==T &
#' detect_table1==F. Then we calculate the number of true positives, the sum of
#' the rows where detect_table2==T & detect_table1==T.
#'
#'
#' @param capHistFile file name of capture history table (TODO pass as
#'   data.frame instead of file).
#' @param season a character indicating the season or month for which to
#'   estimate call density
#' @param snrTruncationThreshold - Exclude rows with SNR below this threshold
#'   (in dB)
#' @param gtColName - Column name that contains the ground truth detections (0
#'   for false positive and 1 for true positive) Default='detect_table1'.
#' @param testColName - Column name of the detections under investigation (i.e.
#'   that contains the detections from which to calculate false positives).
#'   Default='detect_table2')
#' @param snrColName - Name of column that contains SNR estimates (only used if
#'   thresholding by SNR)
#'
#' @returns fdr, list containing false disovery rate, c, and it's CV. The false
#'   discovery rate c, is equal to (1-precision) for detector2 assuming
#'   detector1 is ground truth for detections within the specified season
#' @export
#'
falseDiscoveryRate <- function(ch,
                               season='year',
                               snrTruncationThreshold=NULL,
                               gtColName = 'detect_table1',
                               testColName = 'detect_table2',
                               snrColName = 'snr'
                               ){

  if (!is.null(snrTruncationThreshold)){
    snr <- ch[,snrColName]
    ch <- subset(ch, snr >= snrTruncationThreshold)
  }

  ch <- subsetByTimeCode(ch,ch$t,season)

  FP <- sum(ch[,testColName] & !ch[,gtColName])
  TP <- sum(ch[,testColName] & ch[,gtColName])
  c <- FP/(FP+TP)
  n <- FP + TP

  cv.c <- c_CV(n, c)

  fdr = list(c=c, cv.c = cv.c)

  return(fdr)
}


#' WARNING: This function is not supported, and probably does not do what
#' whatever you were hoping it might do. Perhaps consider function c_CV instead.
#' CV of false discovery rate using Cochran approximation
#'
#' @param hourlyFalsePosFile Name of csv file contianing hourly estimates of
#'   false positive rate.
#' @param season string or number corresponding to the time of year for which
#'   call densities should be estimated (and data subsetted).
#'
#' @returns CV of false discovery rate using Cochran approximation
#' @seealso c_CV
#' @export
#'
falseDiscoveryCV <- function(hourlyFalsePosFile,season){
  #  ### $CV_c$ (Cochran approx.);
  #
  #  First calculate using annotated library
  #
  fp <- read.csv(hourlyFalsePosFile)
  fp <- subsetByTimeCode(fp,lubridate::dmy_hms(fp$dt0), season)
  fp$fdr <- fp$num_fp/fp$num_detections
  x<-fp$fdr # false discovery rate
  wgts<-as.numeric(fp$num_detections) # w = n. of detections
  var.wtd.mean.cochran(x,wgts)

  #D.Harris: "Then, from the variance, you can estimate the CV for each c
  #estimate. This would be the square root of the variance, divided by the
  #square root of the number of samples (so estimate the standard error, then
  #divide by
  #the mean false discovery rate, to get the CV."

  SEc=sqrt(var.wtd.mean.cochran(x,wgts))/sqrt(sum(!is.na(x)))
  SEc

  CV.c=SEc/Hmisc::wtd.mean(x,wgts)
  CV.c
}

#' WARNING: This function is not supported, and probably does not do what
#' whatever you were hoping it might do. Use functions falseDiscoverRate and
#' c_CV instead.
#' False discovery rate from inspection of every Nth detection
#' @param falsePositiveXlsx Excel spreadsheet with timestamp of false positives.
#' This spreadsheet requires two columns: one called UTC and the other called
#' 'True Positive Rate.' These should be in a worksheet called 'conference'.
#' The file used by this function must already have had SNR truncation applied
#' (i.e. this function has no means of appling an snrTruncationThreshold)
#'
#' @param season Month or season over which to subset the data. Months can be
#' 01-12, and seasons can be 'summer','autumn','winter','spring', or 'year'.
#' @seealso falseDiscoveryRate, c_CV
#' @export
#'
falseDiscoveryRateFromNth<- function(falsePositiveXlsx,season='year'){
  ffp <- readxl::read_xlsx(falsePositiveXlsx,sheet = 'conference')

  ffp$dt0 <-lubridate::ymd_hms(ffp$UTC)
  ffp$season <- time2season(ffp$dt0)
  ffp$month <- time2monthCode((ffp$dt0))

  ffp<- subsetByTimeCode(ffp,ffp$dt0,season)
  if (season!= 'year'){
    if (season == 'summer' || season=='autumn' ||
        season=='winter' || season=='spring'){

      fp<-  ffp %>%
        dplyr::mutate(date_col = season ) %>%
        dplyr::group_by(date_col) %>%
        dplyr::summarize(detection.num = dplyr::n(),
                         true.pos.num = sum(`True positive rate`),
                         false.pos.num = detection.num-true.pos.num,
                         c = false.pos.num/(false.pos.num+true.pos.num)
        )
      c <- fp$c
    }else{ # not a season or year, so must be a month

      fp<-  ffp %>%
        dplyr::mutate(date_col = month ) %>%
        dplyr::group_by(date_col) %>%
        dplyr::summarize(detection.num = dplyr::n(),
                         true.pos.num = sum(`True positive rate`),
                         false.pos.num = detection.num-true.pos.num,
                         c = false.pos.num/(false.pos.num+true.pos.num)
        )
      c <- fp$c
    }
  } else {
    c = sum(!ffp$`True positive rate`)/dim(ffp)[1]
  }

  #  Now calculate CV of c where the analyst has inspected every 50th detection.
  #
  #  Here we group the detections by day within each season in order to estimate
  #  CVs.

  # Sample period here is date (24 hours of a single day of the year) I think?
  fp<-  ffp %>%
    dplyr::mutate(date_col = lubridate::date(dt0) ) %>%
    dplyr::group_by(date_col) %>%
    dplyr::summarize(detection.num = dplyr::n(),
                     true.pos.num = sum(`True positive rate`),
                     false.pos.num = detection.num-true.pos.num,
                     c = false.pos.num/(false.pos.num+true.pos.num)
    )

  x <- fp$c
  wgts <- fp$detection.num

  SEc=sqrt(var.wtd.mean.cochran(x,wgts))/sqrt(sum(!is.na(x)))
  SEc

  CV.c=SEc/Hmisc::wtd.mean(x,wgts)
  CV.c

  return(data.frame(c,CV.c))
}


#' Nc_CV: Coefficient of Variation (CV) of the number of detected calls (Nc).
#'
#' CV.Nc is calculated from the variance of Nc, which depends on the probability
#' of detection (pa). Here, Nc represents the number of detected calls, and
#' independence between detections is assumed.
#'
#' Rationale for this calculation Variance of a binomial process, pa, is: 1)
#' sigma^2= n * pa * (1-pa)
#'
#' We require n, the number of trials conducted, but we have Nc, the number of
#' positive detections. So we correct for false positives, then scale true
#' positives by pa to get total number of trials
#'
#' 2) n =  (Nc*(1-c)/pa; # n: the number of trials conducted We substitute the
#' right hand side back into 1) to get: 3) sigma^2 = (Nc*(1-c)/pa*pa*(1-pa);
#' dividing by and multiplying by pa is the same as multiplying by one, leaving
#' 4) sigma^2 = Nc * (1-c) * (1-pa)
#'
#' @param Nc - Number of calls detected in the dataset (see function
#'   countDetections())
#' @param pa - Average probability of detection in the study area (see function
#'   pDetInArea())
#' @param c - false discovery rate (proportion of false positive detections
#'   divided by the total number predicted positive)
#'
#' @returns cv.Nc - coefficient of variation of the number of detected calls
#' @export
#'
Nc_CV <- function(Nc,pa,c){
  var.Nc = Nc * (1-c) * (1-pa) # variance of NC

  # Calculate the standard deviation of Nc
  sd.Nc = sqrt(var.Nc)

  # Calculate the coefficient of variation of Nc where (Nc*(1-c)) is mean of the
  # binomial process
  cv.Nc <- sd.Nc / (Nc*(1-c))

  return(cv.Nc)
}


#' c_CV: Coefficient of Variation (CV) of false discovery rate (c).
#'
#' CV.c is calculated from the variance of c, which depends on the probability
#' of false alarm (c), a binomial process.
#'
#' Rationale for this calculation Variance of a binomial process, pa, is: 1)
#' sigma^2= n * c * (1-c)
#'
#' We require n, the number of trials conducted (i.e. how many detections were
#' inspected when estimating c) and the estimate of c
#'
#'
#' @param n - Number of trials used to estimate c
#' @param c - false discovery rate (proportion of false positive detections
#'   divided by the total number predicted positive)
#'
#' @returns cv.c - coefficient of variation of the false discovery rate
#' @export
#'
c_CV <- function(n,c){
  var.c = n * c * (1-c) # variance of c

  # Calculate the standard deviation of Nc
  sd.c = sqrt(var.c)

  # Calculate the coefficient of variation of c. Here n*c is mean of the
  # binomial process (number of trials multiplied by number of false positives)
  cv.c <- sd.c / (n*c)

  return(cv.c)
}



#'pa_CV: Coefficient of Variation (CV) of the probability of detection (pa).
#'
#'CV.pa is calculated using the standard deviation of the means and the sum of
#'the standard deviations of the transects. This approach assumes independence
#'between transects and homogeneity of error. The summation of SDs assumes a
#'conservative model.
#'
#'
#'@param pa.all.transects - matrix containing two columns and same number of
#'  rows as number of transects. The first column contains the mean pa and the
#'  second column contains the standard deviation of pa for each transect. The
#'  final column contains the overall mean and sd.
#'
#'@returns CV.pa - CV of the probability of deteciton in the area
#'@export
pa_CV <- function(pa.all.transects){
# The above function, pDetInArea writes results to a bunch of files
# Load the file that we need that has the p_a in them, and ignore the others

no.transects <- dim(pa.all.transects)[1]-1;

pa.all.transects<- as.data.frame(pa.all.transects)

names(pa.all.transects)<-c('Mean','SD')

# Radials are in rows
y=pa.all.transects[1:no.transects,]

## IMPORTANT: to check if the var, SE and CV calculation is correct ##
varp.pa<-((sd(y$Mean)/sqrt(no.transects))^2)+((sum(y$SD)/no.transects)^2)
sep.pa<-sqrt(varp.pa)
CV.pa<-sep.pa/mean(y$Mean)
return(CV.pa)
}

#' Coefficient of Variation (CV) for call density (Dc)
#'
#' The uncertainty (CV) in the estimated call density (Dc) is calculated using
#' the Delta Method, which combines the CV of the individual parameters (Nc, c
#' and pa).
#'
#' Here, the uncertainties of Nc, c, and pa are assumed to be independent, given
#' the fundamental premise of the Delta Method in this application.
#'
#' @param CV.Nc - CV of number of calls Nc
#' @param CV.pa - CV of probability of detection in the area, pa
#' @param CV.c - CV of the false discovery rate, c
#'
#' @returns CV.Dc - coefficient of variation of the call density
#' @export
#'
Dc_CV<- function(CV.Nc,CV.pa,CV.c){
  CV.Dc<-sqrt((CV.Nc^2)+(CV.pa^2)+(CV.c^2))
  return(CV.Dc)
}


#' Call Density Estimate
#' title: "Call Density"
#' author: "Brian Miller"
#' date: "2022-10-12"
#'
#' This function estimates call density with all parameters included as function
#' arguments. This is in contrast to the previous operation, i.e. where a
#' data.frame of parameters was used to indicate which files to load from disk.
#' Here no parameter file/data.frame is used, so all data and parameters are
#' must already be objects in memory;
#'
#' Density of blue whale D-calls (BmD) in long-term moored recording dataset.
#' Apply methods of Castro, Harris, et al (in prep) to obtain call density.
#'
#' $D_c = \\frac{N_c*(1-c)}{kTP_a{\\pi}w^2}$
#
#' where:
#' $D_c$ is call density
#' $N_c$ is number of calls
#' $c$ is false discovery rate
#' $k$ is number of sensors (here always 1)
#' $T$ is the duration of data analysed
#' $P_a$ is the probability of detection in the study area
#' ${\\pi}w^2$ is the study area (in km\\^2)
#'
#' @param Nc - Total count of true positive detections for call density estimate
#' @param capHistTab - Capture history table used to derive false positive rate.
#'   Can also be used to derive probability of detection as a function of SNR if
#'   snrDetFun is not specified.
#' @param SL - List containing distribution of source level parameters.
#'   SL must contain named elements named: mean, sd, and sampleSize.
#' @param TL - Data.frame of transmission losses for Monte Carlo simulation. The
#'   first column must contain the ranges, the remaining columns contain
#'   transmission losses (in dB) for each radial transect at that range.
#' @param A - Scalar indicating the area for the call density estimate
#'   (default=1)
#' @param modelType - Either 'glm','gam','scam',or 'vglm' indicating which type
#'   of model to use for estimating probability of detection as a function of
#'   SNR (default='scam')
#' @param numKnots  - number of knots to use wor modelling detection function vs
#'   SNR (only if model types 'gam' or 'scam') default = 5
#' @param output.resolution.m Spatial resolution along radial transects for
#'   estimating probability of detection in area. Default = 100
#' @param outerloop - Number of bootstrap iterations for Monte Carlo simulation
#' @param transectFile - Optional name of text file to store probability of
#'   detection at each output.resolution.m step (potentially yielding a large file)
#' @param simResultsFile - Optional name of text file to store Monte Carlo
#'   Simulation parameters and results
#' @param paFile - Optional name of text file to store mean probability of
#'   detection per transect and overall area
#' @param season TimeCode specifying month, season, or year for outputs
#' @param snrDetFun OPTIONAL linear-model like structure (GLM,GAM,SCAM,etc)
#'    specifying the SNR-detection function to use. If this is not included,
#'    then the SNR-detection function will be derived from the capture history
#'    table.
#' @param truncationDistance scalar or matrix of truncation distances.
#'    If a matrix is provided, then the dimensions should be 1xN with N being
#'    the same as the number of transects
#' @param snrTruncationThreshold scalar SNR in dB below which the probability of
#'    detection will be forcibly set to zero.
#' @param NL list or data.frame containing distribution of noise level
#'   parameters. Must contain items named: mean, sd, and sampleSize (similar to
#'   SL).
#' @param siteCode Optional string containing a code or label associated with
#'   the site
#' @param k Number of sensors (default=1)
#' @param densityResultsFile Optional name of csv file where final call density
#'   results will be written as a data.frame
#'
#' @returns data.frame containing call density inputs, results, and CVs
#'
#' @importFrom stats rnorm sd vcov
#' @importFrom utils read.csv write.table
#' @importFrom magrittr %>%
#' @export
#'
cde <- function (Nc,
                 capHistTab,
                 SL,
                 TL,
                 T = 1,
                 A = 1,
                 k=1,
                 season='year',
                 ### Monte carlo simulation params
                 snrDetFun=NULL,
                 NL=NULL,
                 modelType = 'scam',
                 numKnots = 5,
                 output.resolution.m=100,
                 outerloop = 10,
                 transectFile='',  # file output names
                 simResultsFile='',
                 paFile='',
                 ###
                 truncationDistance=Inf,
                 snrTruncationThreshold=-Inf,
                 siteCode='',
                 densityResultsFile=''
                 ){
  # Check inputs and create outputs
  # Store all parameters for call-density estimation in data frame called 'p'
  #
  # InputCheck------------------------------------------------------------------



  # Number of calls
  # Nc--------------------------------------------------------------------------
  # Nc <- countDetections(p$detectorParams$fullYearDetectionCsv,season,
  #                       snrTruncationThreshold=snrTruncationThreshold)

  # Duration of time monitored
  #T----------------------------------------------------------------------------
  # T <- deploymentDuration(p$detectorParams$fullYearEffortFile, season)

  #  ### Calculate study area
  # A---------------------------------------------------------------------------
  # A = studyArea(p$w,truncationDistance)

  # c, CV_c---------------------------------------------------------------------
  fdr <- falseDiscoveryRate(capHistTab, season, snrTruncationThreshold)

  # CV.c and c are packed in a list, so unpack
  CV.c <- fdr$cv.c
  c <- fdr$c


  #   ### $p_a$ (Overall probability of detection)
  # This is the most complicated part and takes the longest time
  # Pa--------------------------------------------------------------------------
  SNRinfo <- capHist2snrInfo(capHistTab,season)

  # Read SL distribution parameters from overall parameters.
  # SL <- data.frame(mean=p$slParams$slMean,
  #                  sd=p$slParams$slStd,
  #                  sampleSize=p$slParams$slSampleSize)

  # TL<-utils::read.csv(p$tlParams$tlFile)

  # Check whether user supplied an snr detection function or estimate from data
  # NB: This needs to occur AFTER any SNR truncation.
  if (is.null(snrDetFun)){
    # Estimate snrDetFun from SNRinfo
    snrDetFun <- fitSNRdetectionFunc(
      subset(SNRinfo,SNR>=snrTruncationThreshold),
      modelType=modelType, numKnots)
  }

  # If user has not specified NL distribution (data.frame with columns mean, sd,
  # samplesize), then extract this information it from SNRinfo.
  # NB: This needs to occur with untruncated SNRInfo (i.e. prior to truncation).
  if (is.null(NL)){
    NL <- nlFromSnrInfo(SNRinfo, snrDetFun)
  }

  pDetResults <- pDetInArea(snrDetFun, SL, TL,  NL, # Sonar equation inputs
                            output.resolution.m=p$output.resolution.m,
                            outerloop=p$outerloop,
                            truncationDistance=truncationDistance,
                            snrTruncationThreshold= snrTruncationThreshold,
                            transectFile,  # file output names
                            simResultsFile,
                            paFile)

  # The above function, pDetInArea writes results to a bunch of files
  # Load the file that we need that has the p_a in them, and ignore the others
  # pa.all.transects <- read.csv(p$paFile, header = TRUE, sep = ' ');
  pa.all.transects <- pDetResults$perTransectMeanSD
  no.transects <- dim(pa.all.transects)[1]-1;
  pa <- pa.all.transects[no.transects+1,1]

  #CV_Pa------------------------------------------------------------------------
  # The line below will read the per-transect pa from file.
  # pa.all.transects <- read.csv(p$paFile, header = TRUE, sep = ' ');

  # Per-transect mean and SD of pa now returned in the list of pDetResults, so
  # no need to read or write it to a file
  pa.all.transects <- pDetResults$perTransectMeanSD
  CV.pa <- pa_CV(pa.all.transects)

  #CV_Nc------------------------------------------------------------------------
  CV.Nc <- Nc_CV(Nc,pa,c)

  #  ### Finally, estimate Dc (Density of calls) [calls $h^{-1} km^{-2}$]
  # Dc --------------------------------------------------------------------------
  Dc <- (Nc * (1-c) )/( k * A * pa * T )
  Dc     # per h per km^2
  Dc*1e3 # per h per 1000 km^2?

  # CV_Total--------------------------------------------------------------------
  CV.Dc <- Dc_CV(CV.Nc,CV.pa,CV.c)

  #  ## Collate into data frame and write to csv file
  # Format output---------------------------------------------------------------

  result <- data.frame(season, siteCode,Nc,c,k,T,A,pa,
                       SL$mean,SL$sd,NL$mean,NL$sd,modelType,
                       CV.Nc,CV.c,CV.pa,Dc,CV.Dc)
  names(result) <- c('season','siteCode','Nc','c','k','T','A','pa',
                     'SLmean','SLsd','NLmean','NLsd','modelType',
                     'CV.Nc','CV.c','CV.pa','Dc','CV.Dc')

  if (!is.null (densityResultsFile)){
    utils::write.csv(result, file = densityResultsFile, row.names=F)
  }
  return(result)
}

