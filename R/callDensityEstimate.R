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
cde <- function (p,season, snrDetFun=NULL, truncationDistance=Inf,
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

  #   ### Nc, T, c, CV_c
  #
  # Nc, T, c, CV_c--------------------------------------------------------------
  Nc <- countDetections(p$detectorParams$fullYearDetectionCsv,season,
                        snrTruncationThreshold=snrTruncationThreshold)

  T <- deploymentDuration(p$detectorParams$fullYearEffortFile, season)

  c <- falseDiscoveryRate(p,season)
  CV.c <- falseDiscoveryCV(p,season)

  #  alternative measure of FDR where analyst inspected every nth detection
  if (p$useSeparateFPdata){
    c <- falseDiscoveryRateFromNth(p$manualFPfileName,season)
    # CV.c and c are packed in a data frame, so unpack
    CV.c <- c$CV.c
    c <- c$c
  }


  #   ### $p_a$ (Overall probability of detection)
  #
  # Pa--------------------------------------------------------------------------
  SNRinfo <- capHist2snrInfo(p$capHistFile,season)

  # If user has not specified NL distribution (data.frame with columns mean, sd,
  # samplesize), then extract this information it from SNRinfo
  if (is.null(NL)){
      NL <- SNRinfo %>% dplyr::summarise(mean=mean(NoiseRL,na.rm = TRUE),
                                     sd=sd(NoiseRL,na.rm = TRUE),
                                     sampleSize=dplyr::n()-sum(is.na(NoiseRL)))
  }

  SL <- data.frame(mean=p$slParams$slMean, sd=p$slParams$slStd, sampleSize=p$slParams$slSampleSize)

  TL<-utils::read.csv(p$tlParams$tlFile)

  # Check whether user supplied an snr detection function or estimate from data
  if (is.null(snrDetFun)){
  # Estimate snrDetFun from SNRinfo
    snrDetFun <- fitSNRdetectionFunc(SNRinfo,
                                modelType=p$modelType, p$numKnots)
  }

  pDetInArea(snrDetFun, SL, TL,  NL, # Sonar equation inputs
             p$transectFile, p$simResultsFile, p$paFile, # file output names
             output.resolution.m=p$output.resolution.m, outerloop=p$outerloop,
             truncationDistance=truncationDistance,
             snrTruncationThreshold= snrTruncationThreshold)

  # The above function, pDetInArea writes results to a bunch of files
  # Load the file that we need that has the p_a in them, and ignore the others
  pa.all.transects <- read.csv(p$paFile, header = TRUE, sep = ' ');
  no.transects <- dim(pa.all.transects)[1]-1;
  pa <- pa.all.transects[no.transects+1,1]
  pa

  #  ### $CV_{P_a}$ (see Harris, 2012 for details);
  #
  #  Approach is to calculate CV using mean and SD from each radial?
  #
  #CV_Pa------------------------------------------------------------------------
  # Coef. of variation - Pa

  names(pa.all.transects)<-c('Mean','SD')

  # Radials are in rows
  y=pa.all.transects[1:no.transects,]

  ## IMPORTANT: to check if the var, SE and CV calculation is correct ##
  ## Per year:

  varp.pa<-((sd(y$Mean)/sqrt(no.transects))^2)+((sum(y$SD)/no.transects)^2)
  sep.pa<-sqrt(varp.pa)
  CV.pa<-sep.pa/mean(y$Mean)
  CV.pa

  #  ### $CV_{N_c}$ (see Harris, 2012 for details);
  # Confirm how to calculate SD and CV (?)
  # CV_Nc-----------------------------------------------------------------------
  varn=Nc*(pa)*(1-pa)
  SDn=sqrt(varn)
  CV.Nc=SDn/Nc

  #  ### Calculate study area
  #
  # At most baisc study area is a circle of pi*w^2 radius. If different
  # truncation distances are used for each transect, then need to estimate
  # the area of each transect as sectors of a circle of r=truncation distance
  # A---------------------------------------------------------------------------

  A = pi*p$w^2
  if (any(truncationDistance < p$w)     ){
     if (length(truncationDistance) == 1) {
       A = pi*truncationDistance^2
     } else {
       A = sum(pi*truncationDistance^2/length(truncationDistance))
     }


  }

  #  ### Finally, estimate Dc (Density of calls) [calls $h^{-1} km^{-2}$]
  #
  # Dc --------------------------------------------------------------------------

  Dc <- (Nc * (1-c) )/( p$k * A * pa * T )
  Dc     # per h per km^2
  Dc*1e3 # per h per 1000 km^2?

  #  ### CV_Total (Delta Method)
  #
  # CV_Total--------------------------------------------------------------------

  CV.Dc<-sqrt((CV.Nc^2)+(CV.pa^2)+(CV.c^2))
  CV.Dc


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
#'   ('summer','autumn','winter','spring'), month name/abbreviation, or 'year'.
#' @param snrTruncationThreshold <double> detections with SNR less than
#'   snrTruncation threshold will be excluded from count
#'
#' @returns <numeric> number of detections in the detectionFile within specified
#'   season and >= snrTruncationThreshold
#' @export
#'
countDetections <- function(detectionFile,season,snrTruncationThreshold=-Inf) {
  # detectionFile must be
  det <- read.delim(detectionFile, sep = '\t', header = TRUE)

  if (snrTruncationThreshold != -Inf){ # Only if truncation requested
    if ( "snr" %in% colnames(det) ){
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
  Nc
}

#  ### $T$ (time of deployment)
#
#  Sum the duration of the audio that has been analysed
#
## deploymentDuration-----------------------------------------------------------
deploymentDuration <- function(fullYearEffortFile, season){
  # File generated from Matlab wavFolderInfo
  wavInfo <- read.csv(fullYearEffortFile)
  wavInfo <- subsetByTimeCode(wavInfo,mat2Rdate(wavInfo$startDate),season)

  T <- sum(wavInfo$duration)/3600 # in hours
  T
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
#' First, we count the number of FP that the AI produced from the capture
#' history table. These are the sum of the rows where detect_table2==T &
#' detect_table1==F. Then we calculate the number of true positives, the sum of
#' the rows where detect_table2==T & detect_table1==T.
#'
#' Alternatively, we could use the double observer mark-recapture to estimate
#' the total number of calls in the dataset (see Miller et al 2022 - AI bests
#' human observer - RSEC). This would require incorporating the Huggins/RMark
#' script that can be found in the supplement to the Miller et al 2022 paper.
#'
#' @param p data.frame containing top-level parameters for a call density
#'   estimate
#' @param season a character indicating the season or month for which to
#'   estimate call density
#'
#' @returns c, the false discovery rate (1-precision) for detector2 assuming
#'   detector1 is ground truth for detections within the specified season
#' @export
#'
falseDiscoveryRate <- function(p,season){

  effort <- readxl::read_excel(p$annotationParams$effortFile);
  effort$season <- time2season(effort$StartTime)
  effort$month <- time2monthCode(effort$StartTime)
  ch <- read.csv(file = p$capHistFile)
  ch <- capHistTimeSeason(ch)

  subsetByTimeCode(effort,effort$StartTime,season)
  subsetByTimeCode(ch,ch$t,season)

  FP <- sum(ch$detect_table2 & !ch$detect_table1)
  TP <- sum(ch$detect_table2 & ch$detect_table1)
  c <- FP/(FP+TP)
  c
}

#' CV of false discovery rate using Cochran approximation
#'
#' @param p data.frame containing all parameters for call density estimate
#' @param season string or number corresponding to the time of year for which
#'   call densities should be estimated (and data subsetted).
#'
#' @returns CV of false discovery rate using Cochran approximation
#' @export
#'
falseDiscoveryCV <- function(p,season){
  #  ### $CV_c$ (Cochran approx.);
  #
  #  First calculate using annotated library
  #


  fp <- read.csv(p$hourlyFalsePosFile)
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

#' False discovery rate from inspection of every Nth detection
#' @param falsePositiveXlsx Excel spreadsheet with timestamp of false positives.
#' This spreadsheet requires two columns: one called UTC and the other called
#' 'True Positive Rate.' These should be in a worksheet called 'conference'.
#'
#' @param season Month or season over which to subset the data. Months can be
#' 01-12, and seasons can be 'summer','autumn','winter','spring', or 'year'.
#'
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
