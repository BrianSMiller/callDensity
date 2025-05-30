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
#'
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
#'
nlFromSnrInfo <- function(snrInfo, snrDetFun){
NL = snrInfo %>% dplyr::summarise(mean=mean(NoiseRL,na.rm = TRUE),
                             sd=sd(NoiseRL,na.rm = TRUE),
                             sampleSize=dplyr::n()-sum(is.na(NoiseRL)))

# Use the SNR intercept estimate from the VGLM to add back into the NL
# distribution. Have resulted to a clunky kludge since I can't figure out how to
# get the intercept on the response scale. Instead, just predict the
# probabilities on the response scale for a grid of SNR values and lookup the
# median value.
snrGrid = data.frame(SNR=seq(from=-20, to=20, by=0.1))
probs <- predict(snrDetFun,newdata=snrGrid,type='response')

# VGLMs are handled differently than other models
if (any(class(snrDetFun)=='vglm')){
  predAny <- VGAM::predict(snrDetFun,newdata=snrGrid,type="response",
                           type.fitted='onempall0')
  predDet2 <- VGAM::predict(snrDetFun,newdata=snrGrid,type="response")

  # VGLMs can have multiple observers, so we need to know which of these
  # to use for probability of detection.
  index = ifelse(is.null(snrDetFun@extra$whichObserver),
                 dim(predDet2)[2],
                 which(colnames(predDet2)==snrDetFun@extra$whichObserver) )

  predDet2 <- predDet2[,index]
  probs = apply(cbind(predDet2,predAny),1,prod)
}

snrOffset <- snrGrid[which.min(abs(probs-0.5)),]

NL$mean <- NL$mean+snrOffset

return(NL)

}

#' Convert a capture history FILE of multiple detections into an SNR.Info file
#'
#' @param capHistFile A capture history file containing detections and SNR info
#' @param season A timeCode corresponding to months, seasons, or year
#'
#' @export
capHist2snrInfo <- function(snr,season='year'){

  # remove false positives
  snr <- snr[snr$detect_table1==1,]

  #Check for duplicate lines:
  # class(snr) # dataframe
  dup<-duplicated(snr$key,fromLast = TRUE)

  ### Z Call: If only work with z calls
  #snr<-subset(snr, classification == "BmAntZ")

  # CallRL and NoiseRL here in dB, but do NOT need to be re 1 uPa (09-Nov-22)
  Detected<-snr$detect_table2
  CallRL<-as.numeric(snr$signalRMSdB) #dB
  NoiseRL<-as.numeric(snr$noiseRMSdB) #dB

  # SNR<-as.numeric(snr$signalRMSdB-snr$noiseRMSdB) # SNR measured as in DH
  SNR<-CallRL-NoiseRL # to check
  SNRinfo<-data.frame(Detected,CallRL,NoiseRL,SNR,snr$t,snr$month,snr$season)
  names(SNRinfo)[5:7]<- c('t','month','season')
  # Define some sample sizes needed inside the loop

  # Check for SNR=Inf: put in a warning to user here, since infinite SNR is likely
  # indicative of a problem with the inputs
  lines<-c(which(SNRinfo$SNR=="Inf", arr.ind = TRUE))
  if ( length(lines) >0) { # Run this line if inf occurs
    warning( sprintf('%g of the SNR were infinite:\n',length(lines)) )
    SNRinfo <- SNRinfo[-lines,]
  }

  if (season=='year' || season=='00'){
    return(SNRinfo)
  } else if (season == 'summer' || season=='autumn' ||
             season=='winter' || season=='spring'){
    SNRinfo <- SNRinfo[SNRinfo$season==season,]
  } else {
    SNRinfo <- SNRinfo[SNRinfo$month==season,]
  }

  return(SNRinfo)
}

#' Convert capture history DATA.FRAME into the 'SNRinfo' format used
#' by the callDensity package
#' capHistTosnrInfo
#'
#' @param capHistTab - Capture history table of detections
#'
#' @returns - SNRInfo data.frame containing SNR, detections, RL, NL
#' @export
#'
capHistTosnrInfo <- function(capHistTab){

  # Assume detect_table1 is ground truth, so only keep rows where
  # detect_table1==1 (includes both true and false positives though)
  capHistTab <- capHistTab[capHistTab$detect_table1==1, ]

  # Detected refers to the automated detector, here detect_table2
  Detected<-capHistTab$detect_table2

  # Add time rows for compatibility with callDensity format
  t <- capHistTab$datetime
  season <- time2season(t)
  month <- time2monthCode(t)

  # CallRL and NoiseRL here in dB
  CallRL<-as.numeric(rowMeans(
    subset(capHistTab,select=c('signalRMSdB1','signalRMSdB2')) ,na.rm=T) )#dB
  NoiseRL<-as.numeric(rowMeans(
    subset(capHistTab,select=c('noiseRMSdB1','noiseRMSdB2')) ,na.rm=T) )#dB

  # SNR<-as.numeric(capHistTab$signalRMSdB-capHistTab$noiseRMSdB) # SNR measured as in DH
  SNR<-CallRL-NoiseRL # to check
  SNRinfo<-data.frame(Detected,CallRL,NoiseRL,SNR,t,month,season)

  # Check for SNR=Inf: put in a warning to user here, since infinite SNR is
  # likely indicative of a problem with the inputs
  lines<-c(which(SNRinfo$SNR=="Inf", arr.ind = TRUE))
  if ( length(lines) >0) { # Run this line if inf occurs
    warning( sprintf('%g of the SNR were infinite:\n',length(lines)) )
    SNRinfo <- SNRinfo[-lines,]
  }

  return(SNRinfo)
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
fitSNRdetectionFunc <- function(SNRinfo, modelType='gam', numKnots=3){
  modelType <- tolower(modelType)

  if (modelType == 'glm'){
    res.1=stats::glm(Detected~(SNR), data = SNRinfo,
              family = binomial(link="logit")) # or: family quasibinomial
  }
  if (modelType == 'scam'){
      res.1<-scam::scam(Detected~s(SNR, k=numKnots,bs="mpi"),
                        data=SNRinfo,family=binomial)
  }
  if (modelType == 'gam'){
      res.1<-mgcv::gam(Detected~s(SNR, k=numKnots), data=SNRinfo,
                       family=binomial)
  }
    return(res.1)
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
                       yColNames = c('detect_observer1','detect_observer2'),
                       whichObserver='detect_observer2'){


  res.1 <- VGAM::vglm(as.matrix(SNRinfo[,yColNames]) ~SNR,
                posbernoulli.t(parallel.t = TRUE~0), data=SNRinfo )
  res.1@extra$whichObserver <- whichObserver

  return(res.1)
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
fitSNRvgam <- function(SNRinfo,
                       yColNames = c('detect_observer1','detect_observer2'),
                       whichObserver='detect_observer2'){

  res.1 <- VGAM::vgam(as.matrix(SNRinfo[,yColNames]) ~ s(SNR),
                      posbernoulli.t, data=SNRinfo,
                      na.action = 'na.omit')
  res.1@extra$whichObserver <- whichObserver

  return(res.1)
}


#' Plot SNR-detection function (ggplot2)
#' @param res.1 An SNR-detection function from fitSNRdetectionFunc
#'
#' @export
showSNRdetectionFunc <-  function(res.1){

  p1 <- marginaleffects::plot_predictions(res.1, condition='SNR',
                                  type='response',conf_level=0.95)+
    ggplot2::labs(x = c("SNR (dB)"), y = "P(detecting a call)" , title = NULL)+
      ggplot2::theme(legend.position="bottom",
                     legend.key.height = unit(0.25,"cm"),
                     legend.key.width = unit(0.25,"cm"),
                     legend.text = element_text(size = 6))
    summary(res.1)
    return(p1)
}

#' Model SNR-detection function with season as a factor (deprecated?)
#'
#' @param SNRinfo Data.frame containing columns SNR, Detected, and season
#'
#' @param season Factor containing a timeCode
#' @param useGLM Boolean, if true fit a GLM instead of GAM
#' @param numKnots Number of knots to use in the GAM
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

#' Load specified density results text files in a directory
#'
#' @returns Data.frame of all call density results in the specified path
#'
#' @param path Location (path) of text files containing call density results
#'
#' @param siteCode Prefix of _density files to be loaded
#'
#' @export
collectDensityResults <- function (path, siteCode){
  files <- dir(path=path, pattern = paste0('^density_',siteCode),
               full.names = TRUE)
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
    commonTable <- subset(d[1,], select=c('k','w','SLmean','SLsd'))
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
#' @param truncationDistance - Matrix of truncation distances [1 x Nr]. Nr is
#'   the number of radial transects in the simulation (i.e. the same as the
#'   number of TL profiles).
#'
#' @returns Area of the simulation
#' @export
#'
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
#' @param capHistFile - name of the csv or txt file that contains capture
#'   histories. TODO: Document the required columns/file format better
#'
#' @returns capture history table as data.frame
#' @export
#'
readCapHist <- function(capHistFile){
  ch <- read.csv(file = p$capHistFile,na.strings = c('NaN','NA'))
  ch <- capHistTimeSeason(ch)
  return(ch)
}

