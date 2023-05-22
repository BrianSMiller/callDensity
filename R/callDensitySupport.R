library("stats",quietly = T, warn.conflicts = F)
library("dplyr",quietly = T, warn.conflicts = F)
library("ggplot2",quietly=T,warn.conflicts = F)
library("kableExtra", quietly=T, warn.conflicts = F)
library("lubridate", quietly = T, warn.conflicts = F)
library("Hmisc",quietly = T, warn.conflicts = F)

#' Convert matlab datenum to R POSIXct
#'
#' @param x List of Matlab datenums
#' @returns posixCT corresponding to the Matlab datenum that was input
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
#'
#' @returns Parameter file with file names corresponding to the correct timeCode
#'
#' @export
defaultOutputFileNames <- function(p,season){

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
  p$SLmean <- 189;
  p$SLsd <- 8.0;
  p$SLsamplesize <- 350;

  # Whether to use GLM (faster) or GAM (better)
  p$useGLM <- TRUE;
  # Whether to constrained gam from scam or plain mgcv::gam
  p$useSCAM <- TRUE;
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
  p$paFile <- paste0("Pa_", p$siteCode, "_", season, ".txt")
  p$transectFile <- paste0("pDet_",p$siteCode,'_',season,'_transects.csv')
  p$simResultsFile <- paste0('simResult_', p$siteCode, '_',season,'.txt')

  ## Overall call density output file name
  p$densityResultsFile <- paste0('density_',p$siteCode, '_', season,'.csv')

  return(p)
}

#' Mean and standard deviation of noise levels by month, season, or year.
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

#' Convert a capture history table of multiple detections into an SNR.Info file
#'
#' @param capHistFile A capture history file containing detections and SNR info
#' @param season A timeCode corresponding to months, seasons, or year
#'
#' @export
capHist2snrInfo <- function(capHistFile,season='year'){
  snr<-utils::read.csv(capHistFile,na.strings = c('NaN','NA'))
  snr <- capHistTimeSeason(snr)

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
#' @param useGLM Boolean - if true a GLM will be fit instead of a GAM
#' @param useSCAM Boolean - if true a shape constrained model will be fit instead of a GAM
#' @param numKnots Number of knots in the GAM or SCAM
#'
#' @importFrom stats binomial coef predict quasibinomial rnorm sd vcov
#'
#' @export
fitSNRdetectionFunc <- function(SNRinfo, useGLM=TRUE, useSCAM=TRUE, numKnots=3){
  # Initial glm:

  if (useGLM==TRUE){
    res.1=stats::glm(Detected~(SNR), data = SNRinfo,
              family = binomial(link="logit")) # or: family quasibinomial
  } else {
    if (useSCAM==TRUE){
      res.1<-scam::scam(Detected~s(SNR, k=numKnots,bs="mpi"),
                        data=SNRinfo,family=binomial)
    } else {
      res.1<-mgcv::gam(Detected~s(SNR, k=numKnots), data=SNRinfo,
                       family=quasibinomial) # Test to include k=3
    }
    return(res.1)
  }
}

#' Plot SNR-detection function (ggplot2)
#' @param res.1 An SNR-detection function from fitSNRdetectionFunc
#'
#' @export
showSNRdetectionFunc <-  function(res.1){

  require(ggplot2)
  require(marginaleffects)
  p1 <- marginaleffects::plot_cap(res.1, condition='SNR',
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
#' @export collect.density.results
collect.density.results <- function (path, siteCode){
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
#' @export collect.density.results
show.density.input.table <-function(d, separateCommon=TRUE){

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
#' @export density.results.table
density.results.table<- function (d){

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
#' @export detection.rate.plot
detection.rate.plot <- function(d){

  gRate <- ggplot2::ggplot(data=d, aes(x=season, y=Nc/T, fill=season))+
    ggplot2::geom_bar(stat="identity",show.legend=F)+
    # scale_y_continuous(limits=c(0,15))+
    # geom_errorbar(aes(ymin=CI.low, ymax=CI.high), width=0.5, position=position_dodge(0.9))+
    ggplot2::ylab( TeX("Mean detection rate ($ calls \\cdot h^{-1})" ) )+
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
#' @export detection.rate.corrected.plot
detection.rate.corrected.plot <- function(d){

  gRateCorrected <- ggplot2::ggplot(data=d, aes(x=season, y=(Nc*c)/T, fill=season))+
    ggplot2::geom_bar(stat="identity",show.legend=F)+
    # scale_y_continuous(limits=c(0,15))+
    # geom_errorbar(aes(ymin=CI.low, ymax=CI.high), width=0.5, position=position_dodge(0.9))+
    ggplot2::ylab( TeX("Mean detection rate ($ calls \\cdot h^{-1})" ) )+
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
#' @export density.plot
density.plot <- function(d){

  gDen <- ggplot2::ggplot(data=d, aes(x=season, y=Dc, fill=season))+
    ggplot2::geom_bar(stat="identity",show.legend=F)+
    ggplot2::geom_errorbar(aes(ymin=CI.low, ymax=CI.high), width=0.5,)+
    # scale_y_continuous(limits=c(0,0.25))+
    ggplot2::ylab(TeX("Call density ($ calls \\cdot h^{-1} \\cdot 1000 \\cdot km^{-2}$)" ) )+
    ggplot2::xlab('')+
    ggplot2::labs(tag="(B)") +
    ggplot2::theme_minimal()+
    ggplot2::theme(plot.tag = element_text(),
          plot.tag.position = c(0.225, 0.975))
}


