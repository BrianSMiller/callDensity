
library(dplyr,quietly = T, warn.conflicts = F)

source('probabilityOfDetectionInArea.R')

# Convert matlab datenum to R POSIXct
mat2Rdate <- function (x) as.POSIXct((x - 719529)*86400, origin = "1970-01-01",
                                     tz = "UTC")

# Lookup the season for a given POSIXct, x. Seasons are defined by 3 month 
# periods from Summer=(Dec,Jan,Feb); Autumn=(Mar,Apr,May), Winter=(Jun,Jul,Aug), 
# Spring=(Sep,Oct,Nov)
time2season <- function(x) {
  mon<-lubridate::month(x)
  seasons <- rep(c("winter","spring","summer","autumn"), each = 3)
  names(seasons) <- month.abb[c(6:12,1:5)]
  season <- factor(seasons[month.abb[mon]], 
                   levels = c("summer","autumn","winter","spring"))
}

# Lookup the season for a given POSIXct, x. Seasons are defined by 3 month 
# periods from Summer=(Dec,Jan,Feb); Autumn=(Mar,Apr,May), Winter=(Jun,Jul,Aug), 
# Spring=(Sep,Oct,Nov)
time2monthCode <- function(x) {
  mon<-lubridate::month(x)
  monthCodes <- sprintf("%02d",mon)
  monthCodes <- factor(monthCodes, 
                   levels = c('01','02','03','04','05','06',
                              '07','08','09','10','11','12'))
}



# Add R datetime and season to a capture history table
capHistTimeSeason <- function(x){
  x$t <- mat2Rdate(x$t0_table1)
  x$t[is.na(x$t)] <- mat2Rdate(x$t0_table2[is.na(x$t)])
  x$season <- time2season(x$t)
  x$month <- time2monthCode(x$t)
  # }
  return(x)
}

# Add R datetime and season to a capture history table (version 2)
capHistTimeSeason2 <- function(x){
  x$t <- mat2Rdate(x$t0)
  x$season <- time2season(x$t)
  x$month <- time2monthCode(x$t)
  return(x)
}

# Given data.frame of call density parameters, p, create default parameters and
# output file-names for Monte-Carlo modelling and estimating call-density. 
# Append all of these parameters to the data.frame and return it.
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
  p$tlFile <- paste0('TL_', siteCode, '_woa2018_meanSSP_', 
                     p$output.resolution.m, 'm_res_', season, '.csv')
  p$tlFile <- ifelse(file.exists(p$tlFile), p$tlFile,
                     paste0('TL_', siteCode, '_woa2018_meanSSP_', season, '.csv') )
  
  ### Outputs
  ## Output file names of Monte-Carlo model of Pa
  p$paFile <- paste0("Pa_", p$siteCode, "_", season, ".txt")
  p$transectFile <- paste0("pDet_",p$siteCode,'_',season,'_transects.csv')
  p$simResultsFile <- paste0('simResult_', p$siteCode, '_',season,'.txt')
  
  ## Overall call density output file name
  p$densityResultsFile <- paste0('density_',p$siteCode, '_', season,'.csv')

  return(p)
}

noiseLevelDistribution <- function(nlFile,season='year'){
  nl <- read.csv(nlFile)
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

capHist2snrInfo <- function(capHistFile,season='year'){
  snr<-read.csv(capHistFile,na.strings = c('NaN','NA'))
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


fitSNRdetectionFunc <- function(SNRinfo, useGLM=TRUE, useSCAM=TRUE, numKnots=3){
  # Initial glm:

  if (useGLM==TRUE){
    res.1=glm(Detected~(SNR), data = SNRinfo,
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
  
showSNRdetectionFunc <-  function(res.1){
  
  # require(ggplot2)
  # require(ggeffects)
  # 
  # preds <- ggpredict(res.1,terms=c("SNR[all]"), pretty=FALSE, ci.lvl = 0.95)
  # p1<-plot(preds, rawdata=TRUE, dot.alpha=0.01) + #facet_wrap(~1,nrow = 1)+
  #   geom_hline(yintercept = .50)+
  #   labs(x = c("SNR (dB)"), y = "P(detecting a call)" , title = NULL)+ 
  #   theme(legend.position="bottom",legend.key.height = unit(0.25,"cm"), 
  #         legend.key.width = unit(0.25,"cm"),
  #         legend.text = element_text(size = 6))
  # # coord_cartesian(ylim=c(0,1),xlim=xlim)
  # # print(p1)
  # summary(res.1)
  # return(p1)
  require(ggplot2)
  require(marginaleffects)
  p1 <- plot_cap(res.1, condition='SNR', type='response',conf_level=0.95)+
    labs(x = c("SNR (dB)"), y = "P(detecting a call)" , title = NULL)+ 
      theme(legend.position="bottom",legend.key.height = unit(0.25,"cm"),
            legend.key.width = unit(0.25,"cm"),
            legend.text = element_text(size = 6))
    # coord_cartesian(ylim=c(0,1),xlim=xlim)
    # print(p1)
    summary(res.1)
    return(p1)
}

fitSNRbySeason <- function(SNRinfo, season=year, useGLM=TRUE, numKnots=3){
  # Initial glm:
  
  if (useGLM==TRUE){
    res.2=glm(Detected~(SNR*season), data = SNRinfo,
              family = binomial(link="logit")) # or: family quasibinomial 
  } else {
    res.2<-mgcv::gam(Detected~s(SNR,by=season, k=numKnots), data=SNRinfo,
               family=quasibinomial) # Test to include k=3
  }
  
  preds <- ggpredict(res.2,terms=c("SNR[all]","season[all]"), pretty=FALSE, ci.lvl = 0.95)
  p2<-plot(preds, rawdata=TRUE, dot.alpha=0.01) + #facet_wrap(~1,nrow = 1)+
    geom_hline(yintercept = .50)+
    labs(x = c("SNR (dB)"), y = "Probability of detecting a call" ,
         title = NULL, color=NULL)+ 
    theme(legend.position="bottom",legend.key.height = unit(0.25,"cm"), 
          legend.key.width = unit(0.25,"cm"),
          legend.text = element_text(size = 10))
  # coord_cartesian(ylim=c(0,1),xlim=xlim)
  # print(p2)
  # summary(res.2)# preds <- ggpredict(res.1,terms = c("SNR [-10:15]"))
  return(p2)
}