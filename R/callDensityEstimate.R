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
#'     Usually created by calling function defaultOutputFileNames
#' @param season TimeCode specifying month, season, or year for outputs
#' @returns data.frame containing call density inputs, results, and CVs
#'
#' @importFrom stats rnorm sd vcov
#' @importFrom utils read.csv write.table
#' @importFrom magrittr %>%
#' @export
#'
cde <- function (p,season){

  # Require that parameters have already been defined and are in data.frame, p
  if (!exists("p")){
    warning(c('Required parameter variable, p, not found.
  Please define all call density parameters in data.frame, p, before running
  this script.\n'))
  }

  if (!exists("season")){
    warning(c(  'Required parameter variable, season, not found.
  Please define a variable called season, and assign it a value of
  summer, autumn, winter, spring, or year.\n\n
  Now assuming season<-year'))
    season <-'year'
  }

  # Check inputs and create outputs
  # Store all parameters for call-density estimation in data frame called 'p'
  #
  ## ----------------------------------------------------------------------------------------------------

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

  #
  # ### $N_c$ (number of detections).
  #
  # $T$ and $N_c$ were calculated based, respectively, in the raw data duration and number of automatic detections correspondent to the full dataset.
  #
  ## ----------------------------------------------------------------------------------------------------

  det <- read.csv(p$fullYearDetectionCsv,sep=',')
  det$t <- as.POSIXct(det$UTC, tz='UTC')
  det$season <- time2season(det$t)
  det$month <- time2monthCode(det$t)

  det<- subsetByTimeCode(det,det$t,season)

  Nc <- dim(det)[1]
  Nc

  #
  #  ### $T$ (time of deployment)
  #
  #  Sum the duration of the audio that has been analysed
  #
  ## ----------------------------------------------------------------------------------------------------

  wavInfo <- read.csv(p$fullYearEffortFile) # File generated from Matlab wavFolderInfo
  # wavInfo$season <- time2season( mat2Rdate(wavInfo$startDate) )
  wavInfo <- subsetByTimeCode(wavInfo,mat2Rdate(wavInfo$startDate),season)

  T <- sum(wavInfo$duration)/3600 # in hours
  T

  #  ### Estimate $c$ (False discovery rate)
  #
  #  False discovery rate, c, is defined as: FP/(FP+TP), where FP is the number
  #  of false-positives and TP is the number of true positives
  #  (https://en.wikipedia.org/wiki/Precision_and_recall; TODO: find a primary
  #  literature citation to use instead of citing wikipedia).
  #
  #  Here we estimate false-positive rate of the automated detector from the
  #  manually annotated dataset. The data-file for this is the capture-history
  #  table that contains reconciled annotated, and automated detections.
  #
  #  First, we count the number of FP that the AI produced from the capture
  #  history table. These are the sum of the rows where detect_table2==T &
  #  detect_table1==F. Then we calculate the number of true positives, the sum
  #  of the rows where detect_table2==T & detect_table1==T.
  #
  #  Alternatively, we could use the double observer mark-recapture to estimate
  #  the total number of calls in the dataset (see Miller et al 2022 - AI bests
  #  human observer - RSEC). This would require incorporating the Huggins/RMark
  #  script that can be found in the supplement to the Miller et al 2022 paper.
  ##
  #----------------------------------------------------------------------------------------------------

  effort <- readxl::read_excel(p$effortFile);
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

  #  ### $CV_c$ (Cochran approx.);
  #
  #  First calculate using annotated library
  #
  ## ----------------------------------------------------------------------------------------------------

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

  #  The above calculation uses the annotated library, and requires that you
  #  have defined effort periods for which you've measured the false discovery
  #  rate.
  #
  #  Below is an alternative measure of false-discovery-rate where the analyst
  #  has inspected every nth detection to determine if it is true or false
  #  positive.
  ## ----------------------------------------------------------------------------------------------------

  if (p$useSeparateFPdata){
    ffp <- readxl::read_xlsx(p$manualFPfileName,sheet = 'conference')

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
  }
  c

  #  Now calculate CV of c where the analyst has inspected every 50th detection.
  #
  #  Here we group the detections by day within each season in order to estimate
  #  CVs.
  ## ----------------------------------------------------------------------------------------------------
  if (p$useSeparateFPdata){

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
  }

  #  ### $P_a$ (Overall probability of detection)
  #
  ## ---- echo=TRUE, result='hide', message=F, warning=F, cache=TRUE-------------------------------------

  SNRinfo <- capHist2snrInfo(p$capHistFile,season)
  NL <- SNRinfo %>% dplyr::summarise(mean=mean(NoiseRL), sd=sd(NoiseRL),
                                     sampleSize=dplyr::n())
  SL <- data.frame(mean=p$SLmean, sd=p$slStd, sampleSize=p$SLsamplesize);
  TL<-utils::read.csv(p$tlFile)
  # Initial GLM or GAM:
  snr.det.fun <- fitSNRdetectionFunc(SNRinfo,
                                     useGLM=p$useGLM, useSCAM=p$useSCAM, p$numKnots)

  pDetInArea(snr.det.fun, SL, TL,  NL, # Sonar equation inputs
             p$transectFile, p$simResultsFile, p$paFile, # file output names
             useGLM=p$useGLM, useSCAM=p$useSCAM, numKnots=p$numKnots, # detector SNR curves
             output.resolution.m=p$output.resolution.m, outerloop=p$outerloop)


  ## ----------------------------------------------------------------------------------------------------
  pa.all.transects <- read.csv(p$paFile, header = TRUE, sep = ' ');
  no.transects <- dim(pa.all.transects)[1]-1;
  pa <- pa.all.transects[no.transects+1,1]
  pa

  #  ### $CV_{P_a}$ (see Harris, 2012 for details);
  #
  #  Approach is to calculate CV using mean and SD from each radial?
  #
  ## ----------------------------------------------------------------------------------------------------
  # Coef. of variation - Pa

  names(pa.all.transects)<-c('Mean','SD')


  # Radials from 0-315 degrees are in rows 1:8
  y=pa.all.transects[1:no.transects,]


  ## IMPORTANT: to check if the var, SE and CV calculation is correct ##
  ## Per year:

  varp.pa<-((sd(y$Mean)/sqrt(no.transects))^2)+((sum(y$SD)/no.transects)^2)
  sep.pa<-sqrt(varp.pa)
  CV.pa<-sep.pa/mean(y$Mean)
  CV.pa

  #  ### $CV_{N_c}$ (see Harris, 2012 for details);
  #
  ## ---- warning=FALSE----------------------------------------------------------------------------------
  varn=Nc*(pa)*(1-pa)
  varn

  # Confirm how to calculate SD and CV (?)

  SDn=sqrt(varn)
  SDn
  CV.Nc=SDn/Nc
  CV.Nc

  #  ### Finally, estimate Dc (Density of calls) [calls $h^{-1} km^{-2}$]
  #
  ## ----------------------------------------------------------------------------------------------------

  Dc <- (Nc * (1-c) )/( p$k * pi*p$w^2 * pa * T )
  Dc     # per h per km^2
  Dc*1e3 # per h per 1000 km^2?

  #  ### CV_Total (Delta Method)
  #
  ## ----------------------------------------------------------------------------------------------------

  CV.Dc<-sqrt((CV.Nc^2)+(CV.pa^2)+(CV.c^2))
  CV.Dc


  #  ## Collate into data frame and write to csv file
  ## ----------------------------------------------------------------------------------------------------

  result <- data.frame(season, p$siteCode,Nc,c,p$k,T,p$w,pa,
                       p$SLmean,p$SLsd,NL$mean,NL$sd,p$useGLM,
                       CV.Nc,CV.c,CV.pa,Dc,CV.Dc)
  names(result) <- c('season','siteCode','Nc','c','k','T','w','pa',
                     'SLmean','SLsd','NLmean','NLsd','useGLM',
                     'CV.Nc','CV.c','CV.pa','Dc','CV.Dc')

  utils::write.csv(result, file = p$densityResultsFile, row.names=F)
  return(result)
}



