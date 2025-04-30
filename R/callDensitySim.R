#' Simulate animals calls using a uniform distribution time and space
#'
#' Generate a uniform distribution of locations in time and space. The space
#' component is in two dimensions (x, y), and the time component is in calendar
#' time (posix.ct).
#'
#' simCallLocation creates n calls within radius R of the receiver (hydrophone).
#' This is useful for testing callDensity package with a known call density
#' (uniform random in circular study area).
#'
#' @param n - Number of calls that will be created.
#'
#'   For reference 1 whale call every 10 s would be a maximum of 3.1536M calls
#'   In a continuous year of recording (i.e. 365*24*60*60/10 = 3.1536e6)
#' @param R - Radius of study area.
#' @param k - Number of sensors (presently always 1)
#' @param minDate - Starting date and time for the simulated data (POSIXct)
#' @param maxDate - Ending/latest date and time for the simulated data (POSIXct)
#'
#' @returns sim - data.frame containing a row for each simulated call and columns for the location and time of each call.
#' @export
#'
#' @examples
#' sim <- simCallLocation(n=1e6, R=1e6, minDate = Sys.time(), maxDate=Sys.time()-(365*24*60*60))
simCallLocation <- function(n=1e6, R=1e6, minDate=Sys.time(), maxDate=startDate+86400){

  # Call density, D_c is n/A
  A = pi*(R/1e3)^2 # Circular study area in km^2
  TrueCallDensity = n/A # calls/km^2

  k <- 1      # number of sensors

  r = R*sqrt(runif(n));
  th = 2*pi*runif(n);
  sim <- data.frame(x = r*cos(th), y = r*sin(th), groundTruth=1);

  duration_s <- as.numeric(difftime(maxDate,minDate,unit="sec"))
  sim$datetime <- minDate + sort(runif(n, 0, duration_s))

  sim$d <- sqrt(sim$x^2+sim$y^2)
  return(sim)

}

#' Simulate acoustic properties of calls
#'
#' @param sim - data.frame containing simulated call locations (from
#'   simCallLocations). Each location must have a distance.
#' @param SL - Source level distribution. SL is a data.frame with columns named
#'   'mean' and 'sd' that describe the mean and standard deviation of the
#'   distribution (normal in dB) from which source levels will be generated for
#'   each simulated call.
#' @param NL - Noise level distribution. NL is a data.frame with columns named
#'   'mean' and 'sd' that describe the mean and standard deviation of the
#'   distribution (normal in dB) of noise levels from which noise levels will be
#'   generated for each simulated call
#' @param TL
#'
#' @returns Simulation data.frame including columns that contain simulated
#'   source levels (SL) noise levels (noiseRMSdB), and transmission losses (TL).
#' @export
#'
#' @examples
simCallAcoustics <- function(sim,
                             SL=data.frame(mean=190, sd=4, sampleSize=350),
                             NL= data.frame(mean=84, sd = 4, sampleSize = n),
                             TL=function(r){20*log(r)}
                             ){
  n = dim(sim)[1]
  sim$SL <- rnorm(n, mean=SL$mean, sd=SL$sd)       # Realised SL for each call

  sim$noiseRMSdB <- rnorm(n, mean=NL$mean, sd=NL$sd)# Noise for each call

  # TODO: update this accept other functions
  sim$TL <- TL(sim$d)

  sim$signalRMSdB <- with(sim,SL-TL)
  sim$snr <- with(sim,SL-TL-noiseRMSdB)

  return(sim)
}

#' Simulate transmission loss (TL) radials following geometric (spherical)
#' spreading law.
#'
#' Create a TL data.frame with a column for range and subsequent columns for
#' transmission loss, TL for radial transects around a recorder. Transects will
#' have a total length of maxRange with increments between points along the
#' transect of rangeStep. The number of radial transects is specified as
#' numTransects such that the angular resolution of the transects will be
#' 360/numTransects.
#'
#' @param maxRange - Length of each radial transect
#' @param rangeStep - increment between points on each transect where TL
#' @param numTransects - Number of radial transects (i.e. columns in result)
#'
#' @returns TL -
#' @export
#'
#' @examples
simTLradials_20logR <- function (maxRange, rangeStep, numTransects){
  range_m = seq(from=5, to=maxRange, by=rangeStep)

  # Spherical spreading for a single transect
  tlTransectSpherical <- 20*log10(range_m)

  # Create additional transects here to better match real-world approaches
  tlTransects <- replicate(numTransects,tlTransectSpherical)

  #TL columns will be named tlXXX where XXX is the angle of that transect
  angleStep <- 360/numTransects
  angles = seq(from=0, to=360-angleStep, by=angleStep)
  colnames(tlTransects)<- paste0("tl",angles)

  TL <- data.frame(range_m)
  TL <- cbind(TL,tlTransects)

  return(TL)
}

#' Simulate a detector for a callDensity simulation.
#'
#' Given parameters of an acoustic detector and a callDensity simulation,
#' simulate the detection of true positive and false positive calls and append
#' this to the simulation.
#'
#' @param detParams - A data.frame of parameters for simulating the acoustic
#'   detector. The simulated detector will follow a logistic curve 'plogis' with
#'   the 'location', and 'scale' specified in the columns of detParams. The
#'   detector will also have a false discovery rate, 'c' specified in detParams
#'   as well as a mean and standard deviation of the distribution of the
#'   false-positive-to-noise-ratio, 'fpMean', and 'fpSD' respectively.
#' @param sim - a simulation data.frame created by simCallLocation and
#'   simCallAcoustics
#'
#' @returns sim - simulation data.frame with additional columns for the detector
#'   status and rows for false positives. Columns added by this function
#'   include:
#'      detect_table - logical indicating which calls were detected
#'      groundTruth - logical indicating whether the call is a true positive (1)
#'        or false positive (0).
#'      group - factor with three levels indicating whether the call is a true
#'        positive, false positive, or false negative.
#' @export
#'
#' @examples
simulateDetector <- function(detParams,sim){

  # Probability of detection for each simulated call
  sim$p_det = plogis(sim$snr,
                     location=detParams$location,
                     scale=detParams$scale)

  # Apply Bernoulli trial to SNRs to see which were detected
  sim$detect_table <- as.logical(rbinom(dim(sim)[1],size=1,prob=sim$p_det))

  # Calculate number of false positives given the number of true positive
  # detections for the detector and the specified false discovery rate, c, for
  # that detector.
  n_tp_det1 <- sum(sim$detect_table)
  n_fp_det1 <- n_tp_det1/(1-detParams$c)-n_tp_det1

  # Store false positives in data.frame with same format as sim
  fp <- data.frame(matrix(ncol=length(sim),nrow=n_fp_det1))
  colnames(fp) <- colnames(sim)

  # groundTruth=0 means that these are false positives (handy to keep in mind for
  # when we include these in sim)
  fp$groundTruth <- FALSE
  fp$detect_table <- TRUE


  duration_s<-as.numeric(difftime(max(sim$datetime),min(sim$datetime),
                                  unit="sec"))

  # Generate the right number of false positives uniformly over same time period
  # as true positives.
  fp$datetime <- minDate + sort(runif(n_fp_det1, 0, duration_s))

  # False positives don't have a location, SL, or TL, but do have NL and SNR
  # Use same distribution of NL as true positives
  fp$noiseRMSdB <- rnorm(n_fp_det1, mean=NL$mean, sd=NL$sd)

  # SNR distribution of false positives broadly matches that of Casey2019 human
  # analysts
  fp$snr <- rnorm(n_fp_det1, mean=detParams$fpMean, sd=detParams$fpSD)

  fp$signalRMSdB <- fp$noiseRMSdB + fp$snr

  # Combine false positives into the simulation
  sim <- rbind(sim,fp)

  # We've generated the right number of false positive detections now and merged
  # these into our simulation.

  # The false positives are missing locations, distances, SL, and p_dets. But
  # not sure that this actually matters in any meaningful way.

  # Add factor column to track true positive, false positivs & missed detections
  sim$group <- factor(ifelse(sim$groundTruth,
                             ifelse(sim$detect_table,
                                    "TruePositive","FalseNegative"),
                             "FalsePositive"),
                      levels=c("FalsePositive","TruePositive","FalseNegative")
  )

  return(sim)

}


