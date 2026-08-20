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
#' @param minDate - Starting date and time for the simulated data (POSIXct)
#' @param maxDate - Ending/latest date and time for the simulated data (POSIXct)
#'
#' @returns sim - data.frame containing a row for each simulated call and columns for the location and time of each call.
#' @export
#'
simCallLocation <- function(n=1e6, R=1e6, minDate=Sys.time(), maxDate=minDate+86400){

  # Call density, D_c is n/A
  A = pi*(R/1e3)^2 # Circular study area in km^2
  TrueCallDensity = n/A # calls/km^2

  k <- 1      # number of sensors

  r = R*sqrt(runif(n));
  th = 2*pi*runif(n);
  sim <- data.frame(x = r*cos(th), y = r*sin(th), groundTruth=1);

  duration_s <- as.numeric(difftime(maxDate,minDate,units="sec"))
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
#' @param TL - a function that takes a single argument, vector r, and returns
#'   transmission losses for the ranges in that vector. The default function is
#'   spherical spreading: Default = function(r){20*log(r)}
#'
#' @returns Simulation data.frame including columns that contain simulated
#'   source levels (SL) noise levels (noiseRMSdB), and transmission losses (TL).
#' @export
simCallAcoustics <- function(sim,
                             SL=data.frame(mean=190, sd=4, sampleSize=350),
                             NL= data.frame(mean=84, sd = 4, sampleSize = n),
                             TL=function(r){20*log10(r)}
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
  n_tp <- sum(sim$detect_table)
  n_fp <- as.integer(n_tp/(1-detParams$c)-n_tp)

  # A false discovery rate of zero is legitimate and yields no false positives.
  # Skip the block rather than build a zero-row data.frame, which cannot be
  # assigned into.
  if (n_fp > 0){

    # Store false positives in data.frame with same format as sim
    fp <- data.frame(matrix(ncol=length(sim),nrow=n_fp))
    colnames(fp) <- colnames(sim)

    # groundTruth=0 means that these are false positives (handy to keep in mind for
    # when we include these in sim)
    fp$groundTruth <- FALSE
    fp$detect_table <- TRUE


    duration_s<-as.numeric(difftime(max(sim$datetime),min(sim$datetime),
                                    units="sec"))

    # Generate the right number of false positives uniformly over same time period
    # as true positives.
    fp$datetime <- min(sim$datetime) + sort(runif(n_fp, 0, duration_s))

    # False positives don't have a location, SL, or TL, but do have NL and SNR
    # Use same distribution of NL as true positives
    fp$noiseRMSdB <- rnorm(n_fp,
                           mean = mean(sim$noiseRMSdB, na.rm=TRUE),
                           sd   = sd(sim$noiseRMSdB, na.rm=TRUE))

    # SNR distribution of false positives broadly matches that of Casey2019 human
    # analysts
    fp$snr <- rnorm(n_fp, mean=detParams$fpMean, sd=detParams$fpSD)

    fp$signalRMSdB <- fp$noiseRMSdB + fp$snr

    # Combine false positives into the simulation
    sim <- rbind(sim,fp)
  }
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

#' Subsample from a simulation at evenly spaced time intervals.
#'
#' @param sim - simulation data.frame containing column a time column named
#'   'datetime'
#' @param minDate - posixCT indicating the start of the first subsample
#' @param maxDate - posixCT indicating the start of the last subsample
#' @param interval - difftime interval between subsamples (e.g. '41 hour')
#' @param duration - numeric indicating the duration (in s) of each subsample
#'
#' @returns subsampledSim - a simulation dataframe containing only rows from the
#'   input that fall within the subsample
#' @export
#'
subsampleSimInTime <- function(sim,
                               minDate=min(sim$datetime),
                               maxDate=max(sim$datetime),
                               interval='41 hour',
                               duration=3600){
  sim$subset <- 0
  subStart <- seq(from=minDate,to=maxDate,by=interval)
  subEnd <- subStart+duration

  for (i in 1:length(subStart)){
    sim$subset <- sim$subset |
      (sim$datetime >= subStart[i] & sim$datetime <= subEnd[i])
  }
  sim<- subset(sim,sim$subset, select=-c(subset))
  return(sim)
}

#' Merge N simulated detectors into a matchbox-native capture history table
#'
#' @description
#' Generalizes the original two-detector merge to any number N >= 2 of
#' simulated detectors, and aligns the output columns with matchbox's own
#' canonical convention (\code{key}, \code{t0}/\code{tEnd},
#' \code{detect_observerK}, \code{<col>_observerK} for every other column,
#' suffixed and \code{NA} where that detector missed the event) -- rather
#' than the ad hoc consolidated (\code{SNR}, \code{signalRMSdB}, ...)
#' columns this function used to build to match a now-superseded MATLAB
#' format. Downstream consolidation across observers (e.g. averaging SNR)
#' is \code{\link{chtToSNRinfo}}'s job, via its own \code{signalCol}/
#' \code{noiseCol} vector-averaging support -- not this function's.
#'
#' Detectors are matched by exact \code{datetime} equality, not the
#' temporal-overlap clustering real matchbox tables use for actual bounded
#' detections. This is deliberate, not a shortcut: every detector here runs
#' against the same underlying simulated call population (via
#' \code{\link{simulateDetector}}), so a true positive's \code{datetime} is
#' bit-identical across every detector that examined it, while each
#' detector's own false positives get independently random \code{datetime}s
#' with a negligible chance of colliding with anyone else's. Real acoustic
#' detections need fuzzy overlap-based matching precisely because their
#' bounds are imprecise; these simulated point-in-time detections don't
#' have that problem to begin with.
#'
#' \code{fLow}/\code{fHigh} (matchbox's frequency-envelope columns) are
#' deliberately omitted: this simulation pipeline never models frequency at
#' all, so there is nothing real to put there. The result is a plain
#' data.frame, so add them yourself later if a specific analysis needs them.
#'
#' @param ... Two or more data.frames, each the output of
#'   \code{\link{simulateDetector}} for one simulated detector, sharing a
#'   \code{datetime} column.
#' @param observerSuffix Character vector, one name per detector in
#'   \code{...}, used to suffix every one of that detector's own columns
#'   (\code{detect_<suffix>}, \code{snr_<suffix>}, \code{groundTruth_<suffix>},
#'   etc.) in the output. Default \code{paste0('observer', seq_along(list(...)))}
#'   -- matchbox-native, matching \code{\link{chtToSNRinfo}}'s own default
#'   \code{observers} naming, so a simulated capture history table needs no
#'   renaming before use. Pass e.g. \code{c('table1','table2')} for the
#'   older \code{detect_table1}/\code{detect_table2} shape some existing
#'   tests and fixtures still key on.
#'
#' @returns A data.frame with one row per matched event: \code{key}
#'   (sequential event id), \code{t0}/\code{tEnd} (both the shared event
#'   time, as a genuine MATLAB datenum via \code{\link{Rdate2mat}} -- matching
#'   real matchbox's own \code{t0}/\code{tEnd} encoding, so
#'   \code{chtToSNRinfo()}'s default \code{timeCol='t0'} and \code{cde()}'s
#'   \code{timeCol=NULL} auto-detection both work unmodified on this table),
#'   \code{datetime} (the original, exact POSIXct, kept for anyone who wants
#'   a directly readable time rather than a datenum), \code{detect_<suffix>}
#'   (logical, per detector, \code{FALSE} where that detector examined the
#'   event and missed it), \code{groundTruth} (coalesced across detectors --
#'   see Details -- present only when the inputs have their own
#'   \code{groundTruth} column, as \code{\link{simulateDetector}}'s output
#'   does), and every other column from each detector's own output, suffixed
#'   \code{_<suffix>} and \code{NA} wherever that detector has no row for
#'   the event at all (this includes a per-observer
#'   \code{groundTruth_<suffix>} alongside the coalesced \code{groundTruth}).
#'
#' @export
simsTocaptureHistoryTable <- function(..., observerSuffix = NULL) {

  detectors <- list(...)
  n <- length(detectors)
  if (n < 2) {
    stop("simsTocaptureHistoryTable: need at least two detectors.", call. = FALSE)
  }

  if (is.null(observerSuffix)) {
    observerSuffix <- paste0('observer', seq_len(n))
  }
  if (length(observerSuffix) != n) {
    stop("simsTocaptureHistoryTable: observerSuffix must have one entry per ",
         "detector (", n, " detectors, ", length(observerSuffix),
         " suffixes given).", call. = FALSE)
  }

  # Suffix every column except the join key before merging, so an N-way
  # Reduce(merge, ...) needs no per-step suffixes= (merge() only supports a
  # suffix pair for exactly two tables at a time, which doesn't generalize).
  suffixed <- lapply(seq_len(n), function(i) {
    d <- detectors[[i]]
    isKey <- names(d) == 'datetime'
    names(d)[!isKey] <- paste0(names(d)[!isKey], '_', observerSuffix[i])
    d
  })

  cht <- Reduce(function(a, b) merge(a, b, by = 'datetime', all = TRUE), suffixed)

  # detect_<suffix>: always clean TRUE/FALSE, matching real matchbox's own
  # contract for this column specifically (every other column is left as
  # NA where missing, per matchbox's "NaN where missed" convention for
  # those -- only the detect flag itself is guaranteed non-NA).
  for (i in seq_len(n)) {
    rawDetectCol <- paste0('detect_table_', observerSuffix[i])
    newDetectCol <- paste0('detect_', observerSuffix[i])
    names(cht)[names(cht) == rawDetectCol] <- newDetectCol
    cht[[newDetectCol]][is.na(cht[[newDetectCol]])] <- FALSE
    cht[[newDetectCol]] <- as.logical(cht[[newDetectCol]])
  }

  # groundTruth is unique among the per-detector columns: it's a property
  # of the event itself (is this a real call at all), not a per-observer
  # measurement like snr/signalRMSdB, or a per-observer outcome like group
  # (a detector's own TruePositive/FalseNegative label, which legitimately
  # does differ between detectors for the same event and is correctly left
  # suffixed-only). Every contributing detector's own copy of groundTruth
  # agrees whenever more than one is present, since they're all evaluating
  # the same underlying simulated event -- so unlike every other column,
  # it's both safe and useful to coalesce into one event-level column here,
  # rather than leaving NA on any event only one detector's own table had a
  # row for (that detector's own fabricated false positive). Kept alongside
  # the per-observer groundTruth_<suffix> columns, not instead of them.
  gtCols <- paste0('groundTruth_', observerSuffix)
  if (all(gtCols %in% names(cht))) {
    cht$groundTruth <- Reduce(function(a, b) ifelse(is.na(a), b, a), cht[gtCols])
  }

  cht <- cht[order(cht$datetime), ]
  cht$key  <- seq_len(nrow(cht))
  # Genuine MATLAB datenum, matching real matchbox's actual t0/tEnd encoding
  # -- not a same-named POSIXct, which chtToSNRinfo()'s default timeCol='t0'
  # and cde()'s timeCol=NULL auto-detection would both misinterpret (they
  # unconditionally run mat2Rdate() on anything named t0). datetime (POSIXct,
  # exact) is left untouched above this point for the merge/matching logic,
  # and remains in the output for anyone who wants a plain readable time.
  cht$t0   <- Rdate2mat(cht$datetime)
  cht$tEnd <- Rdate2mat(cht$datetime)

  rownames(cht) <- NULL
  cht
}

#' SNR histogram of true & false positives, and false negatives for callDensity simulation
#'
#' @param sim - callDensity simulation data.frame
#'
#' @returns - ggplot object of type geom_histogram
#' @importFrom ggplot2 ggplot geom_histogram labs theme
#' @importFrom grid unit
#' @export
#'
plotDetectionDistribution <- function(sim){
  # SNR Distribution (includes false positives)
  ggplot(data=sim, aes(x=snr, group=group, fill=group) )+
    geom_histogram(binwidth=1, position="identity",alpha=0.3)+
    # facet_wrap(group~.,nrow=3, scales='free_y')+
    labs(fill='')+
    theme(legend.text = element_text(size=8))+
    # guides(fill="none")+
    theme_minimal()+
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.9),
          legend.key.size = unit(0.25, 'cm'))
    # scale_fill_viridis_d(direction = -1, option = "plasma")
    # scale_color_brewer(palette = "RdYlBu",type="qual",aesthetics = c('fill'))
}

#' Plot the spatial detection density for a callDensity simulation
#'
#' @param sim - callDensity simulation
#'
#' @returns ggplot object of type geom_bin_2d
#' @importFrom ggplot2 ggplot geom_bin_2d coord_equal
#' @importFrom grid unit
#' @export
#'
plotSpatialDetections <- function(sim){

  # Spatial distribution (excludes false positives)
  ggplot2::ggplot(data=sim, aes(x=x/1e3, y=y/1e3, weight=detect_table) )+
    geom_bin_2d(alpha=1,binwidth=c(10,10))+
    coord_equal()+
    xlab("X location (km)")+
    ylab("Y location (km)")
}


