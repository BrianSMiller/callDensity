# tests/testthat/helper-fixtures.R
#
# Shared fixtures loaded automatically by testthat before any test file runs.
#
# All simulation parameters match the callDensity.Rmd vignette exactly.
# n is reduced from 1e6 to 1e4 for speed; everything else is unchanged.

# ---------------------------------------------------------------------------
# Full simulation pipeline (matches vignette)
# ---------------------------------------------------------------------------

make_sim_data <- function(n = 1e4, seed = 42) {
  set.seed(seed)

  R        <- 1e6   # radius in m
  minDate  <- as.POSIXct("2025-01-01")
  maxDate  <- as.POSIXct("2026-01-01")

  SL <- data.frame(mean = 190, sd = 4, sampleSize = 350)
  NL <- data.frame(mean = 84,  sd = 4, sampleSize = n)

  det1params <- data.frame(
    location = 3,
    scale    = 2,
    func     = "plogis",
    c        = 0.1,
    fpMean   = 0,
    fpSD     = 4
  )

  tlFunc <- function(r) 20 * log10(r)

  sim <- simCallLocation(n = n, R = R, minDate = minDate, maxDate = maxDate)
  sim <- simCallAcoustics(sim, SL, NL, TL = tlFunc)
  sim <- simulateDetector(detParams = det1params, sim)
  sim <- sim[order(sim$datetime), ]

  list(sim = sim, SL = SL, NL = NL, R = R,
       minDate = minDate, maxDate = maxDate)
}

# ---------------------------------------------------------------------------
# SNRinfo data frame for detection function fitting (matches vignette)
# ---------------------------------------------------------------------------

make_snr_data <- function(n_subsample = 1e3, seed = 42) {
  s   <- make_sim_data(seed = seed)
  sim <- s$sim
  set.seed(seed)
  subsample <- sim[sample(nrow(sim), n_subsample), ]
  SNRinfo <- with(subsample,
    data.frame(
      Detected = detect_table,
      CallRL   = SL - TL,
      NoiseRL  = noiseRMSdB,
      SNR      = snr,
      t        = datetime,
      season   = callDensity::time2season(datetime)
    )
  )
  SNRinfo
}

# ---------------------------------------------------------------------------
# TL lookup table and SL/NL for pDetInArea (matches vignette)
# ---------------------------------------------------------------------------

make_tl_lookup <- function(R = 1e6, rangeStep = 100, numTransects = 4) {
  simTLradials_20logR(maxRange = R, rangeStep = rangeStep,
                      numTransects = numTransects)
}

make_sl <- function() data.frame(mean = 190, sd = 4, sampleSize = 350)

make_nl <- function(n_subsample = 1e3, seed = 42) {
  SNRinfo <- make_snr_data(n_subsample = n_subsample, seed = seed)
  SNRinfo |>
    dplyr::summarise(
      mean       = mean(NoiseRL, na.rm = TRUE),
      sd         = sd(NoiseRL,   na.rm = TRUE),
      sampleSize = dplyr::n() - sum(is.na(NoiseRL))
    )
}

# ---------------------------------------------------------------------------
# Two-observer data for vglm tests
# ---------------------------------------------------------------------------
# Requires independent partial detection by each observer, with neither
# achieving perfect separation. Constructed directly from SNRinfo so the
# SNR distribution is realistic; detection flags are assigned by independent
# Bernoulli draws at moderate probability.

make_two_observer_data <- function(n_subsample = 1e3, seed = 42) {
  d <- make_snr_data(n_subsample = n_subsample, seed = seed)
  set.seed(seed + 1)
  p <- plogis(0.3 * (d$SNR - median(d$SNR, na.rm = TRUE)))
  d$detect_observer1 <- rbinom(nrow(d), 1, p * 0.8)
  d$detect_observer2 <- rbinom(nrow(d), 1, p * 0.7)
  # vglm requires at least one detection by each observer
  d <- d[d$detect_observer1 + d$detect_observer2 > 0, ]
  d
}

# --- Fixtures for pDetGivenNL and nlFromDetections ----------------------------
# Deliberately independent of make_sim_data(). The point of these tests is to
# check pDetGivenNL against a simulation that shares no code with it, so they
# must not share a fixture with it either.

#' Spherical spreading TL table, r = 0 and non-finite rows removed.
tlSpherical <- function(R = 1e6, rangeStep = 1000, numTransects = 4) {
  TL <- simTLradials_20logR(maxRange = R, rangeStep = rangeStep,
                            numTransects = numTransects)
  TL[is.finite(rowSums(TL)) & TL[[1]] > 0, ]
}

#' Antarctic blue whale source levels, as used in the vignettes.
testSL <- data.frame(mean = 190, sd = 4, sampleSize = 350)

#' A plain logistic detector. Not a detFun object, so pDetGivenNL's
#' is.function() branch is exercised.
testDetector <- function(snr) plogis(snr, location = 1, scale = 2)

#' Place calls on the water, detect them, return the noise at the detections.
#'
#' Brute force. Every call gets a range, a source level and a noise level, and
#' detection is a Bernoulli trial. No area weighting, no quadrature, no binning.
simulateDetectedNoise <- function(n = 1e6, R = 1e6, nlMean = 84, nlSd = 4,
                                  detector = testDetector, seed = 1) {
  set.seed(seed)
  r   <- R * sqrt(runif(n))                       # uniform in two dimensions
  sl  <- rnorm(n, testSL$mean, testSL$sd)
  nl  <- rnorm(n, nlMean, nlSd)
  snr <- sl - 20 * log10(r) - nl
  det <- runif(n) < detector(snr)
  data.frame(NoiseRL = nl[det], SNR = snr[det])
}


# --- Capture history table for cde tests --------------------------------------
# Two detectors on the same calls, merged the way the vignettes do it. This is
# the only fixture that produces something cde will accept: cde calls
# falseDiscoveryRate() and capHist2snrInfo() on its first two lines, and both
# need the detect_tableN / groundTruthN / consolidated-column structure that
# simsTocaptureHistoryTable builds.
#
# No time subsampling. The vignettes subsample 1 hour in 41 to mimic a real
# annotation effort, which at this n would leave a couple of dozen detections.
# Here the point is to exercise cde, not to be realistic about labour.

make_capture_history <- function(n = 5e4, seed = 42,
                                 det1location = 1, det1scale = 2,
                                 det2location = 2, det2scale = 4,
                                 fdr = 0.1,
                                 nlMean = 84, nlSd = 4) {
  set.seed(seed)

  R       <- 1e6
  minDate <- as.POSIXct("2025-01-01")
  maxDate <- as.POSIXct("2026-01-01")
  Time    <- as.numeric(difftime(maxDate, minDate, unit = "days")) / 365

  SL <- data.frame(mean = 190, sd = 4, sampleSize = 350)
  NL <- data.frame(mean = nlMean, sd = nlSd, sampleSize = n)

  det1params <- data.frame(location = det1location, scale = det1scale,
                           func = "plogis", c = fdr, fpMean = 0, fpSD = 2)
  det2params <- data.frame(location = det2location, scale = det2scale,
                           func = "plogis", c = fdr, fpMean = 0, fpSD = 2)

  tlFunc <- function(r) 20 * log10(r)

  sim     <- simCallLocation(n = n, R = R, minDate = minDate, maxDate = maxDate)
  sim     <- simCallAcoustics(sim, SL, NL, TL = tlFunc)
  simDet1 <- simulateDetector(det1params, sim)
  simDet2 <- simulateDetector(det2params, sim)

  A <- studyArea(R / 1e3)

  list(capHistTab = simsTocaptureHistoryTable(simDet1, simDet2),
       SL     = SL,
       NL     = NL,
       R      = R,
       A      = A,
       Time   = Time,
       Nc     = sum(simDet2$detect_table),
       truePa = mean(simDet2$p_det, na.rm = TRUE),
       trueDc = n / (A * Time))
}
