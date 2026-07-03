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
