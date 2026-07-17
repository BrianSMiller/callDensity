# Tests for SNR truncation.
#
# The claim under test is about the estimand, not the arithmetic:
#
#   Zeroing p below theta makes p_a = E[p * 1(SNR >= theta)]. Paired with an Nc
#   that counts only detections above theta, the q(theta) factor cancels and Dc
#   is the density of ALL calls.
#
#   Setting p to NA instead renormalises, giving p_a = E[p | SNR >= theta], and
#   Dc is then the density of above-threshold calls only, which is q(theta)
#   times smaller and is not what anybody wants.
#
# Following the pDetGivenNL pattern: apply the threshold by subsetting in a
# simulation that shares no code with the integral, and by zeroing p inside the
# integral, then check they agree. A test built only from package internals
# could confirm the code does what it does and still be wrong.

# --- Fixture: a fitted detector -----------------------------------------------
# pDetInArea bootstraps the coefficient vector, so it needs a fitted model and
# will not take a plain closure. Fit a glm on a large sample drawn from
# testDetector: with a logit link and a logistic truth the model is correctly
# specified, so the fit is the truth to within a rounding error, and the
# bootstrap draws are negligible.
#
# The brute force side of each test below evaluates THIS OBJECT rather than the
# closure, so any residual fitting error cancels and the only thing left under
# test is the truncation bookkeeping.

# The predictor MUST be named SNR: pDetInArea builds its newdata as
# data.frame(SNR = ...). A model fitted on a differently named column does not
# error, it silently returns fitted values of the wrong length.
fitLogisticDetector <- function(detector, n = 2e5, seed = 99) {
  set.seed(seed)
  SNR <- runif(n, -40, 40)
  detected <- runif(n) < detector(SNR)
  stats::glm(detected ~ SNR, family = stats::binomial,
             data = data.frame(SNR = SNR, detected = detected))
}

fitTestDetector <- function(...) fitLogisticDetector(testDetector, ...)

predFitted <- function(fit, snr) {
  as.vector(stats::predict(fit, newdata = data.frame(SNR = snr),
                           type = "response"))
}

# --- 1. Brute force: does zeroing give E[p * 1(SNR >= theta)]? ---------------

test_that("pDetInArea with SNR truncation agrees with a brute force simulation", {
  skip_on_cran()

  R     <- 1e6
  TL    <- tlSpherical(R, rangeStep = 1000)
  SL    <- testSL
  NL    <- data.frame(mean = 84, sd = 4, sampleSize = 1e3)
  theta <- 5

  # Brute force. Calls on the water, dice rolled, threshold applied by
  # subsetting. No quadrature, no binning, no package code.
  set.seed(42)
  n    <- 5e5
  r    <- R * sqrt(runif(n))                  # uniform in two dimensions
  sl   <- rnorm(n, SL$mean, SL$sd)
  nl   <- rnorm(n, NL$mean, NL$sd)
  snr  <- sl - 20 * log10(r) - nl
  fit  <- fitTestDetector()
  pAll <- predFitted(fit, snr)

  # The estimand: mean over ALL calls of p, with below-threshold calls
  # contributing zero. Not the mean over the calls that survive the threshold.
  paSim <- mean(pAll * (snr >= theta))

  pa <- pDetInArea(fit, SL, TL, NL,
                   outerloop = 20,
                   snrTruncationThreshold = theta)$perTransectMeanSD

  paInt <- pa[nrow(pa), 1]

  expect_equal(paInt, paSim, tolerance = 0.03)
})

test_that("SNR truncation zeroes rather than renormalises", {
  # The distinction that decides the estimand. If truncated cells were NA they
  # would drop out of the weighted mean's denominator, and p_a would come back
  # close to the conditional mean E[p | SNR >= theta], which is much larger.
  skip_on_cran()

  R     <- 1e6
  TL    <- tlSpherical(R, rangeStep = 1000)
  SL    <- testSL
  NL    <- data.frame(mean = 84, sd = 4, sampleSize = 1e3)
  theta <- 5

  set.seed(7)
  n   <- 3e5
  r   <- R * sqrt(runif(n))
  sl  <- rnorm(n, SL$mean, SL$sd)
  nl  <- rnorm(n, NL$mean, NL$sd)
  snr <- sl - 20 * log10(r) - nl
  fit <- fitTestDetector()
  p   <- predFitted(fit, snr)
  above <- snr >= theta

  paUnconditional <- mean(p * above)   # what zeroing gives
  paConditional   <- mean(p[above])    # what NA would give

  # The two must be far apart, or this test proves nothing.
  expect_gt(paConditional / paUnconditional, 3)

  pa <- pDetInArea(fit, SL, TL, NL,
                   outerloop = 20,
                   snrTruncationThreshold = theta)$perTransectMeanSD
  paInt <- pa[nrow(pa), 1]

  expect_equal(paInt, paUnconditional, tolerance = 0.05)
  expect_false(isTRUE(all.equal(paInt, paConditional, tolerance = 0.2)))
})

test_that("truncating at -Inf changes nothing", {
  skip_on_cran()

  R  <- 1e6
  TL <- tlSpherical(R, rangeStep = 2000)
  SL <- testSL
  NL <- data.frame(mean = 84, sd = 4, sampleSize = 1e3)
  fit <- fitTestDetector()

  set.seed(1)
  paNone <- pDetInArea(fit, SL, TL, NL, outerloop = 10)$perTransectMeanSD
  set.seed(1)
  paInf  <- pDetInArea(fit, SL, TL, NL, outerloop = 10,
                       snrTruncationThreshold = -Inf)$perTransectMeanSD

  expect_equal(paNone[nrow(paNone), 1], paInf[nrow(paInf), 1])
})

test_that("distance truncation still renormalises", {
  # Distance truncation must NOT be zeroed: A is redefined as pi*w^2 to match,
  # so those cells leave the estimand entirely. If the split in pDetInArea ever
  # gets applied to the wrong truncation, p_a here collapses towards zero.
  skip_on_cran()

  R  <- 1e6
  TL <- tlSpherical(R, rangeStep = 1000)
  SL <- testSL
  NL <- data.frame(mean = 84, sd = 4, sampleSize = 1e3)
  w  <- 3e5
  fit <- fitTestDetector()

  set.seed(3)
  pa <- pDetInArea(fit, SL, TL, NL, outerloop = 20,
                   truncationDistance = w)$perTransectMeanSD
  paInt <- pa[nrow(pa), 1]

  # Brute force within w only, no threshold.
  set.seed(3)
  n   <- 3e5
  r   <- w * sqrt(runif(n))
  sl  <- rnorm(n, SL$mean, SL$sd)
  nl  <- rnorm(n, NL$mean, NL$sd)
  snr <- sl - 20 * log10(r) - nl
  paSim <- mean(predFitted(fit, snr))

  expect_equal(paInt, paSim, tolerance = 0.05)
})

# --- 2. countDetections -------------------------------------------------------

test_that("countDetections truncates using the column it was told to use", {
  # Regression test. This hardcoded det$snr, so any other column name returned
  # NULL, subset(det, NULL >= theta) returned zero rows, and Nc came back 0
  # with no error at all.
  det <- data.frame(
    t0  = 738000 + seq_len(100) / 24,   # MATLAB datenums, spread over ~4 days
    SNR = seq(-10, 20, length.out = 100)
  )

  n <- countDetections(det, season = 'year',
                       snrTruncationThreshold = 0, snrColName = 'SNR')

  expect_equal(n, sum(det$SNR >= 0))
  expect_gt(n, 0)
})

test_that("countDetections errors when the named SNR column is absent", {
  det <- data.frame(t0 = 738000 + seq_len(10) / 24,
                    SNR = seq(-5, 5, length.out = 10))
  expect_error(
    countDetections(det, snrTruncationThreshold = 0, snrColName = 'snr'),
    "requires column"
  )
})

test_that("countDetections without truncation counts everything", {
  det <- data.frame(t0 = 738000 + seq_len(50) / 24,
                    snr = seq(-10, 20, length.out = 50))
  expect_equal(countDetections(det), 50)
})

# --- 3. The Nc guard ----------------------------------------------------------

test_that("cde refuses a threshold without confirmation that Nc is truncated", {
  skip_on_cran()
  d <- make_capture_history(n = 2e3)
  TL <- tlSpherical(rangeStep = 5000)

  expect_error(
    cde(Nc = d$Nc, capHistTab = d$capHistTab, SL = d$SL, TL = TL,
        T = d$Time, A = d$A, snrTruncationThreshold = 0),
    "NcIsTruncated"
  )
})

test_that("cde default path is untouched by the guard", {
  skip_on_cran()
  d <- make_capture_history(n = 2e3)
  TL <- tlSpherical(rangeStep = 5000)

  expect_no_error(
    suppressWarnings(suppressMessages(
      cde(Nc = d$Nc, capHistTab = d$capHistTab, SL = d$SL, TL = TL,
          T = d$Time, A = d$A, NL = d$NL, outerloop = 2)
    ))
  )
})

# --- 4. Dc recovery: the point of the whole exercise --------------------------

test_that("truncated cde returns all-call density, not above-threshold density", {
  # If the q(theta) cancellation is wrong, Dc comes back scaled by q(theta),
  # which at this threshold is a factor of several. The oracle detection
  # function and true noise are supplied so that the only thing under test is
  # the truncation bookkeeping.
  skip_on_cran()

  theta <- 3
  d  <- make_capture_history(n = 5e4, seed = 11,
                             det2location = 2, det2scale = 4)
  TL <- tlSpherical(rangeStep = 2000)

  # Observer 2's true curve, fitted as a glm because pDetInArea bootstraps
  # coefficients and will not accept a closure. Correctly specified on a large
  # sample, so this is the truth to within a rounding error. The curve is given
  # rather than estimated because this test is about truncation bookkeeping,
  # not about recovering the detection function.
  trueDet2 <- fitLogisticDetector(function(snr) plogis(snr, location = 2, scale = 4))

  ch <- d$capHistTab
  NcTrunc <- sum(ch$detect_table2 & ch$SNR >= theta, na.rm = TRUE)

  fit <- suppressWarnings(suppressMessages(
    cde(Nc = NcTrunc, capHistTab = ch, SL = d$SL, TL = TL,
        T = d$Time, A = d$A, NL = d$NL,
        snrDetFun = trueDet2,
        snrTruncationThreshold = theta,
        NcIsTruncated = TRUE,
        outerloop = 10)
  ))

  DcFraction <- fit$Dc / d$trueDc
  expect_gt(DcFraction, 0.85)
  expect_lt(DcFraction, 1.15)
})
