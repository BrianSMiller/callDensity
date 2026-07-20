# Shared fixture: spherical spreading, blue whale source levels, a plain
# logistic detector. Deliberately not using helper-fixtures.R, because the
# point of the first test is to be an independent implementation.
tlSpherical <- function(R = 1e6, rangeStep = 1000, numTransects = 4) {
  TL <- simTLradials_20logR(maxRange = R, rangeStep = rangeStep,
                            numTransects = numTransects)
  TL[is.finite(rowSums(TL)) & TL[[1]] > 0, ]
}
testSL       <- data.frame(mean = 190, sd = 4, sampleSize = 350)
testDetector <- function(snr) plogis(snr, location = 1, scale = 2)

test_that("pDetGivenNL agrees with a brute force simulation", {
  # This is the test that can catch a wrong idea rather than a typo.
  # pDetGivenNL reaches the answer by area weighting, quadrature and binning.
  # This reaches it by putting calls on the water and rolling dice. They share
  # no code. If they agree, the bookkeeping is right.
  skip_on_cran()
  set.seed(42)

  R  <- 1e6
  TL <- tlSpherical(R)
  nlFixed <- 84

  n   <- 5e5
  r   <- R * sqrt(runif(n))                      # uniform in 2D
  sl  <- rnorm(n, testSL$mean, testSL$sd)
  snr <- sl - 20 * log10(r) - nlFixed
  pSim <- mean(runif(n) < testDetector(snr))

  pInt <- pDetGivenNL(nlFixed, testDetector, SL=testSL, TL=TL)

  expect_equal(pInt, pSim, tolerance = 0.02)
})

test_that("pDetGivenNL falls one decade per 10 dB under spherical spreading", {
  # The analytic result. Under 20*log10(r) with calls uniform in 2D, 10 dB more
  # noise means 3.16 times less range, so 10 times less area, so 10 times fewer
  # detections. Evaluated away from the saturation shoulder at low noise.
  # This test fails if the area weight is ever changed from r.
  TL <- tlSpherical()
  p  <- pDetGivenNL(c(90, 100, 110), testDetector, SL=testSL, TL=TL)
  expect_equal(p[-3] / p[-1], c(10, 10), tolerance = 0.05)
})

test_that("pDetGivenNL is insensitive to bin width", {
  # Spherical spreading with a 1 km range step is the awkward case: at short
  # range consecutive steps are several dB apart in TL, so the low-TL bins hold
  # one cell or none and rounding to the bin centre actually costs something.
  # Real TL tables are far more densely packed. Tolerance is deliberately loose
  # enough to pass on this case and tight enough to catch a binning bug, which
  # would move things by percent rather than by parts per thousand.
  TL <- tlSpherical()
  expect_equal(pDetGivenNL(84, testDetector, SL=testSL, TL=TL, binWidth = 1),
               pDetGivenNL(84, testDetector, SL=testSL, TL=TL, binWidth = 0.1),
               tolerance = 5e-3)
})

test_that("pDetGivenNL falls monotonically with noise", {
  TL <- tlSpherical()
  p  <- pDetGivenNL(seq(60, 120, by = 5), testDetector, SL=testSL, TL=TL)
  expect_true(all(diff(p) < 0))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("pDetGivenNL is relative to the truncated area", {
  # Truncating removes distant cells, which carry a lot of area and few
  # detections. So the fraction of the remaining area that is audible goes up.
  TL <- tlSpherical()
  expect_gt(pDetGivenNL(84, testDetector, SL=testSL, TL=TL,
                        truncationDistance = 1e5),
            pDetGivenNL(84, testDetector, SL=testSL, TL=TL))
})

test_that("pDetGivenNL accepts a detFun object as well as a plain function", {
  TL      <- tlSpherical()
  snrInfo <- make_snr_data()
  p <- pDetGivenNL(84, fitDetFun(snrInfo, modelType = "glm"), SL=testSL, TL=TL)
  expect_true(p >= 0 && p <= 1)
})
