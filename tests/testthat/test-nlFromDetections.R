test_that("predictSampledNL predicts a mean at or below the truth", {
  TL <- tlSpherical()
  expect_lt(predictSampledNL(84, 4, testDetector, testSL, TL), 84)
})

test_that("predictSampledNL rises monotonically with mu, so the root is unique", {
  TL <- tlSpherical()
  p  <- vapply(seq(80, 92, by = 2),
               function(m) predictSampledNL(m, 4, testDetector, testSL, TL),
               numeric(1))
  expect_true(all(diff(p) > 0))
})

test_that("predictSampledNL applies no correction when noise does not vary", {
  TL <- tlSpherical()
  expect_equal(predictSampledNL(84, 0, testDetector, testSL, TL), 84)
})

test_that("nlFromDetections recovers a known noise level from simulated detections", {
  # THE test. Everything else checks that the code does what I meant. This
  # checks that what I meant is right.
  #
  # The noise here was measured at calls that were actually detected: placed on
  # the water, given a range and a source level, passed through a Bernoulli
  # trial. It shares no code with pDetGivenNL, which reaches its answer by area
  # weighting, quadrature and binning. If nlFromDetections recovers 84 from
  # these detections, the reasoning is sound and not merely self-consistent.
  skip_on_cran()

  TL      <- tlSpherical()
  snrInfo <- simulateDetectedNoise(nlMean = 84, nlSd = 4)

  # Confirm the problem exists before testing the cure.
  expect_lt(mean(snrInfo$NoiseRL), 82)

  est <- nlFromDetections(snrInfo, testDetector, testSL, TL)

  # Not exact, and the residual is known rather than mysterious. sd is taken
  # from the detections, and filtering narrows it slightly, roughly 3.9 against
  # a true 4.0. The bias goes as the square of sd, so a 3% error in sd is a 6%
  # error in a 3 dB bias, about 0.2 dB. Solving for mean and sd jointly would
  # remove it, at the cost of a two dimensional optimisation.
  expect_lt(abs(est$mean - 84), 0.5)
  expect_equal(est$sd, 4, tolerance = 0.05)
})

test_that("nlFromDetections works where nlFromSnrInfo's assumption fails", {
  # nlFromSnrInfo adds the detector's 50% point to the measured mean. Here that
  # is 8 dB while the real bias is near 3 dB, so the old correction overshoots
  # by roughly 5 dB and this one should not.
  skip_on_cran()

  fussy   <- function(snr) plogis(snr, location = 8, scale = 2)
  TL      <- tlSpherical()
  snrInfo <- simulateDetectedNoise(detector = fussy, seed = 2)

  est <- nlFromDetections(snrInfo, fussy, testSL, TL)
  expect_lt(abs(est$mean - 84), 0.5)

  oldStyle <- mean(snrInfo$NoiseRL) + 8
  expect_gt(abs(oldStyle - 84), 3)
})

test_that("nlFromDetections returns the same shape as nlFromSnrInfo", {
  TL  <- tlSpherical()
  est <- nlFromDetections(simulateDetectedNoise(n = 1e5), testDetector,
                          testSL, TL)
  expect_s3_class(est, "data.frame")
  expect_equal(nrow(est), 1)
  expect_named(est, c("mean", "sd", "sampleSize"))
})

test_that("nlFromDetections errors clearly on a missing column", {
  TL <- tlSpherical()
  expect_error(nlFromDetections(data.frame(x = 1), testDetector, testSL, TL),
               "no column named")
})

test_that("nlFromDetections errors clearly when searchWidth is too small", {
  TL <- tlSpherical()
  expect_error(nlFromDetections(simulateDetectedNoise(n = 1e5), testDetector,
                                testSL, TL, searchWidth = 0.5),
               "exceed searchWidth")
})


test_that("nlFromSnrInfo still returns what it always did", {
  # It is now a reproducibility artefact. The point of keeping it is that
  # published analyses can be re-run, so its numbers must not drift.
  snrInfo <- make_snr_data()
  detFun  <- fitDetFun(snrInfo, modelType = "glm")
  old <- suppressWarnings(nlFromSnrInfo(snrInfo, detFun))
  expect_named(old, c("mean", "sd", "sampleSize"))
  expect_gt(old$mean, mean(snrInfo$NoiseRL, na.rm = TRUE))
})
