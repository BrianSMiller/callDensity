# Tests for cdeFromParamFile.
#
# cdeFromParamFile was a near-complete duplicate of cde()'s body: it called
# pa_CV() without transect area weights (cde always passes wt =
# truncationDistance^2), used nlFromSnrInfo directly rather than picking up
# the nlFromDetections default, and had no equivalent of cde's NcIsTruncated
# guard. Being deprecated, it is now a thin wrapper: derive cde's arguments
# from p, and call cde. These tests check the derivation and delegation, not
# the underlying file-reading helpers (countDetections,
# deploymentDurationFromsoundFolderCsv, readCapHist), which are mocked here
# and are otherwise these functions' own concern, not this wrapper's.

test_that("cdeFromParamFile is deprecated and delegates to cde with derived arguments", {
  skip_on_cran()

  fakeNc <- 137
  fakeT  <- 42
  fakeCapHistTab <- data.frame(SNR = c(1, 2, 3), detect_table2 = c(TRUE, TRUE, FALSE))

  testthat::local_mocked_bindings(
    countDetections = function(detectionFile, season, snrTruncationThreshold) {
      expect_equal(detectionFile, "fake_detections.csv")
      expect_equal(snrTruncationThreshold, 3)
      fakeNc
    },
    deploymentDurationFromsoundFolderCsv = function(fullYearEffortFile, season) {
      expect_equal(fullYearEffortFile, "fake_effort.csv")
      fakeT
    },
    readCapHist = function(capHistFile) {
      expect_equal(capHistFile, "fake_capHist.csv")
      fakeCapHistTab
    }
  )

  capturedArgs <- NULL
  testthat::local_mocked_bindings(
    cde = function(...) {
      capturedArgs <<- list(...)
      data.frame(Dc = 0.5, CV.Dc = 0.1)  # dummy result, delegation is what's tested
    }
  )

  tlFile <- tempfile(fileext = ".csv")
  write.csv(data.frame(range_m = c(0, 1000), transect1 = c(100, 110)),
           tlFile, row.names = FALSE)

  p <- list(
    detectorParams = list(fullYearDetectionCsv = "fake_detections.csv",
                          fullYearEffortFile = "fake_effort.csv"),
    capHistFile = "fake_capHist.csv",
    slParams = list(slMean = 190, slStd = 4, slSampleSize = 350),
    tlParams = list(tlFile = tlFile),
    w = 1e6, k = 1, modelType = "scam", numKnots = 5,
    output.resolution.m = 100, outerloop = 10,
    transectFile = NULL, simResultsFile = NULL, paFile = NULL,
    siteCode = "testSite", densityResultsFile = NULL
  )

  expect_warning(
    result <- cdeFromParamFile(p, season = "year", truncationDistance = 1000,
                               snrTruncationThreshold = 3),
    "deprecated"
  )

  expect_equal(capturedArgs$Nc, fakeNc)
  expect_equal(capturedArgs$T, fakeT)
  expect_identical(capturedArgs$capHistTab, fakeCapHistTab)
  expect_equal(capturedArgs$SL, data.frame(mean = 190, sd = 4, sampleSize = 350))
  expect_equal(capturedArgs$k, 1)
  expect_equal(capturedArgs$modelType, "scam")
  expect_equal(capturedArgs$siteCode, "testSite")
  expect_equal(capturedArgs$truncationDistance, 1000)
  expect_equal(capturedArgs$snrTruncationThreshold, 3)

  # Nc was derived via countDetections using this same snrTruncationThreshold,
  # so it is truncated consistently by construction -- the wrapper must
  # confirm this to cde() itself, not leave the caller to work it out.
  expect_true(capturedArgs$NcIsTruncated)
})

test_that("cdeFromParamFile passes NcIsTruncated = FALSE when no threshold is set", {
  skip_on_cran()

  testthat::local_mocked_bindings(
    countDetections = function(...) 100,
    deploymentDurationFromsoundFolderCsv = function(...) 10,
    readCapHist = function(...) data.frame(SNR = 1)
  )

  capturedArgs <- NULL
  testthat::local_mocked_bindings(
    cde = function(...) { capturedArgs <<- list(...); data.frame(Dc = 1) }
  )

  tlFile <- tempfile(fileext = ".csv")
  write.csv(data.frame(range_m = c(0, 1000), transect1 = c(100, 110)),
           tlFile, row.names = FALSE)

  p <- list(
    detectorParams = list(fullYearDetectionCsv = "x", fullYearEffortFile = "y"),
    capHistFile = "z", slParams = list(slMean = 190, slStd = 4, slSampleSize = 350),
    tlParams = list(tlFile = tlFile), w = 1e6, k = 1, modelType = "scam",
    numKnots = 5, output.resolution.m = 100, outerloop = 10,
    transectFile = NULL, simResultsFile = NULL, paFile = NULL,
    siteCode = "x", densityResultsFile = NULL
  )

  suppressWarnings(
    cdeFromParamFile(p, season = "year", truncationDistance = 1000)
    # snrTruncationThreshold left at its default, -Inf
  )

  expect_false(capturedArgs$NcIsTruncated)
})
