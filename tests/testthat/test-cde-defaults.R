# Tests for cde()'s default NL estimator.
#
# cde() had no direct tests before this file (see pa_CV fix commit and the
# TODO it logged: cde() and its CV.*/pa_CV/c_CV/Dc_CV helpers are otherwise
# untested). This is a narrow start, scoped to the nlFromSnrInfo ->
# nlFromDetections default swap specifically, not an attempt at full coverage.

test_that("cde() runs end to end with the default NL estimator, scam", {
  skip_on_cran()
  d  <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  snrData <- make_snr_data()
  fit <- fitDetFun(snrData, modelType = "scam", numKnots = 5)
  # simsTocaptureHistoryTable() no longer auto-computes a consolidated SNR --
  # needed by falseDiscoveryRate()'s snrTruncationThreshold handling
  # (a separate mechanism from chtToSNRinfo()'s signalCol/noiseCol averaging,
  # which cde() is pointed at explicitly below).
  d$capHistTab$SNR <- rowMeans(d$capHistTab[, c("snr_table1", "snr_table2")], na.rm = TRUE)

  result <- suppressWarnings(suppressMessages(
    cde(Nc = d$Nc, capHistTab = d$capHistTab, snrDetFun = fit,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 5,
        signalCol = c("signalRMSdB_table1", "signalRMSdB_table2"),
        noiseCol  = c("noiseRMSdB_table1",  "noiseRMSdB_table2"))
  ))

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(c(result$pa, result$Dc, result$NLmean, result$NLsd))))
  expect_gt(result$NLsd, 0)
})

test_that("cde() runs end to end with the default NL estimator, vglm", {
  skip_on_cran()
  d  <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  d$capHistTab$SNR <- rowMeans(d$capHistTab[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  adj <- subset(d$capHistTab, detect_table1 | detect_table2)
  fit <- suppressWarnings(fitSNRvglm(adj, c("detect_table1", "detect_table2"),
                                     whichObserver = "detect_table2"))

  result <- suppressWarnings(suppressMessages(
    cde(Nc = d$Nc, capHistTab = d$capHistTab, snrDetFun = fit,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 5,
        signalCol = c("signalRMSdB_table1", "signalRMSdB_table2"),
        noiseCol  = c("noiseRMSdB_table1",  "noiseRMSdB_table2"))
  ))

  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(c(result$pa, result$Dc, result$NLmean, result$NLsd))))
  expect_gt(result$NLsd, 0)
})

test_that("cde()'s default NL is nlFromDetections, not nlFromSnrInfo", {
  # Regression guard for the swap. The two estimators are motivated by
  # different mechanisms (see the comment in cde()) and generally give
  # different point estimates for the same input. Confirms cde() actually
  # calls the new one rather than silently keeping the old default.
  skip_on_cran()
  d  <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  snrData <- make_snr_data()
  fit <- fitDetFun(snrData, modelType = "scam", numKnots = 5)
  d$capHistTab$SNR <- rowMeans(d$capHistTab[, c("snr_table1", "snr_table2")], na.rm = TRUE)

  SNRinfo <- chtToSNRinfo(d$capHistTab, groundTruth = "table1", observers = "table2",
                          signalCol = c("signalRMSdB_table1", "signalRMSdB_table2"),
                          noiseCol  = c("noiseRMSdB_table1",  "noiseRMSdB_table2"))

  nlOld <- nlFromSnrInfo(SNRinfo, fit)
  nlNew <- nlFromDetections(SNRinfo, fit, d$SL, TL)

  result <- suppressWarnings(suppressMessages(
    cde(Nc = d$Nc, capHistTab = d$capHistTab, snrDetFun = fit,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 5,
        signalCol = c("signalRMSdB_table1", "signalRMSdB_table2"),
        noiseCol  = c("noiseRMSdB_table1",  "noiseRMSdB_table2"))
  ))

  expect_equal(result$NLmean, nlNew$mean, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(result$NLmean, nlOld$mean, tolerance = 1e-6)))
})
