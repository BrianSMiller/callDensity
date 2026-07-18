# Tests for the detect_table1-as-ground-truth convention.
#
# cde()'s internal falseDiscoveryRate()/capHist2snrInfo() calls always read
# capHistTab$detect_table1 as ground truth and detect_table2 as the detector
# under investigation. Neither function is told otherwise by cde(), which
# does not expose gtColName/testColName at all.
#
# For an observer-ground (OG) analysis this holds automatically, since one
# trusted observer's own column genuinely is ground truth. For an adjudicated
# capture-recapture (CR) analysis with two fallible observers, neither raw
# observer is ground truth, and detect_table1 must be deliberately overwritten
# with the adjudicator's verdict before the table reaches cde() -- exactly
# what mchToCR() does in the Common Ground manuscript's own analysis script,
# and what the callDensity_snrThreshold vignette's CR section now does too.
#
# Discovered because it wasn't done in an earlier draft of both the vignette
# and notes/truncationSweep.R: cde() silently computed a false discovery rate
# of ~0.55 for a simulation with a true fdr of 0, because it was comparing one
# fallible observer against another rather than against genuine ground truth.
# No error, no warning -- just a wrong number that looked plausible on its
# own. These tests exist so a future regression looks like this again, loudly,
# in CI, rather than being rediscovered by hand.

test_that("CR capHistTab with a raw observer as detect_table1 gives the wrong c", {
  skip_on_cran()
  d  <- make_capture_history(n = 2e4, seed = 7, det1location = 1, det1scale = 1,
                             det2location = 2, det2scale = 4, fdr = 0)
  ch <- d$capHistTab

  # fdr = 0 in simulation: there are no genuine false positives anywhere, so
  # the TRUE false discovery rate for detector 2 is exactly 0.
  trueFdr <- sum(ch$detect_table2 & !ch$groundTruth2) / sum(ch$detect_table2)
  expect_equal(trueFdr, 0)

  # WRONG construction: detect_table1 is left as the raw first observer, not
  # the adjudicator's verdict. This is the mistake this test exists to catch.
  fdrWrong <- falseDiscoveryRate(ch, "year", -Inf)
  expect_gt(fdrWrong$c, 0.1)  # comparing two fallible observers is not ~0

  # RIGHT construction: detect_table1 overwritten with genuine ground truth.
  chRight <- ch
  chRight$detect_table1 <- chRight$groundTruth2
  fdrRight <- falseDiscoveryRate(chRight, "year", -Inf)
  expect_equal(fdrRight$c, 0, tolerance = 1e-8)
})

test_that("cde() gives the wrong Dc for CR data unless detect_table1 is the adjudicated verdict", {
  skip_on_cran()
  d  <- make_capture_history(n = 2e4, seed = 7, det1location = 1, det1scale = 1,
                             det2location = 2, det2scale = 4, fdr = 0)
  ch <- d$capHistTab
  TL <- tlSpherical(rangeStep = 5000)

  adj <- subset(ch, detect_table1 | detect_table2)
  fit <- suppressWarnings(fitSNRvglm(adj, c("detect_table1", "detect_table2"),
                                     whichObserver = "detect_table2"))
  Nc  <- sum(ch$detect_table2)

  resultWrong <- suppressWarnings(suppressMessages(
    cde(Nc = Nc, capHistTab = ch, snrDetFun = fit,
        SL = d$SL, TL = TL, NL = d$NL, T = d$Time, A = d$A, outerloop = 5)
  ))

  chRight <- ch
  chRight$detect_table1 <- chRight$groundTruth2
  resultRight <- suppressWarnings(suppressMessages(
    cde(Nc = Nc, capHistTab = chRight, snrDetFun = fit,
        SL = d$SL, TL = TL, NL = d$NL, T = d$Time, A = d$A, outerloop = 5)
  ))

  # The wrong construction inflates c (comparing two fallible observers, not
  # against ground truth), which pulls Dc away from the truth relative to the
  # correct construction. This doesn't assert a specific direction or
  # magnitude for c itself -- just that using the true ground truth column
  # lands unambiguously closer to the real density than not doing so.
  expect_true(
    abs(resultRight$Dc - d$trueDc) < abs(resultWrong$Dc - d$trueDc)
  )
})
