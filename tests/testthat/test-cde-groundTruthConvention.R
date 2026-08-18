# Tests for the detect_table1-as-ground-truth convention.
#
# cde()'s internal falseDiscoveryRate()/chtToSNRinfo() calls read
# capHistTab$detect_table1 as ground truth and detect_table2 as the detector
# under investigation BY DEFAULT. Since chtToSNRinfo(), cde() also accepts
# groundTruthCol/observerCol to point at different columns (e.g. a native
# matchbox verdict/detect_observerN table needs no detect_table1/2 renaming
# at all) -- but the default remains detect_table1/detect_table2, and the
# regression this file guards against is independent of that: it's about
# WHAT VALUES occupy the ground-truth column, not what it's named.
#
# For an observer-ground (OG) analysis, one trusted observer's own column
# genuinely is ground truth, so the default holds automatically. For an
# adjudicated capture-recapture (CR) analysis with two or more fallible
# observers, no raw observer column is ground truth on its own, and
# whichever column cde() is pointed at (via groundTruthCol, or by
# overwriting detect_table1 the older way) must contain the adjudicator's
# verdict -- exactly what mchToCR() did in the Common Ground manuscript's
# own analysis script, and what chtToSNRinfo()'s groundTruth argument does
# directly now.
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
