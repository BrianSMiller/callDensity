# Tests for the union call density additions: falseDiscoveryRate()'s
# testColName generalized to a vector, the new unionDetections() helper,
# and cde()'s check that observerCol (single vs union) and snrDetFun's own
# whichObserver tag (single vs "any") agree with each other.

# ---------------------------------------------------------------------------
# A small, hand-countable fixture: verdict (ground truth) plus two detector
# columns, chosen so the union FP/TP/Nc can be checked by counting rows by
# hand, not just internal consistency.
#   row  verdict  det1  det2   -- union(det1,det2)?  true positive?
#    1      1       1     0     yes                    yes (TP)
#    2      1       0     1     yes                    yes (TP)
#    3      1       1     1     yes                    yes (TP)
#    4      1       0     0     no  (neither flagged)   -- (not counted either way)
#    5      0       1     0     yes                    no  (FP)
#    6      0       0     1     yes                    no  (FP)
#    7      0       0     0     no                      -- (not counted either way)
# Union positive: rows 1,2,3,5,6 (5 rows). Of these, TP = rows 1,2,3 (3),
# FP = rows 5,6 (2). Union c = FP/(FP+TP) = 2/5 = 0.4.
# Individual det1 positive: rows 1,3,5 -- TP=2 (rows 1,3), FP=1 (row 5), c=1/3.
# Individual det2 positive: rows 2,3,6 -- TP=2 (rows 2,3), FP=1 (row 6), c=1/3.
# ---------------------------------------------------------------------------
make_union_toy_cht <- function() {
  data.frame(
    verdict = c(1, 1, 1, 1, 0, 0, 0),
    detect_observer1 = c(1, 0, 1, 0, 1, 0, 0),
    detect_observer2 = c(0, 1, 1, 0, 0, 1, 0),
    t = as.POSIXct("2025-06-01", tz = "UTC") + 0:6 * 3600
  )
}

test_that("falseDiscoveryRate with a single testColName is unchanged", {
  cht <- make_union_toy_cht()
  fdr1 <- falseDiscoveryRate(cht, gtColName = "verdict", testColName = "detect_observer1")
  expect_equal(fdr1$c, 1/3)
})

test_that("falseDiscoveryRate with a vector testColName computes the union", {
  cht <- make_union_toy_cht()
  fdrUnion <- falseDiscoveryRate(cht, gtColName = "verdict",
                                 testColName = c("detect_observer1", "detect_observer2"))
  expect_equal(fdrUnion$c, 2/5)
})

test_that("union false discovery rate has a larger sample size than either individual one", {
  # Not just a different number -- the union pools strictly more positive
  # calls than either contributor alone (5 vs 3 vs 3 in this fixture), which
  # is the whole basis for its tighter CV.
  cht <- make_union_toy_cht()
  fdr1     <- falseDiscoveryRate(cht, gtColName = "verdict", testColName = "detect_observer1")
  fdr2     <- falseDiscoveryRate(cht, gtColName = "verdict", testColName = "detect_observer2")
  fdrUnion <- falseDiscoveryRate(cht, gtColName = "verdict",
                                 testColName = c("detect_observer1", "detect_observer2"))
  expect_true(fdrUnion$cv.c < fdr1$cv.c)
  expect_true(fdrUnion$cv.c < fdr2$cv.c)
})

test_that("unionDetections requires explicit fullCoverageConfirmed", {
  cht <- make_union_toy_cht()
  expect_error(
    unionDetections(cht, c("detect_observer1", "detect_observer2")),
    "fullCoverageConfirmed"
  )
})

test_that("unionDetections counts events flagged by at least one detector", {
  cht <- make_union_toy_cht()
  # Rows where detect_observer1 | detect_observer2: 1,2,3,5,6 -- 5 rows,
  # regardless of verdict (Nc counts positive predictions, true or false,
  # matching the existing single-detector Nc convention).
  n <- unionDetections(cht, c("detect_observer1", "detect_observer2"),
                       fullCoverageConfirmed = TRUE)
  expect_equal(n, 5)
})

test_that("unionDetections is between max(individual) and sum(individual)", {
  cht <- make_union_toy_cht()
  n1 <- sum(cht$detect_observer1)
  n2 <- sum(cht$detect_observer2)
  nUnion <- unionDetections(cht, c("detect_observer1", "detect_observer2"),
                            fullCoverageConfirmed = TRUE)
  expect_true(nUnion >= max(n1, n2))
  expect_true(nUnion <= n1 + n2)
})

test_that("unionDetections handles NA (a detector with no row for an event) correctly", {
  # A detector's own NA (never examined) alongside the other's TRUE should
  # still count as a union positive -- NA | TRUE is TRUE regardless of the
  # NA, matching R's own | semantics.
  cht <- data.frame(
    detect_observer1 = c(NA, NA, FALSE),
    detect_observer2 = c(TRUE, FALSE, FALSE)
  )
  n <- unionDetections(cht, c("detect_observer1", "detect_observer2"),
                       fullCoverageConfirmed = TRUE)
  expect_equal(n, 1)  # only the first row is a genuine union positive
})

# ---------------------------------------------------------------------------
# cde()'s consistency check between observerCol and snrDetFun's whichObserver
# ---------------------------------------------------------------------------

test_that("cde errors when observerCol is a union but snrDetFun's whichObserver is not 'any'", {
  skip_on_cran()
  d <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  ch <- d$capHistTab
  ch$SNR <- rowMeans(ch[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  adj <- subset(ch, detect_table1 | detect_table2)
  adj$SNR <- rowMeans(adj[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  fit <- suppressWarnings(fitSNRvglm(adj, c("detect_table1", "detect_table2"),
                                     whichObserver = "detect_table2"))

  expect_error(
    cde(Nc = d$Nc, capHistTab = ch, snrDetFun = fit,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 2,
        observerCol = c("detect_table1", "detect_table2")),
    "whichObserver"
  )
})

test_that("cde errors when snrDetFun's whichObserver is 'any' but observerCol is a single column", {
  skip_on_cran()
  d <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  ch <- d$capHistTab
  ch$SNR <- rowMeans(ch[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  adj <- subset(ch, detect_table1 | detect_table2)
  adj$SNR <- rowMeans(adj[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  fitUnion <- suppressWarnings(fitSNRvglm(adj, c("detect_table1", "detect_table2"),
                                          whichObserver = "any"))

  expect_error(
    cde(Nc = d$Nc, capHistTab = ch, snrDetFun = fitUnion,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 2,
        observerCol = "detect_table2"),
    "whichObserver"
  )
})

test_that("cde does not error on the check when observerCol and whichObserver agree (both union)", {
  skip_on_cran()
  d <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  ch <- d$capHistTab
  ch$SNR <- rowMeans(ch[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  adj <- subset(ch, detect_table1 | detect_table2)
  adj$SNR <- rowMeans(adj[, c("snr_table1", "snr_table2")], na.rm = TRUE)
  fitUnion <- suppressWarnings(fitSNRvglm(adj, c("detect_table1", "detect_table2"),
                                          whichObserver = "any"))
  NcUnion <- unionDetections(ch, c("detect_table1", "detect_table2"),
                             fullCoverageConfirmed = TRUE)

  expect_no_error(suppressWarnings(suppressMessages(
    cde(Nc = NcUnion, capHistTab = ch, snrDetFun = fitUnion,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 2,
        observerCol = c("detect_table1", "detect_table2"),
        signalCol = c("signalRMSdB_table1", "signalRMSdB_table2"),
        noiseCol  = c("noiseRMSdB_table1",  "noiseRMSdB_table2"))
  )))
})

test_that("cde's whichObserver check does not fire for glm/gam/scam models", {
  skip_on_cran()
  d <- make_capture_history(n = 2e4, seed = 5)
  TL <- tlSpherical(rangeStep = 5000)
  snrData <- make_snr_data()
  fit <- fitDetFun(snrData, modelType = "scam", numKnots = 5)
  ch <- d$capHistTab
  ch$SNR <- rowMeans(ch[, c("snr_table1", "snr_table2")], na.rm = TRUE)

  # A vector observerCol with a non-vglm model doesn't make conceptual
  # sense either, but that's a different concern -- this check specifically
  # only applies to vglm's whichObserver tag (isS4()), so it must not fire
  # here regardless.
  expect_no_error(suppressWarnings(suppressMessages(
    cde(Nc = d$Nc, capHistTab = ch, snrDetFun = fit,
        SL = d$SL, TL = TL, T = d$Time, A = d$A, outerloop = 2,
        observerCol = c("detect_table1", "detect_table2"),
        signalCol = c("signalRMSdB_table1", "signalRMSdB_table2"),
        noiseCol  = c("noiseRMSdB_table1",  "noiseRMSdB_table2"))
  )))
})
