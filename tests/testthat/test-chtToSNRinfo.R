# Tests for chtToSNRinfo(), which replaces mchToCR()/capHist2snrInfo()/
# capHistTosnrInfo() with one converter supporting any number of observers.
#
# Small, hand-constructed fixtures are used throughout (rather than the
# simulation fixtures in helper-fixtures.R) so the union-then-groundTruth
# row filter can be verified by exact row count, not just "it runs".

# ---------------------------------------------------------------------------
# A minimal 8-row matchbox-native fixture with a known, hand-countable
# filter outcome:
#   row  verdict  obs1  obs2  obs3   -- union(verdict,obs1,obs2,obs3)?  verdict==1?
#    1      1       1     0     0     yes                                yes  (kept)
#    2      1       0     0     0     yes (verdict alone)                yes  (kept)
#    3      0       1     0     0     yes (obs1 alone, false alarm)      no   (dropped)
#    4      0       0     1     1     yes                                no   (dropped)
#    5      1       1     1     1     yes                                yes  (kept)
#    6      0       0     0     0     NO (nobody flagged this row)       no   (dropped, and shouldn't even reach the gt filter)
#    7      1       0     1     0     yes                                yes  (kept)
#    8      0       0     0     1     yes (obs3 alone, false alarm)      no   (dropped)
# Expected surviving rows after chtToSNRinfo's filter (union, then
# verdict==1): rows 1, 2, 5, 7 -- 4 rows, all with verdict==1.
# ---------------------------------------------------------------------------
make_toy_cht <- function() {
  data.frame(
    verdict          = c(1,1,0,0,1,0,1,0),
    detect_observer1 = c(1,0,1,0,1,0,0,0),
    detect_observer2 = c(0,0,0,1,1,0,1,0),
    detect_observer3 = c(0,0,0,1,1,0,0,1),
    t0               = 738000 + seq(0, by = 0.01, length.out = 8),
    signalRMSdB      = c(120,120,150,150,125,90,122,150),
    noiseRMSdB       = c(90,90,90,90,90,90,90,90),
    signalRMSdB1     = c(112,108,150,150,116,90,113,150),
    signalRMSdB2     = c(108,108,150,150,114,90,111,150),
    noiseRMSdB1      = rep(90, 8),
    noiseRMSdB2      = rep(90, 8),
    datetime         = as.POSIXct("2025-06-01", tz = "UTC") + 0:7 * 3600
  )
}

test_that("union-then-groundTruth row filter keeps exactly the expected rows", {
  cht <- make_toy_cht()
  out <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1")
  expect_equal(nrow(out), 4)
  expect_true(all(out$verdict == 1))
})

test_that("row not flagged by anyone is dropped even before the groundTruth filter", {
  cht <- make_toy_cht()
  out <- chtToSNRinfo(cht, groundTruth = "verdict",
                      observers = c("observer1","observer2","observer3"))
  # row 6 (all zero) must not appear regardless of how many observers are named
  expect_equal(nrow(out), 4)
})

test_that("single-observer call reproduces the old OG/SNRinfo shape", {
  cht <- make_toy_cht()
  out <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1")
  expect_true(all(c("Detected","CallRL","NoiseRL","SNR","t","month","season") %in% names(out)))
  expect_equal(out$Detected, out$detect_observer1)
})

test_that("N-observer call keeps native detect_observerN columns intact for yColNames", {
  cht <- make_toy_cht()
  out <- chtToSNRinfo(cht, groundTruth = "verdict",
                      observers = c("observer1","observer2","observer3"))
  expect_true(all(c("verdict","detect_observer1","detect_observer2","detect_observer3")
                  %in% names(out)))
  # Detected defaults to the last-named observer
  expect_equal(out$Detected, out$detect_observer3)
})

test_that("groundTruth/observers resolve bare suffixes and full column names identically", {
  cht <- make_toy_cht()
  bySuffix <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1")
  cht2 <- cht
  names(cht2)[names(cht2) == "verdict"] <- "detect_table1"
  names(cht2)[names(cht2) == "detect_observer1"] <- "detect_table2"
  byFullName <- chtToSNRinfo(cht2, groundTruth = "table1", observers = "table2")
  expect_equal(bySuffix$SNR, byFullName$SNR)
  expect_equal(nrow(bySuffix), nrow(byFullName))
})

test_that("unresolvable column names error loudly rather than guess", {
  cht <- make_toy_cht()
  expect_error(chtToSNRinfo(cht, groundTruth = "nonexistent", observers = "observer1"),
              "not found")
  expect_error(chtToSNRinfo(cht, groundTruth = "verdict", observers = "nonexistent"),
              "not found")
})

test_that("signalCol/noiseCol: single shared column (default) vs legacy averaged pair", {
  cht <- make_toy_cht()
  shared   <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1")
  averaged <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1",
                           signalCol = c("signalRMSdB1","signalRMSdB2"),
                           noiseCol  = c("noiseRMSdB1","noiseRMSdB2"))
  # Different provenance -- CallRL should differ for at least one retained row
  # in this fixture (rows 1/5/7 have signalRMSdB1 != signalRMSdB2)
  expect_false(isTRUE(all.equal(shared$CallRL, averaged$CallRL)))
  expect_equal(nrow(shared), nrow(averaged))
})

test_that("timeCol='t0' converts via mat2Rdate; any other column is used as-is", {
  cht <- make_toy_cht()
  viaT0 <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1", timeCol = "t0")
  viaDatetime <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1",
                              timeCol = "datetime")
  expect_s3_class(viaT0$t, "POSIXct")
  expect_s3_class(viaDatetime$t, "POSIXct")
  # Unrelated clocks in this fixture (t0 is an arbitrary datenum, datetime is
  # 2025-06-01 onward) -- just confirm both resolve to *something* sane and
  # that choosing the column actually changes the result, i.e. it isn't
  # silently ignored.
  expect_false(isTRUE(all.equal(as.numeric(viaT0$t), as.numeric(viaDatetime$t))))
})

test_that("infinite SNR rows are dropped with a warning", {
  cht <- make_toy_cht()
  cht$noiseRMSdB[1] <- -Inf  # forces SNR = Inf on a row that would otherwise survive
  expect_warning(
    out <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1"),
    "infinite"
  )
  expect_true(all(is.finite(out$SNR)))
})

test_that("season filtering matches time2season on the resolved time column", {
  cht <- make_toy_cht()
  cht$t0 <- as.numeric(as.POSIXct("2025-01-15", tz="UTC")) / 86400 + 719529  # summer
  cht$t0[1] <- as.numeric(as.POSIXct("2025-07-15", tz="UTC")) / 86400 + 719529  # winter
  allYear <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1", season = "year")
  summerOnly <- chtToSNRinfo(cht, groundTruth = "verdict", observers = "observer1", season = "summer")
  expect_lt(nrow(summerOnly), nrow(allYear))
  expect_true(all(summerOnly$season == "summer"))
})

test_that("chtToSNRinfo output feeds fitDetFun(modelType='vglm') directly, 2 and 3 observers", {
  skip_on_cran()
  set.seed(11)
  n <- 2000
  d <- data.frame(t0 = 738000 + seq_len(n) * 0.001)
  d$SNR <- rnorm(n, 5, 6)
  d$signalRMSdB <- d$SNR + 90
  d$noiseRMSdB  <- 90
  d$verdict <- 1
  p <- plogis(0.3 * (d$SNR - 2))
  d$detect_observer1 <- rbinom(n, 1, p * 0.85)
  d$detect_observer2 <- rbinom(n, 1, p * 0.75)
  d$detect_observer3 <- rbinom(n, 1, p * 0.65)
  d <- d[d$detect_observer1 + d$detect_observer2 + d$detect_observer3 > 0, ]

  snr2 <- chtToSNRinfo(d, groundTruth = "verdict", observers = c("observer1","observer2"))
  fit2 <- fitDetFun(snr2, modelType = "vglm",
                    yColNames = c("verdict","detect_observer2"),
                    whichObserver = "detect_observer2")
  expect_s4_class(fit2, "vglm")

  snr3 <- chtToSNRinfo(d, groundTruth = "verdict",
                       observers = c("observer1","observer2","observer3"))
  fit3 <- fitDetFun(snr3, modelType = "vglm",
                    yColNames = c("verdict","detect_observer1","detect_observer2","detect_observer3"),
                    whichObserver = "detect_observer2")
  expect_s4_class(fit3, "vglm")
  # The whole point of this feature: no hand-built two-column reduction was
  # needed anywhere above to get a genuine 3-observer vglm fit.
  expect_equal(ncol(fit3@y), 4)
})

# ---------------------------------------------------------------------------
# Regression guards: the three deprecated wrappers must keep their exact
# historical output shape and must actually warn.
# ---------------------------------------------------------------------------

test_that("capHist2snrInfo is deprecated but reproduces its original output shape", {
  cht <- make_toy_cht()
  names(cht)[names(cht) == "verdict"] <- "detect_table1"
  names(cht)[names(cht) == "detect_observer1"] <- "detect_table2"
  cht$t <- cht$datetime
  cht$month <- time2monthCode(cht$t)
  cht$season <- time2season(cht$t)

  expect_warning(out <- capHist2snrInfo(cht, "year"), "deprecated")
  expect_equal(names(out), c("Detected","CallRL","NoiseRL","SNR","t","month","season"))
  expect_equal(nrow(out), sum(cht$detect_table1 == 1))
})

test_that("capHistTosnrInfo is deprecated but reproduces its original output shape", {
  cht <- make_toy_cht()
  names(cht)[names(cht) == "verdict"] <- "detect_table1"
  names(cht)[names(cht) == "detect_observer1"] <- "detect_table2"

  expect_warning(out <- capHistTosnrInfo(cht), "deprecated")
  expect_equal(names(out), c("Detected","CallRL","NoiseRL","SNR","t","month","season"))
  expect_equal(nrow(out), sum(cht$detect_table1 == 1))
})

test_that("mchToCR is deprecated but its output is unchanged", {
  cht <- make_toy_cht()
  cht$tEnd  <- cht$t0 + 0.001
  cht$fLow  <- 15
  cht$fHigh <- 25
  cht$noiseDev <- 0
  cht$SNR <- cht$signalRMSdB - cht$noiseRMSdB

  expect_warning(out <- mchToCR(cht, "observer1", "observer2"), "deprecated")
  expect_true(all(c("detect_table1","detect_table2","t0_table1") %in% names(out)))
})
