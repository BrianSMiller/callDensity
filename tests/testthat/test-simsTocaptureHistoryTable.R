# Tests for simsTocaptureHistoryTable(), which merges N simulated detectors
# into a matchbox-native capture history table. No dedicated test file
# existed for this function before -- it had only ever been exercised
# indirectly via vignettes.

make_n_detectors <- function(n_detectors = 2, n_calls = 200, seed = 1) {
  R <- 1e5
  minDate <- as.POSIXct("2025-01-01")
  maxDate <- as.POSIXct("2025-01-02")
  SL <- data.frame(mean = 190, sd = 4, sampleSize = 350)
  NL <- data.frame(mean = 84, sd = 4, sampleSize = n_calls)

  set.seed(seed)
  sim <- simCallLocation(n = n_calls, R = R, minDate = minDate, maxDate = maxDate)
  sim <- simCallAcoustics(sim, SL, NL, TL = function(r) 20 * log10(r))

  # Spread quality/location so no two detectors are identical
  locs   <- seq(0, 3, length.out = n_detectors)
  scales <- seq(2, 4, length.out = n_detectors)

  lapply(seq_len(n_detectors), function(i) {
    set.seed(seed + i)
    params <- data.frame(location = locs[i], scale = scales[i],
                         func = "plogis", c = 0.1, fpMean = 0, fpSD = 4)
    simulateDetector(params, sim)
  })
}

test_that("default observerSuffix produces matchbox-native column names, N=2", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])
  expect_true(all(c("key", "t0", "tEnd", "detect_observer1", "detect_observer2",
                    "snr_observer1", "snr_observer2",
                    "groundTruth_observer1", "groundTruth_observer2") %in% names(cht)))
  # fLow/fHigh deliberately omitted -- this simulation never models frequency
  expect_false(any(c("fLow", "fHigh") %in% names(cht)))
})

test_that("generalizes to N=3 and N=4 detectors", {
  for (N in c(3, 4)) {
    d <- make_n_detectors(N)
    cht <- do.call(simsTocaptureHistoryTable, d)
    expect_true(all(paste0("detect_observer", seq_len(N)) %in% names(cht)))
    expect_true(all(paste0("snr_observer", seq_len(N)) %in% names(cht)))
  }
})

test_that("explicit observerSuffix overrides the default naming", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]], observerSuffix = c("table1", "table2"))
  expect_true(all(c("detect_table1", "detect_table2", "snr_table1", "snr_table2") %in% names(cht)))
  expect_false(any(c("detect_observer1", "detect_observer2") %in% names(cht)))
})

test_that("shared true calls merge into one row; each detector's own false positives stay separate", {
  d <- make_n_detectors(2, n_calls = 200)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])

  bothPresent <- !is.na(cht$groundTruth_observer1) & !is.na(cht$groundTruth_observer2)
  expect_equal(sum(bothPresent), 200)

  only1 <- !is.na(cht$groundTruth_observer1) & is.na(cht$groundTruth_observer2)
  only2 <- is.na(cht$groundTruth_observer1) & !is.na(cht$groundTruth_observer2)
  expect_equal(nrow(cht), 200 + sum(only1) + sum(only2))

  # A detector's own false positives are, by construction, not true calls
  expect_true(all(cht$groundTruth_observer1[only1] == FALSE))
  expect_true(all(cht$groundTruth_observer2[only2] == FALSE))
})

test_that("detect_<suffix> is always clean TRUE/FALSE, never NA", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])
  expect_false(anyNA(cht$detect_observer1))
  expect_false(anyNA(cht$detect_observer2))
  # A row that exists only because observer2 flagged a false positive must
  # still show observer1 as a clean FALSE, not NA -- observer1 examined
  # (and by construction, in this simulation, always covers) that instant
  # and simply didn't detect anything there.
  only2 <- is.na(cht$groundTruth_observer1) & !is.na(cht$groundTruth_observer2)
  expect_true(all(cht$detect_observer1[only2] == FALSE))
})

test_that("other (non-detect) columns are NA where a detector has no row for the event", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])
  only1 <- !is.na(cht$groundTruth_observer1) & is.na(cht$groundTruth_observer2)
  expect_true(anyNA(cht$snr_observer2[only1]))
  expect_true(anyNA(cht$signalRMSdB_observer2[only1]))
})

test_that("key is sequential and rows are sorted by datetime", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])
  expect_equal(cht$key, seq_len(nrow(cht)))
  expect_false(is.unsorted(cht$datetime))
})

test_that("t0 and tEnd are genuine MATLAB datenums matching datetime, matchbox-style", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])
  expect_equal(cht$t0, Rdate2mat(cht$datetime))
  expect_equal(cht$tEnd, Rdate2mat(cht$datetime))
  # Genuinely usable via chtToSNRinfo()'s default timeCol='t0' and cde()'s
  # timeCol=NULL auto-detection, not just a same-named POSIXct
  expect_true(is.numeric(cht$t0))
  expect_equal(as.numeric(mat2Rdate(cht$t0)), as.numeric(cht$datetime), tolerance = 1)
})

test_that("groundTruth is coalesced into a single event-level column, alongside per-observer copies", {
  d <- make_n_detectors(2)
  cht <- simsTocaptureHistoryTable(d[[1]], d[[2]])

  expect_true("groundTruth" %in% names(cht))
  expect_true(all(c("groundTruth_observer1", "groundTruth_observer2") %in% names(cht)))
  expect_false(anyNA(cht$groundTruth))

  # Whenever both observers have a row, they must agree -- this is what
  # makes coalescing safe rather than just convenient.
  both <- !is.na(cht$groundTruth_observer1) & !is.na(cht$groundTruth_observer2)
  expect_true(all(cht$groundTruth_observer1[both] == cht$groundTruth_observer2[both]))

  # Rows unique to one observer (that observer's own false positive) must
  # still resolve correctly via the other's NA, not propagate NA.
  only1 <- !is.na(cht$groundTruth_observer1) & is.na(cht$groundTruth_observer2)
  only2 <- is.na(cht$groundTruth_observer1) & !is.na(cht$groundTruth_observer2)
  expect_equal(cht$groundTruth[only1], cht$groundTruth_observer1[only1])
  expect_equal(cht$groundTruth[only2], cht$groundTruth_observer2[only2])
})

test_that("errors clearly on fewer than two detectors or mismatched observerSuffix length", {
  d <- make_n_detectors(2)
  expect_error(simsTocaptureHistoryTable(d[[1]]), "at least two detectors")
  expect_error(simsTocaptureHistoryTable(d[[1]], d[[2]], observerSuffix = "onlyone"),
              "one entry per")
})
