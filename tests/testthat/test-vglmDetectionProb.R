# ---------------------------------------------------------------------------
# vglmDetectionProb() -- the consolidated vglm prediction helper used by both
# predictDetFun.vglm() and pDetInArea(). These tests focus on the
# whichObserver = "any" (union/onempall0) mode, and specifically on it
# generalizing correctly to more than two occasions -- the two-observer case
# alone doesn't distinguish "genuinely N-observer" from "happens to work for
# N=2".
# ---------------------------------------------------------------------------

test_that("whichObserver='any' matches hand-rolled onempall0, for N=2..5 occasions", {
  skip_if_not_installed("VGAM")

  newd <- data.frame(SNR = seq(-20, 20, length.out = 50))

  for (N in 2:5) {
    d <- make_n_observer_data(n_subsample = 600, n_observers = N, seed = 42)
    obsCols <- paste0("detect_observer", seq_len(N))
    fit <- suppressWarnings(fitDetFun(d, modelType = "vglm", yColNames = obsCols))

    p_any <- predictDetFun(fit, newdata = newd, whichObserver = "any")$fit
    hand_rolled <- as.numeric(VGAM::predictvglm(fit, type = "response", newdata = newd,
                                                type.fitted = "onempall0"))

    expect_equal(p_any, hand_rolled, info = paste("N =", N))
  }
})

test_that("whichObserver='any' is >= every individual occasion's own curve, for N=2..5", {
  skip_if_not_installed("VGAM")

  newd <- data.frame(SNR = seq(-20, 20, length.out = 50))

  for (N in 2:5) {
    d <- make_n_observer_data(n_subsample = 600, n_observers = N, seed = 42)
    obsCols <- paste0("detect_observer", seq_len(N))
    fit <- suppressWarnings(fitDetFun(d, modelType = "vglm", yColNames = obsCols))

    p_any  <- predictDetFun(fit, newdata = newd, whichObserver = "any")$fit
    p_each <- sapply(obsCols, function(col) {
      predictDetFun(fit, newdata = newd, whichObserver = col)$fit
    })

    expect_true(all(p_any >= p_each - 1e-9), info = paste("N =", N))
  }
})

test_that("onempall0 equals 1 - prod(1 - p_i) across occasions, for N=2..5", {
  skip_if_not_installed("VGAM")
  # This is the algebraic consequence of conditional independence given SNR,
  # which is what a posbernoulli.t model with no behavioural/individual
  # effects assumes. It isn't just a sanity check that the code runs -- it's
  # an independent reconstruction of the union probability from each
  # occasion's own marginal curve, and should hold to floating-point
  # precision if whichObserver="any" is doing what it claims.

  newd <- data.frame(SNR = seq(-20, 20, length.out = 50))

  for (N in 2:5) {
    d <- make_n_observer_data(n_subsample = 600, n_observers = N, seed = 42)
    obsCols <- paste0("detect_observer", seq_len(N))
    fit <- suppressWarnings(fitDetFun(d, modelType = "vglm", yColNames = obsCols))

    p_any  <- predictDetFun(fit, newdata = newd, whichObserver = "any")$fit
    p_each <- sapply(obsCols, function(col) {
      predictDetFun(fit, newdata = newd, whichObserver = col)$fit
    })
    reconstructed <- 1 - apply(1 - p_each, 1, prod)

    expect_equal(p_any, reconstructed, tolerance = 1e-8, info = paste("N =", N))
  }
})

test_that("pDetInArea() with whichObserver='any' produces a higher p_a than any individual detector", {
  skip_if_not_installed("VGAM")

  d <- make_n_observer_data(n_subsample = 600, n_observers = 3, seed = 42)
  obsCols <- paste0("detect_observer", 1:3)
  fit <- suppressWarnings(fitDetFun(d, modelType = "vglm", yColNames = obsCols))

  TL <- simTLradials_20logR(maxRange = 1e6, rangeStep = 5000, numTransects = 4)
  SL <- data.frame(mean = 190, sd = 4, sampleSize = 350)
  NL <- data.frame(mean = 84, sd = 4, sampleSize = 1000)

  set.seed(1)
  pa_individual <- pDetInArea(fit, SL, TL, NL, outerloop = 20, output.resolution.m = 5000)

  fit_union <- fit
  fit_union@extra$whichObserver <- "any"
  set.seed(1)
  pa_union <- pDetInArea(fit_union, SL, TL, NL, outerloop = 20, output.resolution.m = 5000)

  expect_true(pa_union$overall > pa_individual$overall)
  expect_true(pa_union$overall <= 1)
})

test_that("extractSNRinfo() handles whichObserver='any' (union mode), not just a single observer", {
  # Regression: extractSNRinfo() previously tried model@y[, "any"] directly
  # for a union-mode model, which errors ("subscript out of bounds") since
  # "any" isn't a real column of @y -- found while building
  # callDensity_unionDetectors.Rmd's showDetFun() comparison plot, where
  # predictDetFunList()'s default-range calculation calls extractSNRinfo()
  # on every model in a list, including a union one.
  skip_if_not_installed("VGAM")

  d <- make_n_observer_data(n_subsample = 600, n_observers = 3, seed = 42)
  obsCols <- paste0("detect_observer", 1:3)
  fit <- suppressWarnings(fitDetFun(d, modelType = "vglm", yColNames = obsCols))
  fitUnion <- fit
  fitUnion@extra$whichObserver <- "any"

  info <- expect_no_error(callDensity:::extractSNRinfo(fitUnion))
  expect_true(all(c("SNR", "Detected") %in% names(info)))

  # Detected should be the union across all three occasions, not any single
  # one -- so it must be >= each individual occasion's own Detected count,
  # and equal to the OR of all three.
  infoEach <- lapply(obsCols, function(col) {
    f <- fit
    f@extra$whichObserver <- col
    callDensity:::extractSNRinfo(f)
  })
  expectedUnion <- Reduce(`|`, lapply(infoEach, function(x) x$Detected))
  expect_equal(info$Detected, expectedUnion)
})
