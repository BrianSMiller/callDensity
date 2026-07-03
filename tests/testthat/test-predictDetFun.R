# tests/testthat/test-predictDetFun.R

# ---------------------------------------------------------------------------
# predictDetFun -- single model
# ---------------------------------------------------------------------------

test_that("predictDetFun returns a data frame with SNR and fit columns", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  p   <- predictDetFun(fit)
  expect_s3_class(p, "data.frame")
  expect_true(all(c("SNR", "fit") %in% names(p)))
})

test_that("predictDetFun fit values are in [0, 1]", {
  d <- make_snr_data()
  for (mt in c("glm", "gam")) {
    p <- predictDetFun(fitDetFun(d, modelType = mt))
    expect_true(all(p$fit >= 0 & p$fit <= 1), label = mt)
  }
  p <- predictDetFun(fitDetFun(d, modelType = "scam", numKnots = 5))
  expect_true(all(p$fit >= 0 & p$fit <= 1), label = "scam")
})

test_that("predictDetFun returns CI columns when ci = TRUE", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  p   <- predictDetFun(fit, ci = TRUE)
  expect_named(p, c("SNR", "fit", "lower", "upper"), ignore.order = TRUE)
  expect_true(all(p$lower <= p$fit))
  expect_true(all(p$fit   <= p$upper))
})

test_that("predictDetFun omits CI columns when ci = FALSE", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  p   <- predictDetFun(fit, ci = FALSE)
  expect_false("lower" %in% names(p))
  expect_false("upper" %in% names(p))
})

test_that("predictDetFun respects npoints argument", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "gam")
  p   <- predictDetFun(fit, npoints = 50, ci = FALSE)
  expect_equal(nrow(p), 50)
})

test_that("predictDetFun respects newdata argument", {
  d      <- make_snr_data()
  fit    <- fitDetFun(d, modelType = "glm")
  newdat <- data.frame(SNR = c(-5, 0, 5, 10, 15))
  p      <- predictDetFun(fit, newdata = newdat, ci = FALSE)
  expect_equal(nrow(p), 5)
  expect_equal(p$SNR, newdat$SNR)
})

test_that("predictDetFun SNR grid spans the training data range", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  p   <- predictDetFun(fit, ci = FALSE)
  expect_gte(min(p$SNR), min(d$SNR) - 1e-9)
  expect_lte(max(p$SNR), max(d$SNR) + 1e-9)
})

test_that("predictDetFun.default errors informatively for unknown class", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  class(fit) <- "unknownClass"
  expect_error(predictDetFun(fit), regexp = "No predictDetFun\\(\\) method")
})

# ---------------------------------------------------------------------------
# predictDetFunList -- multi-model
# ---------------------------------------------------------------------------

test_that("predictDetFunList returns long-format data frame with model column", {
  d      <- make_snr_data()
  models <- list(
    glm  = fitDetFun(d, modelType = "glm"),
    scam = fitDetFun(d, modelType = "scam", numKnots = 5)
  )
  p <- predictDetFunList(models, ci = FALSE)
  expect_s3_class(p, "data.frame")
  expect_true("model" %in% names(p))
  expect_setequal(unique(p$model), c("glm", "scam"))
})

test_that("predictDetFunList uses a shared SNR grid across models", {
  d      <- make_snr_data()
  models <- list(
    glm = fitDetFun(d, modelType = "glm"),
    gam = fitDetFun(d, modelType = "gam")
  )
  p <- predictDetFunList(models, npoints = 100, ci = FALSE)
  expect_equal(as.integer(table(p$model)[["glm"]]), 100L)
  expect_equal(as.integer(table(p$model)[["gam"]]), 100L)
  expect_equal(p$SNR[p$model == "glm"], p$SNR[p$model == "gam"])
})

# ---------------------------------------------------------------------------
# Known bug: predictDetFunList calls model.frame() which fails for vglm.
# When fixed, change expect_error() to expect_no_error().
# ---------------------------------------------------------------------------
test_that("predictDetFunList with vglm and no newdata [known bug: model.frame fails]", {
  skip_if_not_installed("VGAM")
  d <- make_two_observer_data()
  fit_vglm <- suppressWarnings(
    fitDetFun(d, modelType = "vglm",
              yColNames = c("detect_observer1", "detect_observer2"))
  )
  expect_error(predictDetFunList(list(vglm = fit_vglm), ci = FALSE))
})
