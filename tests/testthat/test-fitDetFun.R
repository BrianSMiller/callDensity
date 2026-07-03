# tests/testthat/test-fitDetFun.R

test_that("fitDetFun returns 'detFun' as first class for scam", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  expect_s3_class(fit, "detFun")
  expect_equal(class(fit)[1], "detFun")
})

test_that("fitDetFun preserves the underlying model class for scam", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  expect_s3_class(fit, "scam")
})

test_that("fitDetFun preserves the underlying model class for gam", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "gam")
  expect_s3_class(fit, "detFun")
  expect_s3_class(fit, "gam")
})

test_that("fitDetFun preserves the underlying model class for glm", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  expect_s3_class(fit, "detFun")
  expect_s3_class(fit, "glm")
})

test_that("fitDetFun sets modelType attribute correctly", {
  d <- make_snr_data()
  expect_equal(attr(fitDetFun(d, modelType = "glm"),  "modelType"), "glm")
  expect_equal(attr(fitDetFun(d, modelType = "gam"),  "modelType"), "gam")
  expect_equal(attr(fitDetFun(d, modelType = "scam", numKnots = 5), "modelType"), "scam")
})

test_that("fitDetFun rejects invalid modelType", {
  d <- make_snr_data()
  expect_error(fitDetFun(d, modelType = "loess"), "arg")
})

test_that("fitDetFun fitted values are in [0, 1]", {
  d <- make_snr_data()
  expect_true(all(fitted(fitDetFun(d, modelType = "glm"))  >= 0))
  expect_true(all(fitted(fitDetFun(d, modelType = "gam"))  >= 0))
  expect_true(all(fitted(fitDetFun(d, modelType = "scam", numKnots = 5)) >= 0))
  expect_true(all(fitted(fitDetFun(d, modelType = "glm"))  <= 1))
  expect_true(all(fitted(fitDetFun(d, modelType = "gam"))  <= 1))
  expect_true(all(fitted(fitDetFun(d, modelType = "scam", numKnots = 5)) <= 1))
})

test_that("coef() and vcov() work on glm detFun objects", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))
})

test_that("coef() and Vp slot work on scam detFun objects", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  expect_true(is.numeric(coef(fit)))
  expect_true(!is.null(fit$Vp))
})

# ---------------------------------------------------------------------------
# Deprecated wrappers
# ---------------------------------------------------------------------------

test_that("fitSNRdetectionFunc emits a deprecation warning", {
  d <- make_snr_data()
  expect_warning(fitSNRdetectionFunc(d, modelType = "glm"),
                 regexp = "deprecated", ignore.case = TRUE)
})

test_that("fitSNRdetectionFunc returns the same result as fitDetFun", {
  d        <- make_snr_data()
  fit_new  <- fitDetFun(d, modelType = "glm")
  fit_depr <- suppressWarnings(fitSNRdetectionFunc(d, modelType = "glm"))
  expect_s3_class(fit_depr, "detFun")
  expect_s3_class(fit_depr, "glm")
  expect_equal(coef(fit_depr), coef(fit_new))
})

test_that("fitSNRvglm emits a deprecation warning", {
  skip_if_not_installed("VGAM")
  d <- make_two_observer_data()
  expect_warning(fitSNRvglm(d), regexp = "deprecated", ignore.case = TRUE)
})

test_that("fitSNRvgam errors with a deprecation warning", {
  expect_warning(expect_error(fitSNRvgam()),
                 regexp = "deprecated", ignore.case = TRUE)
})
