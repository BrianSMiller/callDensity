# tests/testthat/test-showDetFun.R

test_that("showDetFun dispatches correctly on detFun objects", {
  skip_if_not_installed("ggplot2")
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  expect_s3_class(showDetFun(fit, SNRinfo = d), "ggplot")
})

test_that("showDetFun.list returns a ggplot for a named model list", {
  skip_if_not_installed("ggplot2")
  d      <- make_snr_data()
  models <- list(
    glm  = fitDetFun(d, modelType = "glm"),
    scam = fitDetFun(d, modelType = "scam", numKnots = 5)
  )
  expect_s3_class(showDetFun(models), "ggplot")
})

test_that("showDetFun.list colour aesthetic maps to model names", {
  skip_if_not_installed("ggplot2")
  d      <- make_snr_data()
  models <- list(
    glm = fitDetFun(d, modelType = "glm"),
    gam = fitDetFun(d, modelType = "gam")
  )
  p    <- showDetFun(models)
  data <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(length(unique(data$colour)), 2L)
})

test_that("showDetFun density distribution option works", {
  skip_if_not_installed("ggplot2")
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "gam")
  expect_no_error(showDetFun(fit, SNRinfo = d, distribution = "density"))
})

test_that("showDetFun histogram distribution option works", {
  skip_if_not_installed("ggplot2")
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  expect_no_error(showDetFun(fit, SNRinfo = d, distribution = "histogram"))
})
