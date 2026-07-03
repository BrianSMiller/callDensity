# tests/testthat/test-pDetInArea.R
#
# TL, SL, NL, and truncation distance all match the vignette.
# outerloop is reduced to 20-30 for speed.

run_pdet <- function(fit, outerloop = 20, truncDist = 1e6, ...) {
  pDetInArea(
    detFun             = fit,
    SL                 = make_sl(),
    TLlookup           = make_tl_lookup(),
    NL                 = make_nl(),
    outerloop          = outerloop,
    truncationDistance = truncDist,
    ...
  )
}

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("pDetInArea returns a list with the expected elements", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  res <- run_pdet(fit)
  expect_type(res, "list")
  expect_named(res, c("overall", "perTransectMeanSD", "meanOfAllTransects"),
               ignore.order = TRUE)
})

test_that("pDetInArea overall pDet is a single finite numeric in [0, 1]", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  res <- run_pdet(fit)
  expect_length(res$overall, 1L)
  expect_true(is.finite(res$overall))
  expect_gte(res$overall, 0)
  expect_lte(res$overall, 1)
})

test_that("pDetInArea perTransectMeanSD has correct dimensions", {
  # make_tl_lookup() produces 4 transects (matching vignette)
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  res <- run_pdet(fit)
  expect_equal(nrow(res$perTransectMeanSD), 4 + 1L)
  expect_equal(ncol(res$perTransectMeanSD), 2L)
})

test_that("pDetInArea meanOfAllTransects has range_m and pDet columns", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "gam")
  res <- run_pdet(fit)
  expect_s3_class(res$meanOfAllTransects, "data.frame")
  expect_true(all(c("range_m", "pDet") %in% names(res$meanOfAllTransects)))
})

test_that("pDetInArea pDet declines with range (broadly)", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  res <- run_pdet(fit, outerloop = 30)
  mat <- res$meanOfAllTransects
  near <- mean(mat$pDet[mat$range_m <= 1e4],  na.rm = TRUE)
  far  <- mean(mat$pDet[mat$range_m >= 5e5], na.rm = TRUE)
  expect_gt(near, far)
})

# ---------------------------------------------------------------------------
# Works for all supported model types (serial)
# ---------------------------------------------------------------------------

test_that("pDetInArea runs without error for glm", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  expect_no_error(run_pdet(fit))
})

test_that("pDetInArea runs without error for gam", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "gam")
  expect_no_error(run_pdet(fit))
})

test_that("pDetInArea runs without error for scam", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  expect_no_error(run_pdet(fit))
})

# ---------------------------------------------------------------------------
# Truncation distance
# ---------------------------------------------------------------------------

test_that("pDetInArea pDet is NA beyond truncation distance", {
  d       <- make_snr_data()
  fit     <- fitDetFun(d, modelType = "glm")
  trunc_m <- 1e5
  res     <- run_pdet(fit, truncDist = trunc_m)
  mat     <- res$meanOfAllTransects
  expect_true(all(is.na(mat$pDet[mat$range_m > trunc_m])))
})

# ---------------------------------------------------------------------------
# SNR truncation threshold
# ---------------------------------------------------------------------------

test_that("snrTruncationThreshold sets pDet to NA where SNR is below threshold", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  # Use a very high threshold so it bites across the full TL range
  res <- run_pdet(fit, outerloop = 20, snrTruncationThreshold = 1e6)
  # With threshold = 1e6 dB, no real SNR can exceed it, so all pDet should be NA
  expect_true(all(is.na(res$meanOfAllTransects$pDet)))
})

# ---------------------------------------------------------------------------
# TL default argument -- now fixed
# ---------------------------------------------------------------------------

test_that("pDetInArea default truncationDistance references TLlookup not TL", {
  default_expr <- formals(pDetInArea)$truncationDistance
  expect_true(grepl("TLlookup", deparse(default_expr)))
})

# ---------------------------------------------------------------------------
# parallel argument
# ---------------------------------------------------------------------------

test_that("pDetInArea parallel = FALSE produces a valid result", {
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  res <- run_pdet(fit)
  expect_true(is.finite(res$overall))
})

# ---------------------------------------------------------------------------
# Known bug: parallel path fails for all non-vglm models.
# The model slimming code strips slots that predict() needs in worker
# processes. When fixed, change expect_error() to expect_no_error().
# ---------------------------------------------------------------------------

test_that("pDetInArea parallel = TRUE with glm [known bug: slimming]", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  expect_error(
    suppressWarnings(suppressMessages(run_pdet(fit, outerloop = 10, parallel = TRUE)))
  )
})

test_that("pDetInArea parallel = TRUE with scam", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "scam", numKnots = 5)
  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  res <- suppressWarnings(suppressMessages(run_pdet(fit, outerloop = 10, parallel = TRUE)))
  expect_true(is.finite(res$overall))
  expect_gte(res$overall, 0)
  expect_lte(res$overall, 1)
})

test_that("pDetInArea emits a message for parallel non-vglm models", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  d   <- make_snr_data()
  fit <- fitDetFun(d, modelType = "glm")
  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  expect_message(
    try(suppressWarnings(run_pdet(fit, outerloop = 10, parallel = TRUE)),
        silent = TRUE),
    regexp = "parallel=TRUE with non-VGLM"
  )
})
