# Tests for pa_CV.
#
# Context: pa_CV had no tests at all before this file. Between 2025-06-27 and
# 2026-07-17, se_pa was computed as overall_weighted_mean_pa/sqrt(no.transects)
# -- the ESTIMATE divided by a constant, not a standard error of anything. With
# 4 transects that pins CV.pa within a whisker of 1/sqrt(4) = 0.5 regardless of
# how right or wrong p_a actually is. Confirmed against the truncation sweep:
# true p_a bias ran 0% to 68% while the old formula reported CV.pa = 0.50-0.51
# throughout.
#
# These tests check the restored formula two ways: a deterministic synthetic
# case with a hand-computable answer (catches the exact bug, and would fail
# under the old formula -- verified by mutation), and a directional check that
# CV.pa responds to genuine between-transect disagreement, which is the whole
# point of the statistic.

test_that("pa_CV matches a hand-computed answer under equal weights", {
  # Four synthetic transects. Chosen so the arithmetic is checkable by hand
  # rather than merely internally consistent with the function under test.
  pa_t    <- c(0.10, 0.14, 0.09, 0.11)   # per-transect mean pDet
  sd_pa_t <- c(0.01, 0.02, 0.01, 0.02)   # per-transect SD (coefficient/SL/NL uncertainty)

  pa.all.transects <- data.frame(
    Mean = c(pa_t, weighted.mean(pa_t, rep(1, 4))),
    SD   = c(sd_pa_t, NA)  # last row's SD column is unused by pa_CV
  )

  result <- pa_CV(pa.all.transects, wt = rep(1, 4))

  # Hand-computed expectation: SE_Pa = sd(pa_t)/sqrt(n_t), the classical
  # standard error of the mean, NOT mean(pa_t)/sqrt(n_t).
  se_pa_expected  <- sd(pa_t) / sqrt(4)
  sd_term         <- mean(sd_pa_t)
  se_total        <- sqrt(se_pa_expected^2 + sd_term^2)
  expected        <- se_total / mean(pa_t)

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("pa_CV does not reduce to 1/sqrt(no.transects) regardless of the data", {
  # The signature of the bug: se_pa = mean/sqrt(n) makes CV.pa depend almost
  # entirely on no.transects and barely at all on pa_t or sd_pa_t. Two very
  # differently-behaved transect sets should give visibly different CV.pa.
  tight <- data.frame(Mean = c(0.10, 0.101, 0.099, 0.10, 0.10),
                      SD   = c(0.005, 0.005, 0.005, 0.005, NA))
  loose <- data.frame(Mean = c(0.05, 0.20, 0.03, 0.18, 0.09),
                      SD   = c(0.005, 0.005, 0.005, 0.005, NA))

  cvTight <- pa_CV(tight, wt = rep(1, 4))
  cvLoose <- pa_CV(loose, wt = rep(1, 4))

  # Both sets have the same per-transect SD (sd_pa_t) and a similar overall
  # mean; the only thing that differs is how much the transects disagree with
  # EACH OTHER. Under the old bug this made almost no difference (both would
  # sit near 0.5/mean-ish); under the restored formula, loose should be
  # substantially larger.
  expect_gt(cvLoose, cvTight * 3)
})

test_that("pa_CV is NOT the old mean/sqrt(n) formula (mutation guard)", {
  # Explicitly re-implements the 2025-06-27 to 2026-07-17 bug inline and
  # confirms the package's current pa_CV disagrees with it whenever the
  # transects genuinely disagree with each other. If someone reintroduces the
  # bug, this test starts passing when it shouldn't -- so the assertion is
  # framed as "these must differ", which only holds while the fix is in place.
  pa_t    <- c(0.06, 0.22, 0.05, 0.19)
  sd_pa_t <- c(0.01, 0.01, 0.01, 0.01)
  wt      <- rep(1, 4)

  pa.all.transects <- data.frame(Mean = c(pa_t, weighted.mean(pa_t, wt)),
                                 SD   = c(sd_pa_t, NA))

  current <- pa_CV(pa.all.transects, wt = wt)

  # The old, buggy formula, reproduced verbatim from commit 0a9cb2f.
  overall_weighted_mean_pa <- weighted.mean(pa_t, wt / sum(wt))
  se_pa_buggy <- overall_weighted_mean_pa / sqrt(4)
  overall_weighted_mean_SD <- weighted.mean(sd_pa_t, wt / sum(wt))
  buggy <- sqrt(se_pa_buggy^2 + overall_weighted_mean_SD^2) / overall_weighted_mean_pa

  expect_false(isTRUE(all.equal(current, buggy, tolerance = 0.05)))
})

test_that("pa_CV with uniform weights matches unweighted sd()", {
  # Hmisc::wtd.var(x, rep(1,n), normwt=TRUE) must reduce exactly to var(x).
  # This is the generalisation this fix relies on: it must not silently change
  # behaviour for the equal-transect-length case that Harris (2012) and both
  # published applications (BeyondCountingCalls, commonGround-final tags) used.
  set.seed(1)
  pa_t    <- runif(6, 0.05, 0.25)
  sd_pa_t <- runif(6, 0.005, 0.02)
  wt      <- rep(1, 6)

  pa.all.transects <- data.frame(Mean = c(pa_t, weighted.mean(pa_t, wt)),
                                 SD   = c(sd_pa_t, NA))

  result <- pa_CV(pa.all.transects, wt = wt)

  se_pa_expected <- sd(pa_t) / sqrt(6)
  expected <- sqrt(se_pa_expected^2 + mean(sd_pa_t)^2) / mean(pa_t)

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("unequal transect weights change CV.pa sensibly", {
  # Downweighting the transect that disagrees most with the others should
  # reduce the reported spread relative to equal weighting.
  pa_t    <- c(0.10, 0.10, 0.10, 0.40)  # transect 4 is an outlier
  sd_pa_t <- c(0.01, 0.01, 0.01, 0.01)

  equalWt <- rep(1, 4)
  downweightOutlier <- c(1, 1, 1, 0.05)

  paEqual <- data.frame(Mean = c(pa_t, weighted.mean(pa_t, equalWt)),
                        SD   = c(sd_pa_t, NA))
  paDown  <- data.frame(Mean = c(pa_t, weighted.mean(pa_t, downweightOutlier)),
                        SD   = c(sd_pa_t, NA))

  cvEqual <- pa_CV(paEqual, wt = equalWt)
  cvDown  <- pa_CV(paDown,  wt = downweightOutlier)

  expect_lt(cvDown, cvEqual)
})

test_that("a single transect gives a finite CV.pa with no between-transect term", {
  # Guards the no.transects <= 1 branch added to avoid Hmisc::wtd.var on a
  # single point (which would otherwise return NA/error).
  pa.all.transects <- data.frame(Mean = c(0.12, 0.12), SD = c(0.02, NA))

  result <- pa_CV(pa.all.transects, wt = 1)

  expect_true(is.finite(result))
  expect_equal(result, 0.02 / 0.12, tolerance = 1e-6)
})
