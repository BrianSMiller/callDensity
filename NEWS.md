# callDensity 1.1.0

## Bug fixes

* **`pa_CV()` computed the wrong standard error, inflating `CV.pa` and
  therefore `CV.Dc` and every confidence interval derived from it.**
  Between 2025-06-27 and this release, `se_pa` was computed as
  `overall_weighted_mean_pa / sqrt(no.transects)` -- the point estimate
  divided by a constant, not a standard error of anything. With the
  package's typical transect count this pins `CV.pa` within a whisker of
  `1/sqrt(no.transects)` almost regardless of the underlying data. Restored
  to the correct standard error of the transect means, generalised to
  unequal transect weights via `Hmisc::wtd.var` (reduces exactly to the
  original unweighted formula when weights are uniform). Verified by
  mutation testing and against real study data (16 site/observer/detector
  combinations): the buggy `CV.pa` ran 1.5x to 4.6x the corrected value,
  mean 2.3x. **Anyone who has used `cde()` and reported `CV.pa`, `CV.Dc`,
  or confidence intervals on `Dc` since 2025-06-27 should recompute them
  with this version.**

* **SNR truncation (`snrTruncationThreshold`) set detection probability to
  `NA` below the threshold; it should be `0`.** `NA` drops a cell from both
  the numerator and denominator of `pDetInArea`'s range-weighted mean, so
  `p_a` was computed as `E[p | SNR >= theta]` and `cde()` returned the
  density of *above-threshold* calls, silently scaled down from the density
  of all calls by a factor of `q(theta)` (the fraction of calls arriving
  above the threshold). Zeroing instead gives `E[p * 1(SNR >= theta)]`, so
  `q(theta)` cancels against a correspondingly-truncated `Nc` and `cde()`
  now returns the density of all calls, matching its documented behaviour.
  Distance truncation is unaffected -- it still uses `NA`, correctly, since
  cells beyond the truncation distance genuinely leave the study area
  rather than merely becoming undetectable.

* **`countDetections()` hardcoded `det$snr`** despite taking an
  `snrColName` argument, so truncating with any other column name (e.g.
  `"SNR"`, which is what `capHistTab` actually uses) silently returned
  `Nc = 0` with no error. Fixed to use the column the caller names.

* **`cde()` called the deprecated `fitSNRdetectionFunc()` internally**,
  which was already a pure passthrough to `fitDetFun()`. Every user who let
  `cde()` fit its own detection curve got a spurious "`fitSNRdetectionFunc`
  is deprecated" warning for a function they never called. `cde()` now
  calls `fitDetFun()` directly.

## New features

* **SNR truncation is now fully supported**, not experimental. `cde()`
  gains `NcIsTruncated`, a required confirmation whenever
  `snrTruncationThreshold` is set: `cde()` receives `Nc` as a plain number
  and never sees the detections it was counted from, so it cannot truncate
  or verify it itself. An untruncated `Nc` paired with a truncated `p_a`
  would silently inflate `Dc` by roughly `1/q(theta)`. If `snrDetFun` is
  supplied alongside a threshold, `cde()` now also messages a reminder that
  the curve should be fitted on the truncated sample, since it cannot refit
  a model it did not fit.

* **`nlFromDetections()` is now `cde()`'s default noise-level estimator**,
  replacing `nlFromSnrInfo()`. Both correct the same bias (noise levels
  measured at detections are biased low, because detections over-represent
  quiet periods), but `nlFromSnrInfo()` did so by adding back the SNR at
  which the detection function reaches 0.5 -- a property of the detector --
  while the bias itself is a property of the propagation geometry and the
  noise variance. The two coincide only by chance. `nlFromDetections()`
  inverts the actual detection-probability weighting instead, so it
  generalises across detector shapes. Works with `glm`/`gam`/`scam` and
  `vglm` detection functions alike, via the package's own dispatch-aware
  `predictDetFun()`. Pass `NL` explicitly to use a different estimate, as
  before.

* `callDensity_snrThreshold` vignette completed: the `q(theta)` cancellation
  that makes truncation return all-call density, and a `recaptureCoverage()`
  diagnostic for choosing the threshold from where recaptures (not raw
  sample coverage) run out.

## Deprecated

* **`cdeFromParamFile()` is deprecated.** It had drifted into a
  near-duplicate of `cde()`'s body: it called `pa_CV()` without transect
  area weights (`cde()` always weights by `truncationDistance^2`) and had
  no equivalent of `cde()`'s `NcIsTruncated` guard. It is now a thin
  wrapper that derives `cde()`'s arguments from its parameter object and
  calls `cde()` directly, so it inherits every fix above rather than
  carrying its own, separately-unmaintained copy of the calculation. Use
  `cde()` directly going forward.

## Other

* `future.apply`, used directly by `pDetInArea(parallel = TRUE)`, is now
  declared in `Suggests` (previously an undeclared dependency).

# callDensity 1.0.1

Not tagged as a release at the time. Folded into this NEWS file
retroactively for a complete version history; see git history for exact
commits.

## New features

* New unified detection-function interface: `fitDetFun()`, `predictDetFun()`,
  `showDetFun()`, replacing `fitSNRdetectionFunc()`, `fitSNRvglm()`, and
  `fitSNRvgam()` (retained as deprecated wrappers). Supports `glm`, `gam`,
  `scam`, and `vglm` (capture-recapture) models behind one interface, with
  proper S3 dispatch for prediction and plotting -- including for `vglm`
  objects, which `stats::predict()` cannot handle directly.
* `showDetFun()`: plots a single model (with SNR-info overlay, rug, optional
  mirrored density/histogram) or a named list of models for multi-curve
  comparison.
* `pDetInArea()`: parallel execution via `future.apply`; `truncationDistance`
  default fixed to reference `TLlookup`; `snrTruncationThreshold` argument
  added (see 1.1.0 for the fix to how it's applied).

## Deprecated

* `fitSNRdetectionFunc()`, `fitSNRvglm()`, `fitSNRvgam()`, and the old
  `predict`/`show` functions superseded by the unified interface above.

## Testing

* New testthat suite (79 tests) for the detection-function interface.

# callDensity 1.0.0

Initial public release.
