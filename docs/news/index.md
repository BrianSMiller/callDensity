# Changelog

## callDensity 1.2.0

### New features

- **[`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md)
  replaces the three-function
  [`mchToCR()`](https://briansmiller.github.io/callDensity/reference/mchToCR.md)
  /
  [`capHist2snrInfo()`](https://briansmiller.github.io/callDensity/reference/capHist2snrInfo.md)
  /
  [`capHistTosnrInfo()`](https://briansmiller.github.io/callDensity/reference/capHistTosnrInfo.md)
  sprawl with one converter, and adds genuine N-observer
  capture-recapture support that no function in the package had
  before.** These three grew independently as the package’s capabilities
  expanded, without ever settling on one contract: two assume a
  2-observer `detect_table1`/`detect_table2` shape that requires
  renaming native matchbox columns before use;
  [`capHistTosnrInfo()`](https://briansmiller.github.io/callDensity/reference/capHistTosnrInfo.md)
  additionally assumes SNR was measured *before* matching (a
  per-observer `signalRMSdB1`/`signalRMSdB2` pair, averaged) rather than
  *after* (matchbox-native: one shared value per matched event); and
  none of the three can actually feed `fitDetFun(modelType = 'vglm')`,
  despite
  [`capHistTosnrInfo()`](https://briansmiller.github.io/callDensity/reference/capHistTosnrInfo.md)’s
  own documentation saying it could – `vglm`’s `posbernoulli.t` needs
  the raw per-observer 0/1 columns still present in the table, and both
  functions collapse them into a single `Detected` column before
  returning. Every real capture-recapture analysis using this package
  has had to hand-build its VGLM input as a result.

  `chtToSNRinfo(cht, groundTruth, observers, ...)` fixes all three
  problems at once. It defaults to native matchbox column names (`t0`,
  `signalRMSdB`, `noiseRMSdB`, `verdict`, `detect_observerN`), accepts
  any number of `observers` (not just two –
  [`VGAM::posbernoulli.t`](https://rdrr.io/pkg/VGAM/man/posbernoulli.t.html)
  was never actually limited to two occasions, it just hadn’t been
  exercised with more in this package before now), and keeps the raw
  observer columns intact in its output so the result can be passed
  straight to
  `fitDetFun(modelType = 'vglm', yColNames = c(groundTruthCol, observerCols))`
  with nothing hand-built. A single-element `observers` reproduces the
  old 2-observer OG shape exactly, so this isn’t a second code path
  alongside a simpler one – OG is just the N=1 case of the same
  function.

- **[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  gains `groundTruthCol`/`observerCol` arguments**, forwarded to both
  [`falseDiscoveryRate()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryRate.md)
  (which already accepted them, but couldn’t be reached –
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  called it positionally without forwarding) and the new
  [`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md).
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  can now run directly against a native matchbox `capHistTab`
  (`verdict`, `detect_observerN`) with no reduction step, in addition to
  the legacy `detect_table1`/`detect_table2` shape it already supported.
  Defaults are unchanged (`'detect_table1'`/`'detect_table2'`), so no
  existing call site is affected.

- **[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  also gains `snrColName`/`timeCol` arguments**, same reasoning as
  `groundTruthCol`/`observerCol` above and found the same way – by
  actually running
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  against raw simulation data. `snrColName` (default `'SNR'`) is
  forwarded to
  [`falseDiscoveryRate()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryRate.md),
  which already accepted it but wasn’t reachable from
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md);
  raw simulation/detection tables typically carry lowercase `snr`
  instead. `timeCol` (default `NULL`, auto-detecting `'t0'` then `'t'`
  as before) is forwarded to the internal
  [`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md)
  call, for tables using neither (e.g. `'datetime'`).

- **[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  gains opt-in `returnDetFun`/`returnPdetDetail` arguments**, both
  default `FALSE`. When requested, the detection function actually used
  (whether passed in or fit internally) and/or
  [`pDetInArea()`](https://briansmiller.github.io/callDensity/reference/pDetInArea.md)’s
  full result (including the per-transect detail needed for a
  probability-of-detection-vs-range plot) are attached as attributes on
  the returned data.frame (`attr(result, "detFun")`/`"pDetResults"`)
  rather than changing its shape – `result` is always the same plain
  data.frame either way, so
  [`rbind()`](https://rdrr.io/r/base/cbind.html), `result$Dc`, etc. all
  keep working unchanged whether or not these are requested.

- **[`pDetInArea()`](https://briansmiller.github.io/callDensity/reference/pDetInArea.md)
  now returns `allDetFunctions`**, the full range-by-transect
  probability-of-detection grid that was already being computed
  internally but previously only escaped via the optional `transectFile`
  CSV write. Purely additive (a new named list element, not a
  fixed-shape data.frame), so existing callers accessing
  `$overall`/`$perTransectMeanSD`/`$meanOfAllTransects` are unaffected.

- **New
  [`plotPDetRadials()`](https://briansmiller.github.io/callDensity/reference/plotPDetRadials.md)**
  plots probability of detection as a function of both range and azimuth
  around the recorder – a directional detection-range footprint, rather
  than the 1D range-only view `meanOfAllTransects` gives. Ports Brian’s
  MATLAB `plotPDetRadials.m` to a native ggplot2
  [`coord_polar()`](https://ggplot2.tidyverse.org/reference/coord_radial.html)
  plot. Takes
  [`pDetInArea()`](https://briansmiller.github.io/callDensity/reference/pDetInArea.md)’s
  own `allDetFunctions` output directly.

- **Vignettes reworked to use existing package functions instead of
  hand-rolled equivalents**, closing a gap where the package’s own
  capabilities
  ([`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md),
  [`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md),
  [`showDetFun()`](https://briansmiller.github.io/callDensity/reference/showDetFun.md),
  [`predictDetFun()`](https://briansmiller.github.io/callDensity/reference/predictDetFun.md),
  [`plotSpatialDetections()`](https://briansmiller.github.io/callDensity/reference/plotSpatialDetections.md),
  and now
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  new opt-in returns) had grown without the published examples being
  updated to demonstrate them. `callDensity.Rmd` – the flagship
  vignette, same name as the package – previously reimplemented
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  own density formula by hand, three times over (once per
  detection-function model type), silently dropping every uncertainty
  column (`CV.Dc`, `CV.pa`, `CV.c`)
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  computes automatically; now uses
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  directly, and gained two diagnostics it never had before
  ([`showDetFun()`](https://briansmiller.github.io/callDensity/reference/showDetFun.md)’s
  fitted- curve-vs-observed-SNR plot, and
  [`plotPDetRadials()`](https://briansmiller.github.io/callDensity/reference/plotPDetRadials.md)’s
  directional footprint) rather than only ever plotting probability of
  detection vs. range. `parallel_benchmarks.Rmd`,
  `callDensity_CommonGround.Rmd`, `callDensity_coast.Rmd`, and
  `callDensity_snrThreshold.Rmd` updated similarly;
  `callDensity_CommonGround.Rmd`’s capture-recapture section had its own
  `TODO: ... make this an explicit parameter` comment, directly resolved
  by
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  new `groundTruthCol`/`observerCol`.

- **[`simsTocaptureHistoryTable()`](https://briansmiller.github.io/callDensity/reference/simsTocaptureHistoryTable.md)
  defaults to matchbox-native `detect_observer1`/`detect_observer2`**
  column naming (was `detect_table1`/`detect_table2`). New
  `observerSuffix` argument, default `c('observer1','observer2')`; pass
  `c('table1','table2')` for the previous naming. Every vignette’s
  simulated example now demonstrates the same column convention
  [`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md)
  defaults to, rather than teaching the older convention by example
  while the documented default points elsewhere.

### Deprecated

- **[`mchToCR()`](https://briansmiller.github.io/callDensity/reference/mchToCR.md),
  [`capHist2snrInfo()`](https://briansmiller.github.io/callDensity/reference/capHist2snrInfo.md),
  and
  [`capHistTosnrInfo()`](https://briansmiller.github.io/callDensity/reference/capHistTosnrInfo.md)
  are deprecated in favour of
  [`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md).**
  All three still work exactly as before – bodies unchanged
  ([`mchToCR()`](https://briansmiller.github.io/callDensity/reference/mchToCR.md))
  or reimplemented as thin wrappers with identical output
  ([`capHist2snrInfo()`](https://briansmiller.github.io/callDensity/reference/capHist2snrInfo.md),
  [`capHistTosnrInfo()`](https://briansmiller.github.io/callDensity/reference/capHistTosnrInfo.md))
  – so the published Common Ground and Beyond Counting Calls analysis
  scripts keep reproducing against current `main` without modification.
  New analyses should call
  [`chtToSNRinfo()`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md)
  directly; none of the three deprecated functions’
  2-observer-reduction/renaming machinery is needed by any current
  workflow, since neither `fitDetFun(modelType = 'vglm')` nor (as of
  this release)
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  requires it.

  This is a deliberate streamlining, not a stopgap: going forward, one
  function covers OG and CR alike, for any number of observers, on
  native matchbox column names by default – with the old interfaces kept
  working, not removed, so nothing currently in production or published
  breaks.

## callDensity 1.1.1

Reworks the callDensity_snrThreshold vignette around the actual
failing  
scenario the truncation feature addresses, rather than an arbitrary  
example. Drops the observer-ground (OG) design from that vignette  
entirely (covered elsewhere, in callDensity_CommonGround); fixes a
detection-function  
plot that was showing the truncated fit extrapolating past its own  
support while hiding true ground truth in the same region; and adjusts  
the false-positive simulation parameters to match the qualitative
pattern  
in Miller et al. (2026) Figure 3, with an explicit note that the  
resulting truncation-lowers-c effect is not a general result.

## callDensity 1.1.0 (unreleased)

### Bug fixes

- **[`pa_CV()`](https://briansmiller.github.io/callDensity/reference/pa_CV.md)
  computed the wrong standard error, inflating `CV.pa` and therefore
  `CV.Dc` and every confidence interval derived from it.** Between
  2025-06-27 and this release, `se_pa` was computed as
  `overall_weighted_mean_pa / sqrt(no.transects)` – the point estimate
  divided by a constant, not a standard error of anything. With the
  package’s typical transect count this pins `CV.pa` within a whisker of
  `1/sqrt(no.transects)` almost regardless of the underlying data.
  Restored to the correct standard error of the transect means,
  generalised to unequal transect weights via
  [`Hmisc::wtd.var`](https://rdrr.io/pkg/Hmisc/man/wtd.stats.html)
  (reduces exactly to the original unweighted formula when weights are
  uniform). Verified by mutation testing and against real study data (16
  site/observer/detector combinations): the buggy `CV.pa` ran 1.5x to
  4.6x the corrected value, mean 2.3x. **Anyone who has used
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  and reported `CV.pa`, `CV.Dc`, or confidence intervals on `Dc` since
  2025-06-27 should recompute them with this version.**

- **SNR truncation (`snrTruncationThreshold`) set detection probability
  to `NA` below the threshold; it should be `0`.** `NA` drops a cell
  from both the numerator and denominator of `pDetInArea`’s
  range-weighted mean, so `p_a` was computed as `E[p | SNR >= theta]`
  and
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  returned the density of *above-threshold* calls, silently scaled down
  from the density of all calls by a factor of `q(theta)` (the fraction
  of calls arriving above the threshold). Zeroing instead gives
  `E[p * 1(SNR >= theta)]`, so `q(theta)` cancels against a
  correspondingly-truncated `Nc` and
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  now returns the density of all calls, matching its documented
  behaviour. Distance truncation is unaffected – it still uses `NA`,
  correctly, since cells beyond the truncation distance genuinely leave
  the study area rather than merely becoming undetectable.

- **[`countDetections()`](https://briansmiller.github.io/callDensity/reference/countDetections.md)
  hardcoded `det$snr`** despite taking an `snrColName` argument, so
  truncating with any other column name (e.g. `"SNR"`, which is what
  `capHistTab` actually uses) silently returned `Nc = 0` with no error.
  Fixed to use the column the caller names.

- **[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  called the deprecated
  [`fitSNRdetectionFunc()`](https://briansmiller.github.io/callDensity/reference/fitSNRdetectionFunc.md)
  internally**, which was already a pure passthrough to
  [`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md).
  Every user who let
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  fit its own detection curve got a spurious “`fitSNRdetectionFunc` is
  deprecated” warning for a function they never called.
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  now calls
  [`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md)
  directly.

### New features

- **SNR truncation is now fully supported**, not experimental.
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  gains `NcIsTruncated`, a required confirmation whenever
  `snrTruncationThreshold` is set:
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  receives `Nc` as a plain number and never sees the detections it was
  counted from, so it cannot truncate or verify it itself. An
  untruncated `Nc` paired with a truncated `p_a` would silently inflate
  `Dc` by roughly `1/q(theta)`. If `snrDetFun` is supplied alongside a
  threshold,
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  now also messages a reminder that the curve should be fitted on the
  truncated sample, since it cannot refit a model it did not fit.

- **[`nlFromDetections()`](https://briansmiller.github.io/callDensity/reference/nlFromDetections.md)
  is now
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  default noise-level estimator**, replacing
  [`nlFromSnrInfo()`](https://briansmiller.github.io/callDensity/reference/nlFromSnrInfo.md).
  Both correct the same bias (noise levels measured at detections are
  biased low, because detections over-represent quiet periods), but
  [`nlFromSnrInfo()`](https://briansmiller.github.io/callDensity/reference/nlFromSnrInfo.md)
  did so by adding back the SNR at which the detection function reaches
  0.5 – a property of the detector – while the bias itself is a property
  of the propagation geometry and the noise variance. The two coincide
  only by chance.
  [`nlFromDetections()`](https://briansmiller.github.io/callDensity/reference/nlFromDetections.md)
  inverts the actual detection-probability weighting instead, so it
  generalises across detector shapes. Works with `glm`/`gam`/`scam` and
  `vglm` detection functions alike, via the package’s own dispatch-aware
  [`predictDetFun()`](https://briansmiller.github.io/callDensity/reference/predictDetFun.md).
  Pass `NL` explicitly to use a different estimate, as before.

- `callDensity_snrThreshold` vignette completed: the `q(theta)`
  cancellation that makes truncation return all-call density, and a
  `recaptureCoverage()` diagnostic for choosing the threshold from where
  recaptures (not raw sample coverage) run out.

### Deprecated

- **[`cdeFromParamFile()`](https://briansmiller.github.io/callDensity/reference/cdeFromParamFile.md)
  is deprecated.** It had drifted into a near-duplicate of
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  body: it called
  [`pa_CV()`](https://briansmiller.github.io/callDensity/reference/pa_CV.md)
  without transect area weights
  ([`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  always weights by `truncationDistance^2`) and had no equivalent of
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  `NcIsTruncated` guard. It is now a thin wrapper that derives
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)’s
  arguments from its parameter object and calls
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  directly, so it inherits every fix above rather than carrying its own,
  separately-unmaintained copy of the calculation. Use
  [`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
  directly going forward.

### Other

- `future.apply`, used directly by `pDetInArea(parallel = TRUE)`, is now
  declared in `Suggests` (previously an undeclared dependency).

## callDensity 1.0.1

Not tagged as a release at the time. Folded into this NEWS file
retroactively for a complete version history; see git history for exact
commits.

### New features

- New unified detection-function interface:
  [`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md),
  [`predictDetFun()`](https://briansmiller.github.io/callDensity/reference/predictDetFun.md),
  [`showDetFun()`](https://briansmiller.github.io/callDensity/reference/showDetFun.md),
  replacing
  [`fitSNRdetectionFunc()`](https://briansmiller.github.io/callDensity/reference/fitSNRdetectionFunc.md),
  [`fitSNRvglm()`](https://briansmiller.github.io/callDensity/reference/fitSNRvglm.md),
  and
  [`fitSNRvgam()`](https://briansmiller.github.io/callDensity/reference/fitSNRvgam.md)
  (retained as deprecated wrappers). Supports `glm`, `gam`, `scam`, and
  `vglm` (capture-recapture) models behind one interface, with proper S3
  dispatch for prediction and plotting – including for `vglm` objects,
  which [`stats::predict()`](https://rdrr.io/r/stats/predict.html)
  cannot handle directly.
- [`showDetFun()`](https://briansmiller.github.io/callDensity/reference/showDetFun.md):
  plots a single model (with SNR-info overlay, rug, optional mirrored
  density/histogram) or a named list of models for multi-curve
  comparison.
- [`pDetInArea()`](https://briansmiller.github.io/callDensity/reference/pDetInArea.md):
  parallel execution via `future.apply`; `truncationDistance` default
  fixed to reference `TLlookup`; `snrTruncationThreshold` argument added
  (see 1.1.0 for the fix to how it’s applied).

### Deprecated

- [`fitSNRdetectionFunc()`](https://briansmiller.github.io/callDensity/reference/fitSNRdetectionFunc.md),
  [`fitSNRvglm()`](https://briansmiller.github.io/callDensity/reference/fitSNRvglm.md),
  [`fitSNRvgam()`](https://briansmiller.github.io/callDensity/reference/fitSNRvgam.md),
  and the old `predict`/`show` functions superseded by the unified
  interface above.

### Testing

- New testthat suite (79 tests) for the detection-function interface.

## callDensity 1.0.0

Initial public release.
