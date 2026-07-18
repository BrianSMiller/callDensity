# Call Density Estimate using auxiliary information from sonar equation

Applies the auxiliary information methods of Castro et. al.
(2024)/Harris (2012) to obtain call density, \\\hat{D_c}\\ using the
cannonical density equation: \$\$D_c =
\frac{N_c(1-c)}{kTA\hat{p}\_a}\$\$ and sonar equations \$\$SNR = RL -
NL\$\$ \$\$RL = SL - TL\$\$ where:

- \\\hat{D_c}\\ is call density

- \\N_c\\ is number of calls

- \\c\\ is false discovery rate

- \\k\\ is number of sensors (here always 1)

- \\T\\ is the duration of data analysed

- \\A\\ is the study area (in km\\2)

- \\\hat{p}\_a\\ is the probability of detection in the study area.

  
\\\hat{p}\_a\\ is estimated via `pDetInArea`, a Monte Carlo simulation
using auxiliary information for the sonar equation, so requires:

- SL - Distribution of source levels (mean, sd, samplesize)

- TL - Transmission loss model

- NL - Distribution of noise levels

- snrDetFun - a model of the relationship between SNR and probability of
  detection either via a capture history table or supplied detection
  function

## Usage

``` r
cde(
  Nc,
  capHistTab,
  SL,
  TL,
  T = 1,
  A = 1,
  k = 1,
  season = "year",
  snrDetFun = NULL,
  NL = NULL,
  modelType = "scam",
  numKnots = 5,
  output.resolution.m = 100,
  outerloop = 10,
  transectFile = NULL,
  simResultsFile = NULL,
  paFile = NULL,
  truncationDistance = max(TL[, 1]),
  snrTruncationThreshold = -Inf,
  NcIsTruncated = FALSE,
  siteCode = "",
  densityResultsFile = NULL,
  parallel = FALSE
)
```

## Arguments

- Nc:

  Total count of true positive detections for call density estimate

- capHistTab:

  Capture history table used to derive false positive rate. Can also be
  used to derive probability of detection as a function of SNR for
  estimating \\p_a\\ if snrDetFun is not specified.

  Must contain columns `detect_table1` and `detect_table2`.
  **`detect_table1` is always read as ground truth** (0 = false
  positive, 1 = true positive) and `detect_table2` as the detector under
  investigation – this is inherited unconditionally from
  [`falseDiscoveryRate()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryRate.md)
  and
  [`capHist2snrInfo()`](https://briansmiller.github.io/callDensity/reference/capHist2snrInfo.md),
  which cde calls internally and which do not expose a way to override
  these column names from here.

  For an observer-ground (OG) analysis, where one trusted observer
  stands in for ground truth, `detect_table1` is that observer's own
  detections, and this holds without any extra work.

  For an adjudicated capture-recapture (CR) analysis, where two
  observers are both fallible and neither is ground truth, this requires
  a deliberate construction step: `detect_table1` must be overwritten
  with the adjudicator's verdict (true call vs false positive) before
  the table is passed to cde, *not* left as one of the two raw observers
  being compared. Getting this wrong does not error – it silently
  computes a false discovery rate for one observer against the other,
  which is not a meaningful quantity and will not match the true rate.
  See
  [`vignette("callDensity_snrThreshold")`](https://briansmiller.github.io/callDensity/articles/callDensity_snrThreshold.md),
  section "Calculate call densities", for a worked example of this
  construction.

- SL:

  List containing distribution of source level parameters. SL must
  contain named elements named: mean, sd, and sampleSize.

- TL:

  Data.frame of transmission losses for Monte Carlo simulation. The
  first column must contain the ranges, the remaining columns contain
  transmission losses (in dB) for each radial transect at that range.

- A:

  Scalar indicating the area for the call density estimate (default=1)

- k:

  Number of sensors (default=1)

- season:

  TimeCode specifying month, season, or year for outputs. Timecodes
  follow the same form as those from the World Ocean Atlas.

- snrDetFun:

  (Optional) Linear-model like structure (GLM,GAM,SCAM,VGLM,etc)
  specifying the SNR-detection function to use. If this is not included,
  then the SNR-detection function will be derived from the capture
  history table.

- NL:

  (Optional and Experimental/Unsupported) list or data.frame containing
  distribution of noise level parameters. Must contain items named:
  mean, sd, and sampleSize (similar to SL). If not provided, NL will be
  estimated from the capture history table (and snrDetFun).

- modelType:

  Either 'glm','gam','scam',or 'vglm' indicating which type of model to
  use for estimating probability of detection as a function of SNR
  (default='scam')

- numKnots:

  Number of knots to use wor modelling detection function vs SNR (only
  if model types 'gam' or 'scam') default = 5

- output.resolution.m:

  Spatial resolution along radial transects for estimating probability
  of detection in area. Default = 100

- outerloop:

  Number of bootstrap iterations for Monte Carlo simulation

- transectFile:

  (Optional) Output text file name where probability of detection will
  be saved. This includes the probability of detection at each
  output.resolution.m step, so potentially yields a large file (hundreds
  of MB). Default is to not write this information to disk.

- simResultsFile:

  (Optional) Name of output text file where Monte Carlo Simulation
  parameters and results will be saved. Default is to not write this
  information to disk.

- paFile:

  - (Optional) Name of output text file where mean probability of
    detection per transect and overall area will be saved. Default is to
    not write this information to disk.

- truncationDistance:

  scalar or matrix of truncation distances. If a matrix is provided,
  then the dimensions should be 1xN with each column corresponding to a
  radial transect and N being the the total number of transects.

- snrTruncationThreshold:

  (Experimental) scalar SNR in dB below which the probability of
  detection is forcibly set to zero. Use this when the detection
  function is not identified at low SNR, which happens when the two
  observers stop contributing overlapping detections (recaptures) there.
  Set it where the recaptures run out, not where the sample runs out.

  Truncation must be applied consistently to every term that counts
  detections. `c` and the detection function are truncated by cde; `Nc`
  cannot be, because cde receives it as a number, so you must truncate
  it yourself and confirm with `NcIsTruncated`. `A` is unchanged, and
  the noise level is estimated from the untruncated sample, because
  noise is a property of the ocean rather than of the threshold.

  The returned density is that of all calls, not of above-threshold
  calls. The fraction of calls arriving above the threshold, q(theta),
  is supplied by SL, TL and NL, and cancels between `Nc` and `p_a`.
  Truncation therefore buys accuracy at the cost of precision: `Nc`
  shrinks.

- NcIsTruncated:

  Logical confirmation that `Nc` counts only detections with SNR \>=
  `snrTruncationThreshold`. Required (and only used) when
  `snrTruncationThreshold` is finite. cde never sees the detections `Nc`
  was counted from, so it cannot verify this itself, and an untruncated
  `Nc` paired with a truncated `p_a` inflates `Dc` silently.

- siteCode:

  (Optional) string containing a code or label associated with the site.
  Default is ” (i.e. blank).

- densityResultsFile:

  (Optional) name of csv file where final call density results will be
  written as a data.frame

- parallel:

  Logical. If TRUE, the Monte Carlo simulation inside pDetInArea runs in
  parallel via future.apply::future_lapply, using whatever future plan
  is currently active. Set the plan before calling cde(), e.g.
  `future::plan(future.callr::callr, workers = 30)`. Default FALSE
  (serial).

  Parallel execution provides significant speedups for computationally
  expensive models (e.g. VGLM/VGAM detection functions from the VGAM
  package), where each Monte Carlo iteration involves multiple calls to
  `VGAM::predict`. For lightweight models (e.g. scam, gam, glm), the
  per-iteration work is too fast relative to the overhead of dispatching
  to and collecting from parallel workers, and `parallel=TRUE` will
  typically be *slower* than serial execution.

  As a rule of thumb: use `parallel=TRUE` for VGLM-based detection
  functions (`modelType='vglm'`) and `parallel=FALSE` (the default) for
  everything else.

  Requires the future.apply package. If not installed, falls back to
  serial execution with a warning.

## Value

`cde` returns a data.frame containing the results of the call density
estimate \\\hat{D}\_c\\, intermediate results such as \\p_a, c\\,
coefficient of variation of these terms, and the input parameters.

## Details

`cde` estimates call density with all parameters included as function
arguments. This is in contrast to the previous operation, i.e. where a
data.frame of parameters was used to indicate which files to load from
disk. Here no parameter file/data.frame is used, so all data and
parameters are must already be objects in memory;
