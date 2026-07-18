# Call density estimate from a parameter file (deprecated)

Deprecated. Use
[`cde`](https://briansmiller.github.io/callDensity/reference/cde.md)
directly. This was a near-complete duplicate of `cde`'s body, reading
its inputs from a bundled parameter object instead of explicit
arguments. Being a separate copy, it had drifted from `cde`: it called
[`pa_CV()`](https://briansmiller.github.io/callDensity/reference/pa_CV.md)
without transect area weights (`cde` always weights by
`truncationDistance^2`), and it had no equivalent of `cde`'s
`NcIsTruncated` guard. This wrapper derives `cde`'s arguments from `p`
and calls `cde` directly, so it now gets exactly the same area-weighted,
truncation-aware behaviour as `cde` itself, rather than a second,
unmaintained implementation of the same calculation.

## Usage

``` r
cdeFromParamFile(
  p,
  season,
  snrDetFun = NULL,
  truncationDistance = Inf,
  snrTruncationThreshold = -Inf,
  NL = NULL
)
```

## Arguments

- p:

  A data.frame containing parameters for call density estimation.
  Usually created by calling function defaultOutputFileNames. Expected
  to contain: detectorParams\$fullYearDetectionCsv,
  detectorParams\$fullYearEffortFile, capHistFile,
  slParams\$slMean/slStd/slSampleSize, tlParams\$tlFile, w, k,
  modelType, numKnots, output.resolution.m, outerloop, transectFile,
  simResultsFile, paFile, siteCode, densityResultsFile.

- season:

  TimeCode specifying month, season, or year for outputs

- snrDetFun:

  OPTIONAL linear-model like structure (GLM,GAM,SCAM,etc) specifying the
  SNR-detection function to use. If this is not included, then the
  SNR-detection function will be derived from the capture history table.

- truncationDistance:

  scalar or matrix of truncation distances. If a matrix is provided,
  then the dimensions should be 1xN with N being the same as the number
  of transects

- snrTruncationThreshold:

  scalar SNR in dB below which the probability of detection will be
  forcibly set to zero.

- NL:

  data.frame containing distribution of noise level parameters. This
  data.frame must contain the rows mean, sd, and sampleSize (similar to
  SL).

## Value

data.frame containing call density inputs, results, and CVs

## Details

\$D_c = \fracN_c\*(1-c)kTP_a\piw^2\$

where: \$D_c\$ is call density \$N_c\$ is number of calls \$c\$ is false
discovery rate \$k\$ is number of sensors (here always 1) \$T\$ is the
duration of data analysed \$P_a\$ is the probability of detection in the
study area \$\piw^2\$ is the study area (in km\\2)
