# Probability of detection in the study area at a fixed noise level

Answers the question that `pDetInArea` answers, but with the noise level
held still rather than drawn from a distribution: if the noise were
always exactly `nl`, what fraction of the calls in the study area would
be detected?

## Usage

``` r
pDetGivenNL(
  nl,
  detFun,
  SL,
  TL,
  truncationDistance = max(TL[[1]]),
  nSLnodes = 41,
  binWidth = 0.25
)
```

## Arguments

- nl:

  Vector of noise levels in dB.

- detFun:

  Either a detFun object from `fitDetFun`, or a plain function of SNR
  returning probability of detection.

- SL:

  List or data.frame containing the distribution of source levels, with
  elements named mean and sd.

- TL:

  Data.frame of transmission losses. The first column contains ranges in
  metres, the remaining columns contain TL in dB for each radial
  transect at that range. Same format as `cde` expects.

- truncationDistance:

  Scalar or vector of truncation distances in metres. If a vector, one
  value per transect. Cells beyond the truncation distance carry no
  weight, and the returned probability is relative to the truncated
  area.

- nSLnodes:

  Number of quadrature nodes used to average over the source level
  distribution. Default 41.

- binWidth:

  Width in dB of the transmission loss bins. The only approximation in
  this function. Default 0.25, which is conservative: results are
  typically stable to four significant figures at 1 dB.

## Value

Numeric vector of the same length as `nl`, each between 0 and 1.

## Details

Evaluated over a grid of noise levels this gives a curve, and that curve
is what makes the noise levels measured at detections biased. Quiet
periods are over-represented among detections by exactly the ratio the
curve describes. See `nlFromDetections` for the use, and the noiseLevels
vignette for the derivation.

Transmission loss enters only as a lookup table, so this makes no
assumption about the propagation model. Under spherical spreading the
curve falls by one decade per 10 dB, but that is a consequence rather
than an input.
