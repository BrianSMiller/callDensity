# Estimate the noise level distribution from noise measured at detections

Noise levels measured at detections are biased low, because detections
over-represent quiet periods. This undoes that bias.

## Usage

``` r
nlFromDetections(
  snrInfo,
  snrDetFun,
  SL,
  TL,
  truncationDistance = max(TL[[1]]),
  nlColumn = "NoiseRL",
  searchWidth = 25,
  ...
)
```

## Arguments

- snrInfo:

  Table of SNR information containing a column of noise level
  measurements in dB.

- snrDetFun:

  Detection function, as passed to `pDetInArea`.

- SL:

  Source level distribution, with elements named mean and sd.

- TL:

  Transmission loss table. First column ranges in metres, remaining
  columns TL in dB per radial transect.

- truncationDistance:

  Scalar or one value per transect, in metres.

- nlColumn:

  Name of the noise level column. Default "NoiseRL".

- searchWidth:

  Width in dB of the interval searched above the measured mean. The bias
  cannot be negative, so the search runs upwards only.

- ...:

  Passed to `pDetGivenNL`, e.g. `binWidth`.

## Value

Data.frame with one row and columns mean, sd and sampleSize, the same
format as `nlFromSnrInfo` and `noiseLevelDistribution`.

## Details

Replaces `nlFromSnrInfo`, which corrected the same bias by adding the
SNR at which the detection function reaches 0.5. That quantity is a
property of the detector. The bias is a property of the propagation
geometry and the noise variance. The two coincide only by chance. See
the noiseLevels vignette.

The standard deviation is taken from the measurements directly.
Filtering by detection shifts the mean but barely narrows the
distribution, so the measured standard deviation is close to the truth
even though the measured mean is not. This leaves one unknown, found by
[`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html).
