# Predict the mean noise level you would measure at detections

Given a candidate true noise distribution, predict the mean of the noise
levels that would be measured at detected calls. Those measurements are
biased low, because a quiet period gives a high SNR and a high SNR gives
a detection, so detections over-represent quiet periods.

## Usage

``` r
predictSampledNL(
  mu,
  sigma,
  detFun,
  SL,
  TL,
  truncationDistance = max(TL[[1]]),
  nNodes = 121,
  ...
)
```

## Arguments

- mu:

  Candidate true mean noise level in dB.

- sigma:

  Standard deviation of the true noise distribution in dB.

- detFun:

  Detection function. See `pDetGivenNL`.

- SL:

  Source level distribution, with elements named mean and sd.

- TL:

  Transmission loss table. See `pDetGivenNL`.

- truncationDistance:

  Scalar or one value per transect, in metres.

- nNodes:

  Number of nodes used to integrate over the noise distribution.

- ...:

  Passed to `pDetGivenNL`, e.g. `binWidth`.

## Value

Scalar. Always less than or equal to `mu`.

## Details

The weighting is `pDetGivenNL`: if the chance of detection at 80 dB is
2.5 times the chance at 84 dB, then 80 dB periods appear among the
detections 2.5 times more often than they occur.
