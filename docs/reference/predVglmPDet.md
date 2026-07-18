# Predict probability of detection from a VGLM SNR-detection function

Helper that combines the per-observer detection probability and the
"any-observer detected" probability from a fitted VGLM into a single
marginal probability of detection at each supplied SNR.

## Usage

``` r
predVglmPDet(...)
```

## Arguments

- snrDetFun.vglm:

  A fitted VGLM SNR-detection function (e.g. from
  [`fitSNRvglm`](https://briansmiller.github.io/callDensity/reference/fitSNRvglm.md))
  carrying `@extra$whichObserver`.

- snrs:

  A data.frame with at least one column named `SNR`, giving the SNR
  values at which to predict the probability of detection.

- whichObserver:

  Character. Name of the observer column to use for predictions.
  Defaults to `snrDetFun.vglm@extra$whichObserver` (the observer baked
  into the fitted model by
  [`fitSNRvglm`](https://briansmiller.github.io/callDensity/reference/fitSNRvglm.md)).

## Value

Numeric vector of probabilities of detection (one per row in `snrs`) for
the chosen observer.
