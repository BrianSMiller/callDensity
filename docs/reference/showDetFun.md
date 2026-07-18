# Plot SNR detection function with observed SNR distributions

Visualises a fitted SNR detection function together with the
distribution of detected and missed events. By default draws a rug along
the y=1 (detected) and y=0 (missed) reference lines. Optionally overlays
mirrored density curves or histograms above and below those lines.

## Usage

``` r
showDetFun(model, ...)
```

## Arguments

- model:

  Fitted detection model from fitDetFun().

- SNRinfo:

  Data frame used to fit the model. Must contain columns `SNR` and
  `Detected` (logical or 0/1).

- distribution:

  One of `"none"` (default), `"density"`, or `"histogram"`. Controls
  whether a mirrored distribution is drawn above y=1 (detected) and
  below y=0 (missed). A rug is always drawn unless `rug = FALSE`.

- rug:

  Logical. If `TRUE` (default), draws tick marks along the y=1 line for
  detected events and along y=0 for missed events. For large datasets, a
  random subsample of `rug.max` points is used.

- rug.max:

  Maximum number of rug ticks per group. Default 500. Set to `Inf` to
  show all points (slow for n \> 2000).

- mirror.height:

  Maximum height of mirrored distributions as a fraction of the y-axis
  range. Distributions sit outside the \[0,1\] box. Default 0.20.

- rug.alpha:

  Alpha for rug ticks. Default 0.4.

- npoints:

  Number of SNR grid points for the fitted curve. Default 300.

- show.counts:

  Logical. If `TRUE` (default), adds a subtitle reporting the number of
  detections out of total observations.

## Value

A ggplot object.
