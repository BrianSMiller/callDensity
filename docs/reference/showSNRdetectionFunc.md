# Plot an SNR detection function with mirrored SNR distributions

Shows the fitted detection function together with the distribution of
detected and missed true positives. Detected calls are plotted above the
detection curve and missed calls below.

## Usage

``` r
showSNRdetectionFunc(...)
```

## Arguments

- model:

  A fitted detection model returned by fitSNRdetectionFunc(),
  fitSNRvglm(), or fitSNRvgam().

- SNRinfo:

  Original data.frame used to fit the model.

- whichObserver:

  Observer column for vglm/vgam models. Defaults to the observer stored
  in model@extra\$whichObserver.

- distribution:

  One of "density", "histogram", or "none".

- mirror.height:

  Height of the mirrored distributions.

- npoints:

  Number of SNR values for prediction.

## Value

A ggplot object.
