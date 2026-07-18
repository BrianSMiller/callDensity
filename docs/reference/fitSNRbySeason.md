# Model SNR-detection function with season as a factor (deprecated?)

Model SNR-detection function with season as a factor (deprecated?)

## Usage

``` r
fitSNRbySeason(SNRinfo, season = year, useGLM = TRUE, numKnots = 3)
```

## Arguments

- SNRinfo:

  Data.frame containing columns SNR, Detected, and season

- season:

  Factor containing a timeCode

- useGLM:

  Boolean, if true fit a GLM instead of GAM

- numKnots:

  Number of knots to use in the GAM
