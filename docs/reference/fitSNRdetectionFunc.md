# Fit an SNR-detection function

SNR detection functions are binomial models of the form detected ~ SNR.
These can be a GAM (default), GLM, or SCAM.

## Usage

``` r
fitSNRdetectionFunc(SNRinfo, modelType = "gam", numKnots = 3)
```

## Arguments

- SNRinfo:

  A data.frame containing detection information in columns Detected and
  SNR

- modelType:

  one of either 'GAM', 'GLM', or 'SCAM'. ModelType determines which type
  of model will be used to estimate the detection-SNR function.

- numKnots:

  Number of knots in the GAM or SCAM
