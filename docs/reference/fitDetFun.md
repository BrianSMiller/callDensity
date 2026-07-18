# Fit an SNR detection function

Fits a detection function relating probability of detection to SNR.

Supported model types are:

- glm

- gam

- scam

- vglm (capture-recapture)

Fits a detection function relating probability of detection to SNR.

Supported model types are:

- glm

- gam

- scam

- vglm (capture-recapture)

## Usage

``` r
fitDetFun(
  SNRinfo,
  modelType = c("gam", "glm", "scam", "vglm"),
  numKnots = 3,
  yColNames = c("detect_observer1", "detect_observer2"),
  whichObserver = NULL
)
```

## Arguments

- SNRinfo:

  data.frame containing SNR and detection information.

- modelType:

  character.

- numKnots:

  spline basis dimension for GAM/SCAM.

- yColNames:

  observer columns for capture-recapture models.

- whichObserver:

  observer used for prediction.
