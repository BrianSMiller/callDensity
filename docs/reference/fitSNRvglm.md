# Fit an SNR-detection function with a closed population capture-recapture glm

VGLM SNR-detection functions are positive bernoulli models of the form
(y1,y2,...yN) ~ SNR. Here yi is a column of detections from the ith
observer (with a 1 for detected and 0 for not detected by that
observer). The assumptions for these models is that detection
probability depends on heterogeneity from SNR, observers are
independent, and there are no false positive detections included in the
data for either observer. Models are fitted using the VGAM package (Yee,
Stoklosa & Huggins 2015).

Yee, Thomas W., Jakub Stoklosa, and Richard M. Huggins. “The VGAM
Package for Capture-Recapture Data Using the Conditional Likelihood.”
Journal of Statistical Software 65 (June 1, 2015): 1–33.
https://doi.org/10.18637/jss.v065.i05.

## Usage

``` r
fitSNRvglm(
  SNRinfo,
  yColNames = c("detect_observer1", "detect_observer2"),
  whichObserver = "detect_observer2"
)
```

## Arguments

- SNRinfo:

  A data.frame containing detection information in columns Detected and
  SNR

- yColNames:

  a list of strings containing the column names for each set of observer
  detections (default: 'detect_observer1','detect_observer2')

- whichObserver:

  column name of the observer to use for predictions
