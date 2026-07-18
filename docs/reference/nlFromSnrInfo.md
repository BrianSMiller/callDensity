# Estimate noise level distribution for MC simulation from the noise measurements included in the snrInfo file, plus the SNR detection function.

Estimate noise level distribution for MC simulation from the noise
measurements included in the snrInfo file, plus the SNR detection
function.

## Usage

``` r
nlFromSnrInfo(snrInfo, snrDetFun)
```

## Arguments

- snrInfo:

  - Table of SNR information containing column named NoiseRL with noise
    level measurements in dB

- snrDetFun:

  - SNR detection function (e.g. from fitSNRdetectionFunc or fitSNRvglm)

## Value

Data.frame with 1 row and 3 columns containing parameterised
distribution of noise levels. Column names are mean, sd, and sampleSize.
Distribution assumed to be normal.
