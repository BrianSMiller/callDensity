# WARNING: This function is not supported, and probably does not do what whatever you were hoping it might do. Perhaps consider function c_CV instead. CV of false discovery rate using Cochran approximation

WARNING: This function is not supported, and probably does not do what
whatever you were hoping it might do. Perhaps consider function c_CV
instead. CV of false discovery rate using Cochran approximation

## Usage

``` r
falseDiscoveryCV(hourlyFalsePosFile, season)
```

## Arguments

- hourlyFalsePosFile:

  Name of csv file contianing hourly estimates of false positive rate.

- season:

  string or number corresponding to the time of year for which call
  densities should be estimated (and data subsetted).

## Value

CV of false discovery rate using Cochran approximation

## See also

c_CV
