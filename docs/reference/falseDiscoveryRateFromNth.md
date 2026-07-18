# WARNING: This function is not supported, and probably does not do what whatever you were hoping it might do. Use functions falseDiscoverRate and c_CV instead. False discovery rate from inspection of every Nth detection

WARNING: This function is not supported, and probably does not do what
whatever you were hoping it might do. Use functions falseDiscoverRate
and c_CV instead. False discovery rate from inspection of every Nth
detection

## Usage

``` r
falseDiscoveryRateFromNth(falsePositiveXlsx, season = "year")
```

## Arguments

- falsePositiveXlsx:

  Excel spreadsheet with timestamp of false positives. This spreadsheet
  requires two columns: one called UTC and the other called 'True
  Positive Rate.' These should be in a worksheet called 'conference'.
  The file used by this function must already have had SNR truncation
  applied (i.e. this function has no means of appling an
  snrTruncationThreshold)

- season:

  Month or season over which to subset the data. Months can be 01-12,
  and seasons can be 'summer','autumn','winter','spring', or 'year'.

## See also

falseDiscoveryRate, c_CV
