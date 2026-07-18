# Compute pretty SNR limits across observer columns

Finds the minimum and maximum SNR values across all matching observer
columns and returns pretty axis limits.

## Usage

``` r
computeSNRLims(df, snr_prefix = "snr_observer")
```

## Arguments

- df:

  A data frame containing SNR observer columns.

- snr_prefix:

  Either a single character string pattern identifying SNR columns, or a
  character vector of explicit column names. Defaults to
  `"snr_observer"`.

## Value

A numeric vector of length 2 giving `c(min, max)` pretty limits.

## Examples

``` r
if (FALSE) { # \dontrun{
computeSNRLims(ap)
computeSNRLims(d, snr_prefix = c("snr_observer1", "snr_observer2"))
} # }
```
