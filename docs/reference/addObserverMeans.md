# Compute row-wise mean across observer columns

For each entry in `prefixes`, finds matching columns and appends a new
column to `d` containing the row-wise mean, ignoring `NA`s.

## Usage

``` r
addObserverMeans(
  d,
  prefixes = list(SNR = "snr_observer", noiseRMSdB = "noiseRMSdB_observer", signalRMSdB =
    "signalRMSdB_observer")
)
```

## Arguments

- d:

  A data frame containing observer metric columns.

- prefixes:

  A named list or named character vector where names become the new
  column names and values are either a single string pattern or a
  character vector of explicit column names. Defaults to
  `list(SNR = "snr_observer", noiseRMSdB = "noiseRMSdB_observer", signalRMSdB = "signalRMSdB_observer")`.

## Value

The input data frame `d` with one additional column per entry in
`prefixes`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: pattern-based
d <- addObserverMeans(d)

# Mixed: pattern-based and explicit
d <- addObserverMeans(d, prefixes = list(
  SNR        = "snr_observer",
  noiseRMSdB = c("noiseRMSdB_observer1", "noiseRMSdB_observer2")
))
} # }
```
