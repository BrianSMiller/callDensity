# Pivot SNR observer columns to long format

Helper used by the SNR plotting functions. Selects columns matching
`snr_prefix` and pivots them to long format, attaching observer labels.

## Usage

``` r
pivotSNR(df, snr_prefix = "snr_observer", obs_labels, extra_cols = NULL)
```

## Arguments

- df:

  A data frame containing SNR observer columns.

- snr_prefix:

  Either a single character string pattern identifying SNR columns (e.g.
  `"snr_observer"`), or a character vector of explicit column names.
  Defaults to `"snr_observer"`.

- obs_labels:

  Named character vector mapping numeric observer indices (as character)
  to display labels, e.g. `c("1"="Alice","2"="Bob")`.

- extra_cols:

  Integer or character vector of additional column indices or names to
  retain (e.g. `"i"` for a time index). Defaults to `NULL`.

## Value

A long-format data frame with columns `observer` (factor) and `snr`,
plus any columns specified in `extra_cols`.
