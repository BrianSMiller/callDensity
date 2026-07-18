# Plot SNR time series for the adjudicated subset

Pivots SNR observer columns to long format and produces a scatter plot
of SNR against a time index column, coloured by observer.

## Usage

``` r
plotSNRTimeSeries(
  ap,
  obs_labels,
  snr_prefix = "snr_observer",
  time_col = "i",
  snr_lims = NULL,
  title = "Adjudicated subset: SNR time series by detector",
  ...
)
```

## Arguments

- ap:

  A data frame containing SNR observer columns and a time index column.

- obs_labels:

  Named character vector mapping numeric observer indices (as character)
  to display labels.

- snr_prefix:

  Either a single character string pattern identifying SNR columns, or a
  character vector of explicit column names. Defaults to
  `"snr_observer"`.

- time_col:

  Name of the time index column. Defaults to `"i"`.

- snr_lims:

  Numeric vector of length 2 giving y-axis limits. If `NULL` (default),
  limits are computed from the data using
  [`computeSNRLims`](https://briansmiller.github.io/callDensity/reference/computeSNRLims.md).

- title:

  Plot title. Defaults to
  `"Adjudicated subset: SNR time series by detector"`.

- ...:

  Additional ggplot2 layers (e.g. scales, themes, annotations) passed to
  the plot via `+`.

## Value

A `ggplot` object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
plotSNRTimeSeries(ap, obs_labels)

# Shared limits with extra layer
lims <- computeSNRLims(d)
plotSNRTimeSeries(ap, obs_labels, snr_lims = lims,
  ggplot2::theme(legend.position = "bottom"))
} # }
```
