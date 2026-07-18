# Plot SNR histogram by observer

Pivots SNR observer columns to long format and produces a faceted
histogram of SNR distributions, one panel per observer.

## Usage

``` r
plotSNRHistogram(
  df,
  obs_labels,
  snr_prefix = "snr_observer",
  snr_lims = NULL,
  binwidth = 0.5,
  title = "SNR distribution by detector",
  ...
)
```

## Arguments

- df:

  A data frame containing SNR observer columns.

- obs_labels:

  Named character vector mapping numeric observer indices (as character)
  to display labels.

- snr_prefix:

  Either a single character string pattern identifying SNR columns, or a
  character vector of explicit column names. Defaults to
  `"snr_observer"`.

- snr_lims:

  Numeric vector of length 2 giving x-axis limits. If `NULL` (default),
  limits are computed from the data using
  [`computeSNRLims`](https://briansmiller.github.io/callDensity/reference/computeSNRLims.md).

- binwidth:

  Histogram bin width in dB. Defaults to `0.5`.

- title:

  Plot title. Defaults to `"SNR distribution by detector"`.

- ...:

  Additional ggplot2 layers (e.g. scales, themes, annotations) passed to
  the plot via `+`.

## Value

A `ggplot` object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
plotSNRHistogram(ap, obs_labels,
  title = "Adjudicated subset: SNR distribution by detector")

# Shared limits with extra layer
lims <- computeSNRLims(d)
plotSNRHistogram(ap, obs_labels, snr_lims = lims,
  title = "Adjudicated subset: SNR distribution by detector",
  ggplot2::scale_fill_brewer(palette = "Dark2"))
plotSNRHistogram(d, obs_labels, snr_lims = lims,
  title = "Full dataset: SNR distribution by detector",
  ggplot2::scale_fill_brewer(palette = "Dark2"))
} # }
```
