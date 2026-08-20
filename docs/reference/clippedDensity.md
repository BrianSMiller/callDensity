# Kernel density estimate pre-clipped to a shared x-range

Passing `from`/`to` directly to
[`density()`](https://rdrr.io/r/stats/density.html) keeps the curve
within the shared panel x-limits from the start, rather than computing
the usual wider KDE (which by default extends `cut` bandwidths past the
data) and letting ggplot2 silently drop and warn about the clipped tails
once `scale_x_continuous(limits=)` is applied.

## Usage

``` r
clippedDensity(x, xlim)
```

## Arguments

- x:

  Numeric vector.

- xlim:

  Length-2 numeric vector, the shared panel x-range.

## Value

A data.frame with columns `SNR` and `y`.
