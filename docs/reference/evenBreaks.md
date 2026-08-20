# Evenly spaced histogram breaks spanning exactly xlim

Used instead of [`pretty()`](https://rdrr.io/r/base/pretty.html) so bins
never extend past the shared panel x-limits –
[`pretty()`](https://rdrr.io/r/base/pretty.html)'s round numbers
routinely do, which silently drops the outermost bin (and warns) once
`scale_x_continuous(limits=)` clips it.

## Usage

``` r
evenBreaks(xlim, n = 25)
```

## Arguments

- xlim:

  Length-2 numeric vector, the shared panel x-range.

- n:

  Number of bins. Default 25.

## Value

A numeric vector of `n + 1` break points.
