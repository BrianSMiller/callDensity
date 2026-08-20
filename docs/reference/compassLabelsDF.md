# Compass-label positions for a radial plot

Compass-label positions for a radial plot

## Usage

``` r
compassLabelsDF(rMax, inset = 0.08)
```

## Arguments

- rMax:

  Outer radius (km) of the plotted disc.

- inset:

  Fractional distance beyond `rMax` at which to place the labels.
  Default `0.08` (8% beyond the edge).

## Value

A data.frame with columns `x`, `y`, `label`.
