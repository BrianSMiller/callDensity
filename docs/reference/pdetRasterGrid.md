# Resample (range, azimuth) cells onto a fine square Cartesian grid

Nearest-neighbour resampling, not interpolation: every output pixel
takes the value of its nearest (range, azimuth) sample. Assumes a
constant range step and constant angle step (true of `pDetInArea`'s
output), which is what makes the nearest-cell lookup a closed-form index
calculation rather than a search.

## Usage

``` r
pdetRasterGrid(long, rangeStep, angleStep, rasterRes = 400)
```

## Arguments

- long:

  Data.frame from `pdetToLong`.

- rangeStep:

  Range bin width (km).

- angleStep:

  Azimuth bin width (degrees).

- rasterRes:

  Number of pixels along each axis of the square output grid.

## Value

A data.frame with columns `x`, `y`, `pDet`, one row per pixel (pixels
outside the sampled disc have `pDet = NA`).
