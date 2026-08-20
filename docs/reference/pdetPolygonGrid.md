# Build annular-sector polygon patches for each (range, azimuth) cell

The exact ggplot2 equivalent of MATLAB's
`pcolor(x = r*cos(theta), y = r*sin(theta), z)`: every sampled cell
becomes a quadrilateral (or, with `arcPoints > 2`, a many-sided patch)
at its true Cartesian position. No values are interpolated.

## Usage

``` r
pdetPolygonGrid(long, rangeStep, angleStep, arcPoints = NULL)
```

## Arguments

- long:

  Data.frame from `pdetToLong`.

- rangeStep:

  Range bin width (km).

- angleStep:

  Azimuth bin width (degrees).

- arcPoints:

  Number of vertices used for the inner/outer arc of each cell. `NULL`
  auto-selects: dense (5 degree steps) when there are 8 or fewer
  azimuths, so coarse wedges still render with a curved outer edge;
  straight chords (`arcPoints = 2`) once azimuths are already fine
  enough to look round regardless.

## Value

A data.frame with columns `x`, `y`, `cellID`, `pDet`, one row per
polygon vertex.
