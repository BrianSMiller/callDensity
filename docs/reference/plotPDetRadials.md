# Plot probability of detection as a directional range footprint

Renders p(det) as a function of range and azimuth around the recorder,
in Cartesian coordinates – the direct ggplot2 equivalent of Brian's
MATLAB `plotPDetRadials.m`, which builds the same picture via
`pcolor(x = r*cos(theta), y = r*sin(theta), z)`. Uses the same default
ggplot2 colour scale and theme as
[`plotSpatialDetections`](https://briansmiller.github.io/callDensity/reference/plotSpatialDetections.md).

## Usage

``` r
plotPDetRadials(
  pDetResults,
  maxRange = NULL,
  method = c("polygon", "raster"),
  arcPoints = NULL,
  rasterRes = 400
)
```

## Arguments

- pDetResults:

  Either the full list returned by
  [`pDetInArea`](https://briansmiller.github.io/callDensity/reference/pDetInArea.md)
  (in which case `$allDetFunctions` is used), or that data.frame
  directly – `range_m` in the first column, one column per transect
  named `tl<angle>` (matching
  [`simTLradials_20logR`](https://briansmiller.github.io/callDensity/reference/simTLradials_20logR.md)'s
  convention) giving p(det) at every range along that transect. `angle`
  is a compass bearing, 0 = true north, increasing clockwise.

- maxRange:

  Optional upper limit (m) for the range axis. Default `NULL` uses the
  full range present in the data.

- method:

  `"polygon"` (default) draws each cell as an exact annular-sector patch
  – no interpolation, the direct equivalent of MATLAB's `pcolor`.
  `"raster"` nearest-neighbour-resamples onto a fine square grid
  instead, which can look smoother at very coarse azimuthal resolution
  but is drawing resampled, not sampled, pixels.

- arcPoints:

  Only used when `method = "polygon"`. Number of vertices used for each
  cell's inner/outer arc. `NULL` (default) auto-selects based on the
  number of azimuths present – see
  [`pdetPolygonGrid`](https://briansmiller.github.io/callDensity/reference/pdetPolygonGrid.md).
  Set explicitly to override, e.g. to force straight chords
  (`arcPoints = 2`) even at coarse azimuthal resolution.

- rasterRes:

  Only used when `method = "raster"`. Number of pixels along each axis
  of the output grid. Default `400`.

## Value

A ggplot object.

## Details

Assumes a constant range step and constant azimuth step (true of
`simTLradials_20logR`'s output and of `pDetInArea`'s own
`output.resolution.m` spacing).
