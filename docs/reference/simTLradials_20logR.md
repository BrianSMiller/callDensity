# Simulate transmission loss (TL) radials following geometric (spherical) spreading law.

Create a TL data.frame with a column for range and subsequent columns
for transmission loss, TL for radial transects around a recorder.
Transects will have a total length of maxRange with increments between
points along the transect of rangeStep. The number of radial transects
is specified as numTransects such that the angular resolution of the
transects will be 360/numTransects.

## Usage

``` r
simTLradials_20logR(maxRange, rangeStep, numTransects)
```

## Arguments

- maxRange:

  - Length of each radial transect

- rangeStep:

  - increment between points on each transect where TL

- numTransects:

  - Number of radial transects (i.e. columns in result)

## Value

TL -
