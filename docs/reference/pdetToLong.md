# Pivot a pDetInArea result into long (range, azimuth, pDet) form

Pivot a pDetInArea result into long (range, azimuth, pDet) form

## Usage

``` r
pdetToLong(pDetResults, maxRange = NULL)
```

## Arguments

- pDetResults:

  Either the full list returned by `pDetInArea`, or the
  `allDetFunctions` data.frame directly.

- maxRange:

  Optional upper limit (m) for the range axis.

## Value

A data.frame with columns `range_km`, `azimuth`, `pDet`.
