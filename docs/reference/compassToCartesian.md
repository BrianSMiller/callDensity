# Convert (range, compass azimuth) to Cartesian (x, y)

Convert (range, compass azimuth) to Cartesian (x, y)

## Usage

``` r
compassToCartesian(range_km, azimuth_deg)
```

## Arguments

- range_km:

  Numeric vector of ranges (km).

- azimuth_deg:

  Numeric vector of compass bearings in degrees, 0 = true north,
  increasing clockwise.

## Value

A list with elements `x` and `y`, in the same units as `range_km`.
