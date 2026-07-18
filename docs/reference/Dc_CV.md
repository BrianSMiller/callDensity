# Coefficient of Variation (CV) for call density (Dc)

The uncertainty (CV) in the estimated call density (Dc) is calculated
using the Delta Method, which combines the CV of the individual
parameters (Nc, c and pa).

## Usage

``` r
Dc_CV(CV.Nc, CV.pa, CV.c)
```

## Arguments

- CV.Nc:

  - CV of number of calls Nc

- CV.pa:

  - CV of probability of detection in the area, pa

- CV.c:

  - CV of the false discovery rate, c

## Value

CV.Dc - coefficient of variation of the call density

## Details

Here, the uncertainties of Nc, c, and pa are assumed to be independent,
given the fundamental premise of the Delta Method in this application.
