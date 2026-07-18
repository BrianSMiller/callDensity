# Confidence interval from density and it's CV

Quick function to convert an abundance with a coefficient of variation
(CV) into a 95% confidence interval.

## Usage

``` r
ciFromCV(a, cv)
```

## Arguments

- a:

  abundance (or density)

- cv:

  coefficient of variation, in decimal format, i.e., a CV of 10% would
  be 0.1.

## Value

list with lower and upper 95% confidence intervals
