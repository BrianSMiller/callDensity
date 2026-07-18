# Mean and standard deviation of noise levels by month, season, or year.

Read a table of noise levels containing a columns t and NL for time and
level, respectively.

## Usage

``` r
noiseLevelDistribution(nlFile, season = "year")
```

## Arguments

- nlFile:

  File name of csv containing noise level times and values

- season:

  TimeCode for season, month, or 'year'

## Value

The mean and standard deviation from the table of noise levels.
