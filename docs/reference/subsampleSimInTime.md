# Subsample from a simulation at evenly spaced time intervals.

Subsample from a simulation at evenly spaced time intervals.

## Usage

``` r
subsampleSimInTime(
  sim,
  minDate = min(sim$datetime),
  maxDate = max(sim$datetime),
  interval = "41 hour",
  duration = 3600
)
```

## Arguments

- sim:

  - simulation data.frame containing column a time column named
    'datetime'

- minDate:

  - posixCT indicating the start of the first subsample

- maxDate:

  - posixCT indicating the start of the last subsample

- interval:

  - difftime interval between subsamples (e.g. '41 hour')

- duration:

  - numeric indicating the duration (in s) of each subsample

## Value

subsampledSim - a simulation dataframe containing only rows from the
input that fall within the subsample
