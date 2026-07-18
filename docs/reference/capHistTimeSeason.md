# Add R datetime, season, and month codes to a capture history table

Add R datetime, season, and month codes to a capture history table

## Usage

``` r
capHistTimeSeason(x)
```

## Arguments

- x:

  Capture history table with column t0_tableN of Matlab datenum's

## Value

Capture history table with column t, season, and month containing
POSIXct, and timeCode factors for season and month.
