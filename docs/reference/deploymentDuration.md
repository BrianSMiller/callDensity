# Calculate deployment duration from recording start times and durations

Calculate deployment duration from recording start times and durations

## Usage

``` r
deploymentDuration(duration, startDate, season = "year")
```

## Arguments

- duration:

  - list of recording durations (in s)

- startDate:

  - list of start dates for each duration

- season:

  - Subset the startDate by this season/month/year (default='year')

## Value

- duration of effort (as numeric)
