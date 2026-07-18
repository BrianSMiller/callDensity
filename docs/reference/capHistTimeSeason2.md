# Add R datetime and season to a capture history table (version 2)

Add R datetime and season to a capture history table (version 2)

## Usage

``` r
capHistTimeSeason2(x)
```

## Arguments

- x:

  A capture history table with column t0 containing a matlab datenum

## Value

A capture history table with original columns plus: t: a posixCT
corresponding to the date and time of t0 season: a factor containing
summer, autumn, winter, spring, or year corresponding to the season of
t0 month: a factor containing '01'-'12' corresponding to the month of t0
