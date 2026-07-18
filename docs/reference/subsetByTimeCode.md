# Subset a data.frame by by months or season

Subset a data.frame by timeCode, where timeCode can be either a 'season'
(summer, autumn, winter, spring), a month characters ('01'-'12'), or
'year', with 'year' the same as including all data

## Usage

``` r
subsetByTimeCode(df, dt, timeCode)
```

## Arguments

- df:

  A data.frame

- dt:

  A posixCT datetime

- timeCode:

  A timecode indicating the season, month, or 'year' for all data
