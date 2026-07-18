# Lookup the season for a given POSIXct, x.

Seasons are defined by 3 month periods from Summer=(Dec,Jan,Feb);
Autumn=(Mar,Apr,May), Winter=(Jun,Jul,Aug), Spring=(Sep,Oct,Nov)

## Usage

``` r
time2season(x)
```

## Arguments

- x:

  POSIXct

## Value

Factor containing seasons summer, autumn, winter, or spring
