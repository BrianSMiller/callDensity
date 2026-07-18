# Simulate animals calls using a uniform distribution time and space

Generate a uniform distribution of locations in time and space. The
space component is in two dimensions (x, y), and the time component is
in calendar time (posix.ct).

## Usage

``` r
simCallLocation(
  n = 1e+06,
  R = 1e+06,
  minDate = Sys.time(),
  maxDate = minDate + 86400
)
```

## Arguments

- n:

  - Number of calls that will be created.

    For reference 1 whale call every 10 s would be a maximum of 3.1536M
    calls In a continuous year of recording (i.e. 365*24*60\*60/10 =
    3.1536e6)

- R:

  - Radius of study area.

- minDate:

  - Starting date and time for the simulated data (POSIXct)

- maxDate:

  - Ending/latest date and time for the simulated data (POSIXct)

## Value

sim - data.frame containing a row for each simulated call and columns
for the location and time of each call.

## Details

simCallLocation creates n calls within radius R of the receiver
(hydrophone). This is useful for testing callDensity package with a
known call density (uniform random in circular study area).
