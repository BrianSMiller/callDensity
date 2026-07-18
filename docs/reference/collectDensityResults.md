# Load specified density results text files in a directory

Load specified density results text files in a directory

## Usage

``` r
collectDensityResults(path, densityString = paste0("^density_", siteCode))
```

## Arguments

- path:

  Location (path) of text files containing call density results

- densityString:

  Prefix of \_density files to be loaded

## Value

Data.frame of all call density results in the specified path
