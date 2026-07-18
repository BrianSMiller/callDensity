# Default names of output files for call density estimation

Given data.frame of call density parameters, p, create default
parameters and output file-names for Monte-Carlo modelling and
estimating call-density. Append all of these parameters to the
data.frame and return it.

## Usage

``` r
defaultOutputFileNames(p, season, outputFolder = ".")
```

## Arguments

- p:

  Data.frame containing parameters for call density estimation

- season:

  TimeCode to filter month, season, or full year

- outputFolder:

  Location/path/folder where output files will be saved default is '.'
  (the current working directory).

## Value

Parameter file with file names corresponding to the correct timeCode
