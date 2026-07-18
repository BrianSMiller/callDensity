# Calculate deployment duration from a Matlab soundFolder csv file.

Calculate deployment duration from a Matlab soundFolder csv file.

## Usage

``` r
deploymentDurationFromsoundFolderCsv(fullYearEffortFile, season = "year")
```

## Arguments

- fullYearEffortFile:

  - csv file containing a row for each time period that was included in
    the study and columns named startDate and duration. startDate should
    be a matlab datenum, and duration should be numeric and in seconds.

- season:

  - default='year'
