# Read a capture history csv file (e.g. created in Matlab)

Required columns for each observer are detect\_, snr\_, signalRMSdB\_,
noiseRMSdB\_, and the underscore is followed by input parameter
observerNames. -detect\_ columns are binary i.e. 1 a detection and 0 for
absence -snr\_, signalRMSdB\_, and noiseRMSdB\_ are all numeric and in
dB -groundTruth is the binary column to use as ground-truth detections.
It can be a detect_observerName or another column. OG capture history
tables have exactly two observers with suffixes,
observerNames=c('table1','table2'). These are the type of tables used in
the 'Beyond Counting Calls' manuscript. CR capture history tables can
have N observers and can have any number of suffixes of any form.
MultiCaptureHistoryTable from the common ground manuscript produces
suffixes of the form: ('observer1','observerN-1','observerN')

## Usage

``` r
readCapHist(
  capHistFile,
  observerNames = c("table1", "table2"),
  groundTruth = "table1",
  whichObserver = "detect_table2"
)
```

## Arguments

- capHistFile:

  - name of the csv or txt file that contains capture histories. TODO:
    Document the required columns/file format better

- observerNames:

  - list of strings that correspond to the suffixes used for each
    observer. Default: c('table1','table2')

- groundTruth:

  - Name of column with binary data used as ground-truth.

- whichObserver:

  - Name of observer for which detection function will be modelled, and
    call density estimated

## Value

capture history table as data.frame
