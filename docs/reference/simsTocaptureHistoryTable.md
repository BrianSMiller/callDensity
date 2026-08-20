# Merge two simulated detectors into a capture history table

Merge two simulated detectors into a capture history table

## Usage

``` r
simsTocaptureHistoryTable(
  subsampleDet1,
  subsampleDet2,
  observerSuffix = c("observer1", "observer2")
)
```

## Arguments

- subsampleDet1, subsampleDet2:

  Output of
  [`simulateDetector`](https://briansmiller.github.io/callDensity/reference/simulateDetector.md)
  for each of the two simulated detectors.

- observerSuffix:

  Character vector of length 2 giving the suffix used for the two
  detection-flag columns in the output. Default
  `c('observer1', 'observer2')` – matchbox-native, matching
  [`chtToSNRinfo`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md)'s
  own default `observers` naming, so a simulated capture history table
  needs no renaming before use. Pass `c('table1', 'table2')` for the
  older `detect_table1`/`detect_table2` shape.

  Only the two detection-flag columns are affected. The per-observer
  `groundTruth`/`snr`/`signalRMSdB`/`noiseRMSdB` columns keep their
  existing `1`/`2` numbering regardless – they're simulation-only
  ground-truth/diagnostic columns with no equivalent in real matchbox
  output (which carries one shared `signalRMSdB`/`noiseRMSdB` per
  matched event, already consolidated below, and no ground truth at
  all), so there's no matchbox convention for `observerSuffix` to align
  them to.
