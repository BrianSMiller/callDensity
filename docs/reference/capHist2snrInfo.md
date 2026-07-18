# Convert a capture history table into an SNRinfo data.frame

Convert a capture history table into an SNRinfo data.frame

## Usage

``` r
capHist2snrInfo(snr, season = "year")
```

## Arguments

- snr:

  A capture history data.frame containing detections and SNR info. Must
  have columns `detect_table1`, `detect_table2`, `signalRMSdB`,
  `noiseRMSdB`, `t`, `month`, and `season`. Passing a file path here is
  deprecated and will error.

  `detect_table1` is always read as ground truth: rows are filtered to
  `detect_table1 == 1` before anything else happens. For an
  observer-ground (OG) analysis this is the trusted observer's own
  column. For an adjudicated capture-recapture (CR) analysis, where
  neither raw observer is ground truth, `detect_table1` must be
  overwritten with the adjudicator's verdict before calling this
  function – see
  [`cde`](https://briansmiller.github.io/callDensity/reference/cde.md)'s
  `capHistTab` documentation for the same requirement and a pointer to a
  worked example.

- season:

  A timeCode corresponding to months, seasons, or 'year' (the default,
  which returns all rows).

## Value

SNRinfo data.frame with columns `Detected`, `CallRL`, `NoiseRL`, `SNR`,
`t`, `month`, and `season`, filtered to the requested timeCode.
