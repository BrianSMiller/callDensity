# Title

Title

## Usage

``` r
countDetections(
  detectionFile,
  season = "year",
  snrTruncationThreshold = -Inf,
  snrColName = "snr"
)
```

## Arguments

- detectionFile:

  - Text file or data.frame containing detections for the full dataset.
    File format must be tab separated with header, and must contain the
    following columns: t0: Matlab datenum snr: Signal to noise ratio of
    that detection (required only if snrTruncationThreshold is other
    than default value of -Inf)

- season:

  - Time period over which to subset detections. Can follow World Ocean
    Atlas numeric time codes (0-16), or be name of season
    ('summer','autumn','winter','spring'), month name/abbreviation, or
    'year' (default).

- snrTruncationThreshold:

  detections with SNR less than snrTruncation threshold (in dB) will be
  excluded from count

- snrColName:

  - name of column containing SNR values (only used if SNR
    truncation\>-Inf)

## Value

number of detections in the detectionFile within specified season and
\>= snrTruncationThreshold
