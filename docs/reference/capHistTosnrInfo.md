# Convert capture history DATA.FRAME into the 'SNRinfo' format used by the callDensity package

**Deprecated.** A thin wrapper around
[`chtToSNRinfo`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md),
kept so the published Common Ground and Beyond Counting Calls analyses
keep reproducing exactly against current `main`. New code should call
`chtToSNRinfo` directly.

## Usage

``` r
capHistTosnrInfo(capHistTab)
```

## Arguments

- capHistTab:

  - Capture history table of detections. Must have
    `detect_table1`/`detect_table2`, a `signalRMSdB1`/ `signalRMSdB2`
    and `noiseRMSdB1`/`noiseRMSdB2` pair (averaged, matching original
    behaviour), and `datetime`.
