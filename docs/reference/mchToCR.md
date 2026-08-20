# Convert a multi-observer capture history table to a two-observer CR table

Extracts a two-observer (table1, table2) capture-recapture capture
history table from a multi-observer source table by selecting only rows
where at least one of the two named observers detected the event, and
renaming the relevant columns to the `_table1`/`_table2` suffixes
expected by the rest of the callDensity package. The `verdict` column
from the source is used as the ground-truth detection for table1.

## Usage

``` r
mchToCR(d, table1suffix, table2suffix)
```

## Arguments

- d:

  A multi-observer capture history data.frame containing per-event
  columns `t0`, `tEnd`, `fLow`, `fHigh`, `SNR`, `signalRMSdB`,
  `noiseRMSdB`, `noiseDev`, `verdict`, plus per-observer detection
  columns named `detect_<observerSuffix>`.

- table1suffix:

  Character. Suffix of the observer columns to remap to the `_table1`
  role (treated as ground truth via `verdict`).

- table2suffix:

  Character. Suffix of the observer columns to remap to the `_table2`
  role (the detector being evaluated).

## Value

A capture history data.frame in the two-observer (table1, table2) format
expected by
[`cde`](https://briansmiller.github.io/callDensity/reference/cde.md),
with time/season columns added by
[`capHistTimeSeason`](https://briansmiller.github.io/callDensity/reference/capHistTimeSeason.md).

## Deprecated

No current analysis needs this reduction.
[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)'s
only use of `capHistTab` is
[`falseDiscoveryRate()`](https://briansmiller.github.io/callDensity/reference/falseDiscoveryRate.md),
which already accepts `gtColName`/`testColName` directly against native
`detect_observerN` columns (once
[`cde()`](https://briansmiller.github.io/callDensity/reference/cde.md)
forwards them – see its own argument list).
`fitDetFun(modelType = 'vglm')` accepts any number of native-named
observer columns directly via `yColNames`, with no reduction to two
required at all –
[`VGAM::posbernoulli.t`](https://rdrr.io/pkg/VGAM/man/posbernoulli.t.html)
is not limited to two occasions. This function's body is unchanged and
will keep working, kept only so the published Common Ground and Beyond
Counting Calls analysis scripts keep reproducing exactly. New work
should skip the 2-observer reduction entirely and use
[`chtToSNRinfo`](https://briansmiller.github.io/callDensity/reference/chtToSNRinfo.md)
with a `observers` vector of whatever length is needed.
