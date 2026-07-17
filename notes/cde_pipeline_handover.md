# CDE Pipeline — Format Inventory and Handover Notes
<!-- Prepared: 2026-07-03 -->
<!-- Context: callDensity v1.0.0, bsnr v0.3.1-beta, matchbox (MATLAB), CopernicusPyRAM v0.1.0-alpha -->
<!-- Purpose: Handover to new session for holistic pipeline design discussion -->

---

## Background

The call density estimation (CDE) pipeline has accumulated technical debt
from iterative development across multiple projects. The core tools are now
stable enough to reason about the interfaces between them cleanly:

- **matchbox** (MATLAB) — acoustic event matching and adjudication → CHT
- **bsnr** (MATLAB) — SNR estimation from audio
- **CopernicusPyRAM** (Python) — transmission loss from GLORYS12 SSPs
- **callDensity** (R) — detection function fitting and call density estimation

This document inventories the current formats at each interface, identifies
friction points, and frames the key design question for a clean pipeline.

---

## Pipeline overview

```
Audio recordings
    └─ detector (PAMGuard / Koogu / LFDCS)
           └─ MATLAB converters (cnnTableToDetection, lfdcsTableToDetection, ...)
                  └─ detection table (common matchbox/bsnr input format)
                         ├─ [bsnr] ← SNR estimation (option: before matchbox)
                         └─ matchbox adjudication
                                └─ CHT (capture history table)
                                       ├─ [bsnr] ← SNR estimation (option: after matchbox)
                                       └─ callDensity (R)
                                              ├─ fitDetFun → snrDetFun
                                              ├─ cde()
                                              └─ pDetInArea()
                                                     └─ CopernicusPyRAM TL
```

SL (source level) is always external — supplied by the analyst from prior
measurements. It never flows through this pipeline.

---

## Format inventory

### 1. Detection table (post-converter, pre-matchbox)

Common format after MATLAB converter scripts. Required fields for matchbox
and bsnr input:

```
t0, tEnd        detection start/end (MATLAB datenum or datetime; datetime preferred)
fLow, fHigh     frequency band in Hz
soundFolder     path to wav files
channel         recording channel index
duration        detection duration in seconds
```

Each detector's native format differs; converters normalise to this schema.
Converters: `cnnTableToDetection.m`, `lfdcsTableToDetection.m`,
`standardizeDetectorLocale.m` (in `callDensitySupport`).

### 2. matchbox CHT output

One row per matched event. All input columns are carried forward with
per-observer suffixes:

```
key                      event identifier (1..nEvents)
t0, tEnd                 event envelope (min start, max end over observers)
fLow, fHigh              event envelope (min fLow, max fHigh over observers)
detect_observerK         logical, true if observer K detected this event
<col>_observerK          every column of observer K's input table, suffixed
                         NaN/missing where observer K did not detect the event
```

SNR columns (`snr_observerK`, `signalRMSdB_observerK`, `noiseRMSdB_observerK`)
are present only if bsnr was run **before** matchbox and those results were
included in the input detection tables.

### 3. bsnr output (vector/table mode)

Indexed to input annotation table. Output rows align with input rows and
can be joined directly to a detection table or CHT.

```
snr                      SNR in dB
signalRMSdB              signal RMS level in dB
noiseRMSdB               noise RMS level in dB
noiseVar                 noise power variance
```

When `params.calibration` is provided, also includes:
```
signalBandLevel_dBuPa    calibrated signal level (dB re 1 µPa)
noiseBandLevel_dBuPa     calibrated noise level  (dB re 1 µPa)
```

bsnr can be run at two points in the pipeline (see friction points below).

### 4. CopernicusPyRAM TL output

```
data/output/TL_{SITE}_{SPECIES}_GLORYS12_{SEASON}_{FREQ}Hz.csv
    range_m      range in metres
    t1..tN       TL in dB, one column per radial transect

data/output/metadata/TL_{SITE}_{SPECIES}_GLORYS12_{SEASON}_{FREQ}Hz_metadata.csv
```

This is the format `pDetInArea()` and `cde()` expect for the `TL` argument.

### 5. callDensity inputs

#### SL — source level (always external)
```r
list(mean = ..., sd = ..., sampleSize = ...)
```
Supplied by the analyst from prior measurements. Never derived from the
pipeline.

#### NL — noise level (derived or supplied)
```r
list(mean = ..., sd = ..., sampleSize = ...)
# or
data.frame(mean = ..., sd = ..., sampleSize = ...)
```
Derived internally by `nlFromSnrInfo()` from SNRinfo if not supplied.

#### SNRinfo — OG format (single detector)
Output of `capHist2snrInfo()`. Used by `fitDetFun()` and internally by
`cde()`.

```
Detected     0/1, automated detector result
CallRL       signal RMS in dB
NoiseRL      noise RMS in dB
SNR          CallRL - NoiseRL
t            datetime
month        month code
season       season code
```

Derived from CHT by filtering to `detect_table1 == 1` (ground truth rows
only). Implicitly single-detector.

#### SNRinfo — CR format (two-observer capture-recapture)
Output of `capHistTosnrInfo()`. Used by `fitDetFun(..., modelType = "vglm")`.

```
detect_table1, detect_table2
signalRMSdB1, signalRMSdB2
noiseRMSdB1, noiseRMSdB2
datetime, month, season
```

SNR averaged across observers. Derived from matchbox CHT via `mchToCR()`
then `capHistTosnrInfo()`.

#### capHistTab — as consumed by `cde()`

Two-observer format with `detect_table1`/`detect_table2` column naming.
Consumed by `falseDiscoveryRate()` and `capHist2snrInfo()` internally.
Must be in this naming convention — matchbox multi-observer CHT cannot go
directly into `cde()` without reshaping via `mchToCR()` first.

```
detect_table1    ground truth observer
detect_table2    automated detector being evaluated
signalRMSdB      signal level in dB
noiseRMSdB       noise level in dB
t, month, season
```

---

## Current converter functions

| Function | Input | Output | Notes |
|---|---|---|---|
| `readCapHist()` | CSV file path | CHT data.frame | Adds time/season columns. Has undocumented TODOs. |
| `mchToCR()` | multi-observer CHT | two-observer CR CHT | Selects two observers, renames to table1/table2. Uses `verdict` as ground truth. |
| `capHist2snrInfo()` | two-observer CHT | SNRinfo (OG) | Season filtering. Expects `signalRMSdB` (single column). |
| `capHistTosnrInfo()` | two-observer CHT | SNRinfo (CR) | No season filter. Expects `signalRMSdB1`, `signalRMSdB2` (per-observer). Averages across observers. |

---

## Friction points

### 1. SNR join timing is uncontrolled

bsnr can be run before or after matchbox adjudication:

- **Before**: SNR estimated per raw detection. Matched events inherit
  per-observer SNR via the `_observerK` suffix mechanism.
- **After**: SNR estimated per matched/adjudicated event. One SNR per row,
  shared across observers.

Both are scientifically valid but produce different column structures in the
CHT. Downstream converters assume one or the other implicitly and will
silently produce wrong results if given the other.

### 2. Two converter functions with overlapping purpose

`capHist2snrInfo` and `capHistTosnrInfo` both convert CHT → SNRinfo but
expect different column schemas and behave differently:

| | `capHist2snrInfo` | `capHistTosnrInfo` |
|---|---|---|
| SNR columns expected | `signalRMSdB` (single) | `signalRMSdB1`, `signalRMSdB2` |
| Season filtering | yes | no |
| Time/season source | `snr$t`, `snr$month`, `snr$season` (pre-computed) | `datetime` column via `time2season()` |

Nothing enforces which to call for which input. Choice is implicit.

### 3. matchbox CHT → `cde()` requires undocumented reshaping

The path from matchbox output to `cde()` input is:

```
matchbox CHT (multi-observer)
    → mchToCR()              # observer selection + rename
        → capHistTosnrInfo() # → SNRinfo for fitDetFun
        → capHist2snrInfo()  # called internally by cde()
        → falseDiscoveryRate() # consumes CHT directly
```

`cde()` calls `capHist2snrInfo()` internally, so the `capHistTab` argument
must already be in two-observer `table1`/`table2` format. This is not
documented in `cde()`'s Rd page and will fail silently with wrong column
names.

### 4. `cde()` calls deprecated functions internally

`cde()` calls `fitSNRdetectionFunc()` when no `snrDetFun` is supplied.
That function is now deprecated in favour of `fitDetFun()`. Should be
updated.

### 5. `capHistTab` schema consumed by `cde()` is ambiguous

It is unclear whether `signalRMSdB` in the `capHistTab` passed to `cde()`
should be a single column (post-averaging) or whether per-observer columns
are expected. `falseDiscoveryRate()` and `capHist2snrInfo()` need auditing
to confirm which schema they actually require.

---

## Key design question for next session

**What is the canonical format at the MATLAB → R boundary?**

Currently this is implicit: "a CSV written by matchbox (after optional
`mchToCR`-equivalent reshaping in MATLAB), read by `readCapHist()`."

Three options:

### Option A — matchbox CHT is the canonical handoff
R receives the full multi-observer table and handles all reshaping.
One new R function replaces `mchToCR` + both `capHistTosnrInfo` variants.
Requires matchbox to always write the full suffixed format.

*Pro*: MATLAB side stays simple; all logic in R where it's easier to test.
*Con*: R function needs to handle variable numbers of observers.

### Option B — two-observer CR table is the canonical handoff
MATLAB does observer selection (equivalent to `mchToCR`) before writing CSV.
R receives a clean two-observer table.

*Pro*: Simpler R side.
*Con*: Observer selection logic lives in MATLAB where it's harder to test
and version-control.

### Option C — define an explicit callDensity interchange format
Specify column names, types, and units that matchbox writes to and
callDensity reads from. Independent of internal representations of either.

*Pro*: Clean contract; both sides can evolve independently.
*Con*: Requires changes on both sides; more upfront design work.

**Recommendation going in**: Option A, with the matchbox CHT as the
handoff format and a single R function `chtToSNRinfo(cht, groundTruth,
observer, season)` replacing the current converter proliferation. The
`groundTruth` argument specifies which observer (or `"verdict"`) to treat
as ground truth; `observer` specifies which detector is being evaluated.

---

## Other notes for next session

- `readCapHist()` has two undocumented TODOs (column checking, format
  conversion) and hardcoded defaults. If Option A is adopted, this function
  becomes a thin `read.csv` + column validation wrapper, which is worth
  keeping but simplifying.
- `time2season()` and `time2monthCode()` are used inconsistently —
  sometimes called inside converters, sometimes expected to have been called
  upstream. Should be called once at the boundary.
- The `season` filtering in `capHist2snrInfo()` and `cde()` is a candidate
  for moving upstream to the boundary converter, so all downstream functions
  receive pre-filtered data.
- `plotCaptureHistory` was logged as a planned public utility for
  `callDensitySupport` — useful context for the interface discussion since
  it would also consume the CHT format.
- Practical starting point for new branch: `git checkout -b interface-standardisation`
  from current `main` after pushing the testthat commits.

---

## Session context

- callDensity v1.0.0, published MER
- bsnr v0.3.1-beta, 13/13 tests passing
- matchbox: active MATLAB development, gallery.m in progress
- CopernicusPyRAM v0.1.0-alpha, shared with Fran Castro and Danielle Harris
- testthat suite: 79 tests, all passing as of this session
