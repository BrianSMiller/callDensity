# Convert a capture history table to SNRinfo format

Canonical converter from a capture history table to the 'SNRinfo' format
used throughout callDensity for detection-function fitting
([`fitDetFun`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md))
and
[`cde`](https://briansmiller.github.io/callDensity/reference/cde.md).
Replaces `capHist2snrInfo`/`capHistTosnrInfo` (now deprecated thin
wrappers around this function – see
[`mchToCR`](https://briansmiller.github.io/callDensity/reference/mchToCR.md),
also deprecated, for the older 2-observer-reduction route these
superseded).

Supports both:

- **OG** (observer-ground-truth): a single `observers` column, treated
  as the detector under evaluation against a (possibly imperfect)
  `groundTruth` column. This is the original, degenerate-format shape
  (Danielle Harris' original script): every retained row is treated as a
  real call, and `Detected` records whether one automated detector also
  caught it.

- **CR** (capture-recapture): two or more `observers` columns, fit
  jointly via
  `fitDetFun(modelType = 'vglm', yColNames = c(groundTruthCol, observerCols))`.
  Any number of observers is supported –
  [`VGAM::posbernoulli.t`](https://rdrr.io/pkg/VGAM/man/posbernoulli.t.html)
  is not limited to two occasions, despite most examples in this package
  using two.

Row selection is a union filter across `groundTruth` and all `observers`
columns (flagged by at least one), followed by restriction to
`groundTruth == 1` (i.e. adjudicated/assumed-true calls only). This is
the sample construction `posbernoulli.t` (the CR/vglm case) assumes
generated the data – it conditions on detection by at least one occasion
internally, so handing it a sample that was itself assembled that way is
correct, not circular. The same construction is also correct, and
simpler, for the OG/glm/gam/scam case. What is NOT safe is filtering to
`groundTruth == 1` and then fitting an ordinary `glm` on a `vglm`-style
union-filtered sample – `glm` has no equivalent internal correction for
that selection mechanism, and doing so silently biases the fit (see
`notes/cde_nl_and_truncation_handover.md`'s ablation section for a
worked example of exactly this mistake).

## Usage

``` r
chtToSNRinfo(
  cht,
  groundTruth = "verdict",
  observers,
  signalCol = "signalRMSdB",
  noiseCol = "noiseRMSdB",
  timeCol = "t0",
  season = "year"
)
```

## Arguments

- cht:

  Capture history table. Matchbox-native column names (`t0`,
  `signalRMSdB`, `noiseRMSdB`, `verdict`, `detect_observerN`) are
  assumed by default. Legacy `detect_table1`/`detect_table2` naming also
  works directly without any renaming step – e.g.
  `groundTruth = 'table1'`, `observers = 'table2'` reproduces the old
  `capHist2snrInfo`/`capHistTosnrInfo` row selection exactly (a
  single-element `observers` degenerates algebraically to their
  `groundTruth == 1` filter, since the union filter with a single
  observer is redundant with the ground-truth restriction that follows
  it).

- groundTruth:

  Column name (or bare suffix, resolved as `detect_<groundTruth>` if not
  found directly) to treat as ground truth. Default `'verdict'`.

- observers:

  Character vector of one or more column names (or bare suffixes, same
  resolution as `groundTruth`) for the detector(s) under evaluation. A
  single element reproduces the old OG/SNRinfo shape; two or more is a
  genuine multi-observer CR sample.

- signalCol, noiseCol:

  Column(s) supplying CallRL/NoiseRL. Default
  `'signalRMSdB'`/`'noiseRMSdB'` (matchbox-native: one shared value per
  matched event, SNR measured after matching). A length-2 vector instead
  averages a per-observer pair (legacy shape, SNR measured before
  matching – see `notes/cde_pipeline_handover.md`'s "SNR join timing"
  note).

- timeCol:

  Column to derive `t`/`month`/`season` from. Default `'t0'`
  (matchbox-native, a Matlab datenum, converted via
  [`mat2Rdate`](https://briansmiller.github.io/callDensity/reference/mat2Rdate.md)).
  Any other column is used as-is (assumed already POSIXct-like) – e.g.
  `timeCol = 't'` for the legacy pre-computed-time shape.

- season:

  A timeCode corresponding to months, seasons, or 'year' (the default,
  which returns all rows).

## Value

SNRinfo data.frame: the resolved `groundTruth` and all `observers`
columns (kept intact, unrenamed – ready to pass as `yColNames` to
`fitDetFun(modelType = 'vglm')`), plus `Detected` (the last-named
observer, for the single-observer OG case and backward compatibility),
`CallRL`, `NoiseRL`, `SNR`, `t`, `month`, and `season`.
