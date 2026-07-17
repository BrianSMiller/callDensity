# callDensity: noise levels, CR, and SNR truncation — handover (2026-07-17)

Supersedes the NL estimator handover of 2026-07-16, which contained claims the
tests did not support.

---

## What is settled

### 1. `nlFromDetections` works. `nlFromSnrInfo` does not generalise.

`nlFromSnrInfo` corrects the noise bias by adding the SNR at which the
detection function reaches 0.5. That quantity is a property of the detector.
The bias is a property of the ocean and the geometry.

Measured across a grid of detector locations and scales:

| | mean \|NL error\| | range |
|---|---|---|
| naive (no correction) | 2.98 dB | — |
| `nlFromSnrInfo` | 2.01 dB | −3.99 to +1.86 |
| `nlFromDetections` | **0.58 dB** | −1.35 to +0.98 |

`nlFromSnrInfo`'s error tracks detector 2's location one for one. That is fatal
for a general estimator. It may still be correct for the narrow case it was
written for; that case was never defined and never tested.

`nlFromDetections` inverts the `pDetGivenNL` weighting with `uniroot`. Its
residual tracks detector *scale*, not location, and stays under about 1 dB.

Under CR, with the union curve as the selecting filter, the residual is
−0.27 dB mean, −0.66 to +0.33 across 32 cells.

### 2. `nlFromDetections` needs no new argument.

`snrDetFun` is used for exactly one thing inside it: the weighting in
`predictSampledNL`. It is never used as "the fitted curve". So under CR you
pass `p_union` directly. `pDetGivenNL` has an `is.function()` branch, so a plain
closure works.

`p_union` is already computed inside VGAM. `predictDetFun.vglm` calls
`predictvglm(type.fitted = "onempall0")` to unconditional-ise the per-observer
probability, but does not return it. One call gets it:

```r
pUnionFromVglm <- function(model, snrGrid = seq(-100, 100, by = 0.1)) {
  pu <- VGAM::predictvglm(model, newdata = data.frame(SNR = snrGrid),
                          type = "response", type.fitted = "onempall0")
  stats::approxfun(snrGrid, as.vector(pu), rule = 2)
}
```

### 3. NL is second order. It was never the blocker.

Under CR, estimating the noise rather than being given it costs 0.2% of `Dc` at
fdr 0 and 5% at fdr 0.1. Under OG the same substitution cost 12% to 39%,
because OG has no access to the curve that selected its sample.

The blocker is `p_a`. See below.

---

## What CR fixes, and what it does not

CR is not a code path in `cde`. `cde` always treats `detect_table1` as truth.
OG puts detector 1's raw detections there. CR puts the adjudicated verdict
there. `mchToCR` does exactly this with real data:

```r
names(ch)[names(ch)=='verdict'] <- 'detect_table1'
```

### `c`: CR fixes it completely.

Simulated with **zero** false positives anywhere:

| method | estimated `c` (32 cells) |
|---|---|
| OG | mean 0.46, range 0.02 to 0.79 |
| CR | **0.0000, every cell** |

Simulated with a true false discovery rate of 0.1:

| method | estimated `c` |
|---|---|
| OG | mean 0.51, range 0.12 to 0.82 |
| CR | **0.0999, every cell** |

CR does not merely report zero when zero is true. It recovers a non-zero rate
correctly. This is the single largest result of the session.

Under OG, detector 1 is assumed correct and is not, so every real call detector
2 finds and detector 1 misses is recorded as a false positive. That flows into
`Dc` through `(1-c)`.

### OG's detection function is not the problem.

Given SNR the two detectors are independent, so fitting detector 2's curve on
detector 1's detections is unbiased. With the true noise supplied, OG's `p_a`
is within 2% of truth. OG's entire `Dc` error is the spurious `c`.

### The VGLM is load-bearing, not decoration.

Ablation: hand OG's ordinary binomial model the correct verdict, changing
nothing else. `p_a` comes out at 0.77 against a truth of 0.12.

The CR sample is calls at least one observer flagged. A call both missed is not
in it. At low SNR the sample therefore holds only detector 2's successes, so a
glm sees roughly `p2/(p1+p2)` rather than `p2`. `posbernoulli.t` conditions on
detection by at least one observer, which is how the sample arose.

Fixing the verdict alone makes things **worse** than OG.

### `p_a`: CR does not fix it, and this is now the blocker.

`p_a` is the mean detection probability over all calls. Nearly every call is
far away and faint, so `p_a` is decided almost entirely by the low-SNR tail.
The sample contains only calls somebody heard, which are the loud ones. The
curve is fitted where the data is and evaluated where it is not.

Sweep over `det1scale` at 0.5, 1, 2, 4, 8, ten seeds per cell:

| det1scale | CR paFraction (mean) | oracle |
|---|---|---|
| 0.5 | 1.54 to 1.70 | 0.99 |
| 1 | 1.08 to 1.16 | 0.99 |
| 2 | 1.02 to 1.08 | 1.00 |
| 4 | 1.02 to 1.07 | 1.00 |
| 8 | **1.00 to 1.02** | 1.00 |

`Dc` is biased low by the reciprocal: 0.59 at det1scale 0.5.

### The mechanism is recapture coverage, not sample coverage.

This is the sharpest result and it took several wrong turns to reach.

At `det2scale = 1`, `snrP05` (5th percentile of the CR sample) moves from −1.7
to −16.1 as `det1scale` rises, and the bias collapses. Consistent with sample
coverage.

At `det2scale = 4`, `snrP05` barely moves (−11.6 to −15.9) and the bias still
collapses from 1.70 to 1.00. Detector 2 is shallow, so it reaches down to −11 dB
on its own and the sample **has** data there. But detector 1 never saw those
calls, so there are **no recaptures**, and a capture-recapture model identifies
nothing from a region where only one observer contributes.

**The requirement is that both observers contribute detections across the SNR
range that matters. It is a property of the partner's slope, not its recall.**

Detector 1's recall is knowable in the simulation (it equals `truePa` for the
matching detector 2 configuration). Sorted by recall the error is not even
monotone. Two configurations with recall 0.089 and 0.097 differ six-fold in
error; the one that works is the shallow one.

This refines the Common Ground discussion (Miller et al. 2026, §4), which
attributes the effect to low recall. The mechanism there is right ("did not
provide sufficient information to accurately model the heterogeneity in SNR
across the full range of that covariate... effectively extrapolated"). The
proxy is wrong. Recall is the height of the partner's curve; what matters is
its slope.

Practical consequence: **a mediocre, gradually-degrading partner is better for
CR than a sharp one, however good its recall.** An SNR-insensitive CNN (e.g.
the D-call DNN of Miller et al. 2022) is close to an ideal partner. A hard
energy threshold is close to a worst case.

### `cde`'s Monte Carlo estimator is sound.

The oracle fits the detection function to detector 2's own detections over
*every call that was made*, and is given the true noise. It is not a method:
that sample cannot exist. It is the control, differing from CR in exactly one
respect, SNR coverage.

Oracle `paFraction`: 0.99 to 1.00 in all 10 configurations, 10 seeds each.
Oracle `DcFraction`: 0.99 to 1.01.

`outerloop = 10` is adequate. Same data, same fit, ten replicates of `cde`:
0.112 to 0.121, CV about 3%. `outerloop = 200` gives the same answer. Monte
Carlo noise is second decimal and is not the problem.

### `CV.pa` is blind to all of this.

`CV.pa` sat at **0.501 to 0.504 across all 200 runs**, oracle and CR,
`det1scale` 0.5 to 8, while the actual error in `p_a` ran from 1% to 100%.

`pa_CV` is built from `sd_pa_t`, the per-transect SD of pDet **across cells**.
That is the spatial spread of detection probability, near 1 at the recorder and
near 0 far away. It is set by geometry. It is not a standard error of `p_a`.

This is not a bug: the paper describes `CV.pa` as combining spatial and Monte
Carlo uncertainty, and that is what it does. But it means **the dominant error
in `p_a` is not in the confidence interval**, and nothing in `cde`'s output
moves when `p_a` is 70% wrong. This warrants a conversation with Danielle
Harris and Nat Kelly, not a TODO.

Likely connection: unmodelled heterogeneity in double-observer mark-recapture
distance sampling, and the point-independence treatment of Laake & Borchers
(2004). Not verified. Danielle's ground.

---

## SNR truncation: audit and plan

### Why it is the fix

Zeroing `p` below θ multiplies the extrapolated part of the curve by zero. It
cannot contribute. What does the work instead is q(θ), the fraction of calls
arriving above θ, which comes from SL, TL and NL.

This creates no information about faint calls. It relocates the assumption from
the detector's behaviour at SNRs nobody ever observed, which is unmeasurable in
principle, to the sonar equation, which is already the backbone of the method
and is independently checkable. Cost: `Nc` shrinks, so precision falls.

### The algebra, and why NA vs 0 decides the estimand

Let q(θ) = P(received SNR > θ) for a call uniformly in the area.

- **NA** (current): `p_a = E[p | snr > θ]`, renormalised. Then `Dc_hat = Dc * q(θ)`.
  You get the density of **above-threshold calls**.
- **0**: `p_a = E[p * 1(snr > θ)] = q(θ) * E[p | snr > θ]`. The q(θ) cancels.
  You get the density of **all calls**.

The 2026-07-16 note (line 377 of `callDensity_snrThreshold.Rmd`) proposed
tracking the excluded proportion and scaling back. That proportion is q(θ), and
zeroing obtains it for free. No bookkeeping.

**The two truncations must not be treated alike.** Distance truncation stays NA,
because `A` is redefined as `pi*w^2` to match. SNR truncation must be 0, because
`A` is not redefined.

### Where each term goes

| term | what changes |
|---|---|
| `Nc` | count only detections with measured SNR >= θ |
| `c` | adjudicate only detections with SNR >= θ; FDR among those |
| detection function | fit on the truncated sample; above θ it has recaptures, so it interpolates |
| `p_a` | zero p below θ, do not subset. Distance truncation stays NA. |
| `A` | unchanged |
| NL | estimated from the **un**truncated sample. NL is a property of the ocean. |

### Audit: six points of contact

| # | where | status |
|---|---|---|
| 1 | `pDetInArea:289` sets `NA` not `0` | wrong estimand |
| 2 | `Nc` | not handled at all; `cde` takes it as a number and never truncates it |
| 3 | `countDetections:280` | checks `snrColName %in% colnames(det)` then hardcodes `det$snr` |
| 4 | `falseDiscoveryRate` | works, but defaults `snrColName='SNR'` while `countDetections` uses `'snr'` |
| 5 | `cde:171` fits on truncated `SNRinfo` | correct, **but only when `snrDetFun` is NULL** |
| 6 | `cde:179` estimates NL from untruncated `SNRinfo` | correct, and deliberately commented |

Plus: `cdeFromParamFile` (line 736+) duplicates every one of these.

### The fixes

**1. `pDetInArea:289`** — split the two truncations:

```r
        # Distance truncation: outside w, not part of the estimand -> NA
        predmatrix[is.na(allTrunc[, j])] <- NA
        # SNR truncation: inside the estimand, undetectable by fiat -> 0
        predmatrix[newd$SNR < snrTruncationThreshold] <- 0
```

**2. `Nc` guard in `cde`.** `cde` never sees the detections, so it cannot fix
this. It can refuse:

```r
  if (snrTruncationThreshold != -Inf && !isTRUE(NcIsTruncated)) {
    stop("snrTruncationThreshold is set, so Nc must count only detections with SNR >= the threshold. Pass NcIsTruncated = TRUE to confirm.")
  }
```

**3. `countDetections:280`:**

```r
      det <- subset(det, det[[snrColName]] >= snrTruncationThreshold)
```

Column-name trap: `capHistTab` carries `SNR` (callDensitySim.R:274); raw
detections carry `snr` (callDensitySim.R:75). The current code fails silently:
`det$snr` returns NULL, `subset(det, NULL >= 0)` returns zero rows, `Nc = 0`.

**4. Tests.** `pDetGivenNL` is validated against a half-million-call brute force
simulation sharing no code with the integral, agreeing to 2%. That is the only
test in the pile that can catch a wrong idea rather than a typo. Extend the same
pattern to truncation: apply the threshold by subsetting in the simulation and
by zeroing p in the integral, and confirm they agree. Then a `Dc` recovery test:
truncate at θ, confirm `Dc_hat` returns all-call density, not q(θ) times it.

### Choosing θ

θ should be set where **recaptures** run out, not where the sample runs out.
That is measurable: the SNR range over which both observers contribute
detections.

---

## Open, and not for this branch

- **`nlFromDetections`'s signature is not final.** Under truncation the
  selecting filter becomes p_union **and** the threshold. Do not finalise the NL
  interface before truncation lands.
- **`capHist2snrInfo` cannot feed a VGLM.** It returns `Detected`, `CallRL`,
  `NoiseRL`, `SNR`, `t`, `month`, `season`, dropping the per-observer columns.
  So `cde` can never fit a vglm internally, and every CR analysis must
  hand-build its sample. This is why the real-data vignette has been stuck.
  Interface-standardisation branch.
- **The noise measured at detections does two jobs.** It is a covariate in
  fitting `p(SNR)` and it parameterises the noise distribution `pDetInArea`
  draws from. Those are different quantities with different biases, and tying
  them together is the only reason `nlFromDetections` needs the selecting curve
  at all. `noiseLevelDistribution` reads an independent noise table; where one
  exists, nothing needs correcting. Whether an independent NL can be made
  commensurate with the SNRs the detection function was fitted on is a separate
  and probably more important question: SNR = RL - NL only holds if NL is
  measured the same way RL was, same band, same averaging duration.
- `class(snr)=='character'` in `capHist2snrInfo:299`; use `is.character()`.
- `capHist2snrInfo:311` computes `dup <- duplicated(snr$key)` and never uses it.
  `key` is never created by `simsTocaptureHistoryTable`, so it is
  `duplicated(NULL)`. Fossil of pre-matchbox sequential pairwise matching.
- `future.apply` used but undeclared in DESCRIPTION. Will bite the parallel path.
- `readme.html` and `.github` trip `R CMD check` NOTEs.

---

## Method notes for future sweeps

**Force `RNGkind`, not just the seed.** `future.seed = TRUE` switches the
generator to L'Ecuyer-CMRG, so `set.seed(42)` in a worker seeds a *different*
generator and produces a different realisation. Serial and parallel then
disagree by tens of percent and it looks like a bug. Forcing
`RNGkind("Mersenne-Twister")` inside the simulation function, before
`set.seed()`, makes every cell deterministic end to end, because `cde`'s
internal `mvrnorm` draws inherit the generator too. Verified: serial and
parallel produce identical results.

**One realisation is one ticket.** All cells sharing `set.seed(42)` makes cells
comparable but gives the sweep a single draw in the dimension that matters. Two
draws of the same 16-cell grid gave `paFraction` ranges of 0.65 to 1.49 and 0.82
to 1.16. Seed must be a grid dimension.

**Parallelise the cells, not the inner loop.** `cde(parallel = FALSE)`; cells are
independent and CPU-bound, and scale nearly linearly. 32 cells across 16 workers
runs in under a minute.

**Never `library(callDensity)` in an experiment script.** It loads the installed
build, not the source tree, so the script silently tests code that is not on the
branch. `devtools::load_all()`.

**Report truth beside estimate.** `trueC` next to `c`, `truePa` next to `pa`.
`truePa` depends only on detector 2, so any movement across detector 1
configurations is estimation error. Without that column the sweep is
uninterpretable.
