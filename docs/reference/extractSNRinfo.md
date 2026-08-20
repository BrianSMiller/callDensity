# Recover an SNRinfo-shaped data.frame directly from a fitted model

Every model type
[`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md)
produces is fit against SNR and a detected/missed outcome, so both are
still there to pull back out – this is what makes the separate `SNRinfo`
argument optional.

## Usage

``` r
extractSNRinfo(model)
```

## Arguments

- model:

  A fitted detFun object.

## Value

A data.frame with columns `SNR` and `Detected`.

## Details

glm/gam/scam are fit via a plain formula against a data.frame
(`Detected ~ SNR`), so
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) already
returns exactly those two columns – the same mechanism
[`predictDetFunList()`](https://briansmiller.github.io/callDensity/reference/predictDetFunList.md)
already relies on for its default SNR grid.

vglm's capture-recapture formula uses a multi-observer response matrix
rather than a single `Detected` column, so it needs its own path: SNR
comes from the `@x` design matrix and Detected from the `@y` response
matrix, picking out `whichObserver`'s column – the same slot
`predictDetFun.vglm()` already uses for its own SNR range. Two list
entries can share one underlying vglm fit and differ only in
`@extra$whichObserver` (as in the CommonGround vignette, comparing
observer 1 vs observer 2 on the same model) – this reads whichObserver
fresh from each model object, so that distinction is preserved
correctly.
