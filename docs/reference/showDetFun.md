# Plot SNR detection function with observed SNR distributions

Visualises one or more fitted SNR detection functions together with the
distribution of detected and missed events. Dispatches on the class of
`model`: a single fitted model (glm/gam/scam/vglm/vgam) draws one curve
(`showDetFun.detFun`); a named list of fitted models overlays all of
them, each in its own colour (`showDetFun.list`).

By default draws a rug along the y=1 (detected) and y=0 (missed)
reference lines. Optionally overlays density curves or histograms of the
detected/missed SNR distribution above and below the main panel, via
patchwork.

## Usage

``` r
showDetFun(model, ...)
```

## Arguments

- model:

  A single fitted detection model from
  [`fitDetFun()`](https://briansmiller.github.io/callDensity/reference/fitDetFun.md),
  or a *named* list of them (names become the legend/colour labels).

- ...:

  Passed to `showDetFun.detFun()` or `showDetFun.list()`.

- SNRinfo:

  Data frame used to fit the model(s). Must contain columns `SNR` and
  `Detected` (logical or 0/1). Optional: if omitted (`NULL`, the
  default), it's recovered directly from each fitted model, since that's
  what it was fit on anyway. For the list form, this can also be a named
  list of data frames, one per model (matched by name) – use this only
  if you want to show *different* data than what a model was actually
  fit on. A single shared data.frame is used for every model otherwise.

- distribution:

  One of `"none"` (default), `"density"`, or `"histogram"`. Adds
  detected/missed SNR distributions above/below the main panel. For the
  list form, one histogram/density per model is overlaid, coloured to
  match its curve.

- rug:

  Logical. If `TRUE` (default), draws tick marks along the y=1 line for
  detected events and along y=0 for missed events. For large datasets, a
  random subsample of `rug.max` points is used.

- rug.max:

  Maximum number of rug ticks per group. Default 500. Set to `Inf` to
  show all points (slow for n \> 2000).

- rug.alpha:

  Alpha for rug ticks. Default 0.4.

- panel.heights:

  Relative heights of the top/main/bottom panels when
  `distribution != "none"`. Default `c(1, 4, 1)`.

- hist.alpha:

  Fill alpha for histogram/density panels. Default 0.4 for a single
  model; for a list of models, default `1/length(model)` so overlapping
  fills stay legible.

- palette:

  Optional named vector of colours, one per model (list form only).
  Default `NULL` uses ggplot2's default discrete palette, applied
  identically across all three panels so a model's line, histogram, and
  rug always match – even if that model has, say, zero missed events and
  would otherwise be silently dropped from just the bottom panel's own
  colour scale.

- npoints:

  Number of SNR grid points for the fitted curve(s). Default 300.

- xlim:

  Optional length-2 numeric vector restricting the SNR range shown.
  Applied consistently to every panel (and to the fitted curve's own
  prediction grid) at construction time, so the top/main/bottom panels
  can't drift out of sync the way they would if a scale were changed on
  just one panel after the fact – see Details.

- show.counts:

  Logical, default `TRUE`. For a single model, adds a subtitle reporting
  the number of detections out of total observations. For a list of
  models, a subtitle doesn't scale past one model, so instead each
  legend entry gets its own `(detections/events)` count appended, with
  the legend titled accordingly – e.g. "GLM (62/800)" under a "Model (#
  detections / \# events)" legend title.

## Value

A ggplot object, or a patchwork object when `distribution != "none"`.

## Details

When `distribution != "none"`, the return value is a patchwork object
made of three independently-built ggplot panels (top/main/bottom)
sharing one x-range fixed at construction time – patchwork doesn't keep
panel axes "linked" the way, say, MATLAB's `linkaxes` does; each panel
is just a static plot glued to its neighbours at render time. So
rescaling only one panel after the fact (e.g.
`p[[2]] <- p[[2]] + xlim(...)`) only affects that one panel, and the
others silently keep their original range. Pass `xlim` here instead, so
every panel is built from the same range to begin with.

Adding a plain ggplot layer (`geom_line()`, `theme()`, etc) to that
patchwork object with `+` applies it to whichever panel was supplied
*last* to
[`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
– which is deliberately the main panel here (supplied after top and
bottom, then placed back in the visual middle via a `design` string),
since that's almost always the panel people want to add
curves/annotations to. So `p + geom_line(...)` lands on the main panel
directly – no indexing needed. To target the top or bottom panel
specifically instead, index it directly by its position in
`wrap_plots()`'s argument list (which is top, bottom, main – *not* their
top-to-bottom visual order): `p[[1]] <- p[[1]] + ...` for the top panel,
or `p[[2]] <- p[[2]] + ...` for the bottom one.

`+` only ever reaches that one last-supplied panel, regardless of what's
being added – there's no way to make some additions (like a post-hoc
`xlim()`) broadcast to all three panels while others stay targeted at
just the main one, using `+` for both. For anything that should apply to
every panel at once (an `xlim()` change included), use patchwork's own
broadcast operator `&` instead: `p & ggplot2::xlim(-10, 10)`. That said,
changing the range post hoc this way re-clips already-built panels
rather than rebuilding them, so prefer passing `xlim` up front where
possible.
