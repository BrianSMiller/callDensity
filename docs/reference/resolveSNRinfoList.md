# Resolve SNRinfo into a named list matching modelNames

Accepts a single shared data.frame (used for every model), a named list
of data.frames (one per model, for models fit to different underlying
data), or NULL (recover it from each model itself via
[`extractSNRinfo()`](https://briansmiller.github.io/callDensity/reference/extractSNRinfo.md))
– and always returns the list shape.

## Usage

``` r
resolveSNRinfoList(SNRinfo, model, modelNames)
```
