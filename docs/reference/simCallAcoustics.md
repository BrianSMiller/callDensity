# Simulate acoustic properties of calls

Simulate acoustic properties of calls

## Usage

``` r
simCallAcoustics(
  sim,
  SL = data.frame(mean = 190, sd = 4, sampleSize = 350),
  NL = data.frame(mean = 84, sd = 4, sampleSize = n),
  TL = function(r) {
20 * log10(r)
 }
)
```

## Arguments

- sim:

  - data.frame containing simulated call locations (from
    simCallLocations). Each location must have a distance.

- SL:

  - Source level distribution. SL is a data.frame with columns named
    'mean' and 'sd' that describe the mean and standard deviation of the
    distribution (normal in dB) from which source levels will be
    generated for each simulated call.

- NL:

  - Noise level distribution. NL is a data.frame with columns named
    'mean' and 'sd' that describe the mean and standard deviation of the
    distribution (normal in dB) of noise levels from which noise levels
    will be generated for each simulated call

- TL:

  - a function that takes a single argument, vector r, and returns
    transmission losses for the ranges in that vector. The default
    function is spherical spreading: Default = function(r)20\*log(r)

## Value

Simulation data.frame including columns that contain simulated source
levels (SL) noise levels (noiseRMSdB), and transmission losses (TL).
