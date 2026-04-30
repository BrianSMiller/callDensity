# callDensity
  <!-- badges: start -->
  [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![DOI](https://zenodo.org/badge/638054659.svg)](https://doi.org/10.5281/zenodo.19932510)
  <!-- badges: end -->

### Estimating the density of animal calls using auxiliary information from the sonar equation.

This R package implements the methods described by Harris (2012) and Castro et al (2024) for estimating the density of calls produced by baleen whales from detections on a single fixed hydrophone.

The output of this software is an estiamte of call density: i.e. the number of calls per unit time per unit area. Call densities are similar to detection rates (calls per unit time), but are also standardised by unit area. Furthermore, the call densities produced by this software can also be corrected for detector biases (false positives and probability of detection). The software also produces estimates of CV call density, and each of the components used to estimate call density.

The inputs into this package are a 'capture history table' of detections of calls, and sonar equation information for these detections. Capture history tables include reconciled detections from two or more detectors. The required inputs from the sonar equation includes: 1) Signal-to-noise ratio (SNR) of each detection. 2) Source level, (SL) distribution for the calls (mean & standard deviation) 3) Noise level, (NL) distribution (mean & standard deviation) 4) Transmission losses (TL) for the recording site as a Nx2D radial profiles originating at the recording site and ending at a point at which the probability of detection can be safely assumed to be negligible.

## Installation instructions

Install directly from GitHub:

``` r
devtools::install_github("BrianSMiller/callDensity")
```

Alternatively, clone the repository and install locally:

``` r
devtools::install('path/to/callDensity/')
```

## Usage

Example usage of callDensity can be found in the vignettes that are included in the package documentation. Vignettes can be installed by including build_vignettes = TRUE during the installation:

``` r
devtools::install_github("BrianSMiller/callDensity", build_vignettes = TRUE)
```

## References

1.  Harris, DV. *Estimating Whale Abundance Using Sparse Hydrophone Arrays.* PhD Thesis, University of St Andrews, 2012. <https://research-repository.st-andrews.ac.uk/handle/10023/3463>.
2.  Castro, FR et al.. *Beyond Counting Calls: Estimating Detection Probability for Antarctic Blue Whales Reveals Biological Trends in Seasonal Calling.* Frontiers in Marine Science 11 (2024). <https://doi.org/10.3389/fmars.2024.1406678>.
3.  Harris, DV et al. *Estimating the Detection Probability of Long-Ranging Baleen Whale Song Using a Single Sensor: Towards Density Estimation.* The Journal of the Acoustical Society of America 158, no. 6 (2025): 4582–93. <https://doi.org/10.1121/10.0036892>.
4.  Miller, BS et al. (in press). Common ground: efficient, consistent, observer-independent bioacoustic call density estimation with adjudicated ground truth and capture-recapture detection functions. *Methods in Ecology and Evolution*. doi: [to be assigned]
