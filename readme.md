# callDensity
### Estimating the density of animal calls using auxiliary information from the sonar equation.

This R package implements the methods described by Harris (2012) and Castro et al (in prep) for estimating the density of calls produced by baleen whales from detections on a single fixed hydrophone.

The output of this software is an estiamte of call density: i.e. the number of calls per unit time per unit area. Call densities are similar to detection rates (calls per unit time), but are also standardised by unit area. Furthermore, the call densities produced by this software can also be corrected for detector biases (false positives and probability of detection). The software also produces estimates of CV call density, and each of the components used to estimate call density. 

The inputs into this package are a 'capture history table' of detections of calls, and sonar equation information for these detections. 
Capture history tables include reconciled detections from two or more detectors.
The required inputs from the sonar equation includes:
  1) Signal-to-noise ratio (SNR) of each detection. 
  2) Source level, (SL) distribution for the calls (mean & standard deviation)
  3) Noise level, (NL) distribution (mean & standard deviation)
  4) Transmission losses (TL) for the recording site as a Nx2D radial profiles originating at the recording site and ending at a point at which the probability of detection can be safely assumed to be negligible. 

References

Harris, Danielle V. “Estimating Whale Abundance Using Sparse Hydrophone Arrays.” Thesis, University of St Andrews, 2012. https://research-repository.st-andrews.ac.uk/handle/10023/3463.