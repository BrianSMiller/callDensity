# Predict detection function from fitted SNR models

Unified prediction interface for SNR-based detection models. Returns
fitted detection probabilities over a sequence of SNR values (or
user-supplied newdata), optionally with 95% confidence intervals.

Supports:

- glm, gam, scam (mgcv/scam-style models)

- VGAM vglm / vgam (posbernoulli.t models)

## Usage

``` r
predictDetFun(model, ...)
```

## Arguments

- model:

  Fitted model from: fitSNRdetectionFunc(), fitSNRvglm(), or
  fitSNRvgam()

- newdata:

  Data.frame with column `SNR`. If NULL, a regular grid is generated
  from the model frame.

- ci:

  Logical. If TRUE, compute 95% confidence intervals.

- npoints:

  If `newdata = NULL`, number of SNR grid points.

- nsim:

  Number of simulations for VGAM confidence intervals.

## Value

A data.frame with:

- SNR

- fit

- lower (optional)

- upper (optional)
