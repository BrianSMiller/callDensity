# Predict detection functions for multiple models

Evaluates one or more fitted detection models on a shared SNR grid and
returns a long-format data frame suitable for ggplot or base plotting.

## Usage

``` r
predictDetFunList(models, newdata = NULL, ci = TRUE, npoints = 300, ...)
```

## Arguments

- models:

  Named list of fitted models.

- newdata:

  Optional data.frame with SNR column. If NULL, a grid is created.

- ci:

  Logical; include confidence intervals where available.

- npoints:

  Number of SNR grid points if newdata is NULL.

## Value

data.frame with columns: SNR, model, fit, lower, upper
