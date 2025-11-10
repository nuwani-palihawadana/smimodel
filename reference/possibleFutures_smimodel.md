# Possible future sample paths (multi-step) from `smimodel` residuals

Generates possible future sample paths (multi-step) using residuals of a
fitted `smimodel` through recursive forecasting.

## Usage

``` r
possibleFutures_smimodel(
  object,
  newdata,
  bootstraps,
  exclude.trunc = NULL,
  recursive_colRange
)
```

## Arguments

- object:

  A `smimodel` object.

- newdata:

  The set of new data on for which the forecasts are required (i.e. test
  set; should be a `tsibble`).

- bootstraps:

  Generated matrix of bootstrapped residual series.

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string.

- recursive_colRange:

  The range of column numbers in `newdata` to be filled with forecasts.

## Value

A list containing the following components:

- firstFuture:

  A `numeric` vector of 1-step-ahead simulated futures.

- future_cols:

  A list of multi-steps-ahead simulated futures, where each list element
  corresponds to each 1-step-ahead simulated future in `firstFuture`.
