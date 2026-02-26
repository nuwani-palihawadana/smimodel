# Futures through single season block bootstrapping

Generates possible future sample paths by applying the single season
block bootstrap method.

## Usage

``` r
blockBootstrap(
  object,
  newdata,
  resids,
  preds,
  season.period = 1,
  m = 1,
  num.futures = 1000,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL
)
```

## Arguments

- object:

  Fitted model object.

- newdata:

  Test data set. Must be a data set of class `tsibble`.

- resids:

  In-sample residuals from the fitted model.

- preds:

  Predictions for the test set (i.e. data for the forecast horizon).

- season.period:

  Length of the seasonal period.

- m:

  Multiplier. (Block size = `season.period * m`)

- num.futures:

  Number of possible future sample paths to be generated.

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string.

- recursive:

  Whether to obtain recursive forecasts or not (default - FALSE).

- recursive_colRange:

  If `recursive = TRUE`, The range of column numbers in test data to be
  filled with forecasts.

## Value

A matrix of simulated future sample paths.
