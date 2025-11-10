# Obtaining forecasts on a test set from a fitted `pprFit`

Gives forecasts on a test set.

## Usage

``` r
# S3 method for class 'pprFit'
predict(
  object,
  newdata,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL,
  ...
)
```

## Arguments

- object:

  A `pprFit` object.

- newdata:

  The set of new data on for which the forecasts are required (i.e. test
  set; should be a `tsibble`).

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string. (Since the nonlinear
  functions are estimated using splines, extrapolation is not desirable.
  Hence, if any predictor variable in `newdata` that is treated
  non-linearly in the estimated model, will be truncated to be in the
  in-sample range before obtaining predictions. If any variables are
  listed here will be excluded from such truncation.)

- recursive:

  Whether to obtain recursive forecasts or not (default - `FALSE`).

- recursive_colRange:

  If `recursive = TRUE`, the range of column numbers in `newdata` to be
  filled with forecasts. Recursive/autoregressive forecasting is
  required when the lags of the response variable itself are used as
  predictor variables into the model. Make sure such lagged variables
  are positioned together in increasing lag order (i.e.
  `lag_1, lag_2, ..., lag_m`, `lag_m =` maximum lag used) in `newdata`,
  with no break in the lagged variable sequence even if some of the
  intermediate lags are not used as predictors.

- ...:

  Other arguments not currently used.

## Value

A `tsibble` with forecasts on test set.
