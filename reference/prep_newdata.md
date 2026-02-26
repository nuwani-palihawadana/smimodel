# Prepare a data set for recursive forecasting

Prepare a test data for recursive forecasting by appropriately removing
existing (actual) values from a specified range of columns (lagged
response columns) of the data set. Handles seasonal data with gaps.

## Usage

``` r
prep_newdata(newdata, recursive_colRange)
```

## Arguments

- newdata:

  Data set to be ared. Should be a `tsibble`.

- recursive_colRange:

  The range of column numbers (lagged response columns) in `newdata`
  from which existing values should be removed. Make sure such columns
  are positioned together in increasing lag order (i.e.
  `lag_1, lag_2, ..., lag_m`, `lag_m =` maximum lag used) in `newdata`,
  with no break in the lagged variable sequence even if some of the
  intermediate lags are not used as predictors.

## Value

A `tibble`.
