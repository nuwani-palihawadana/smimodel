# Remove actual values from a data set for recursive forecasting

Appropriately removes existing (actual) values from the specified column
range (lagged response columns) of a given data set (typically a test
set for which recursive forecasting is required).

## Usage

``` r
remove_lags(data, recursive_colRange)
```

## Arguments

- data:

  Data set (a `tibble`) from which the actual lagged values should be
  removed.

- recursive_colRange:

  The range of column numbers in `data` from which lagged values should
  be removed.

## Value

A `tibble`.
