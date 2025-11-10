# Calculate interval forecast coverage

This is a wrapper for the function
[`conformalForecast::coverage`](https://xqnwang.github.io/conformalForecast/reference/coverage.html).
Calculates the mean coverage and the ifinn matrix for prediction
intervals on validation set. If `window` is not `NULL`, a matrix of the
rolling means of interval forecast coverage is also returned.

## Usage

``` r
avgCoverage(object, level = 95, window = NULL, na.rm = FALSE)
```

## Arguments

- object:

  An object of class `bb_cvforecast` or `cb_cvforecast`.

- level:

  Target confidence level for prediction intervals.

- window:

  If not `NULL`, the rolling mean matrix for coverage is also returned.

- na.rm:

  A `logical` indicating whether `NA` values should be stripped before
  the rolling mean computation proceeds.

## Value

A list of class `coverage` with the following components:

- mean:

  Mean coverage across the validation set.

- ifinn:

  A indicator matrix as a multivariate time series, where the \\h\\th
  column holds the coverage for forecast horizon \\h\\. The time index
  corresponds to the period for which the forecast is produced.

- rollmean:

  If `window` is not `NULL`, a matrix of the rolling means of interval
  forecast coverage will be returned.
