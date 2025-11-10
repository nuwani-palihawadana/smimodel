# Calculate interval forecast width

This is a wrapper for the function
[`conformalForecast::width`](https://xqnwang.github.io/conformalForecast/reference/width.html).
Calculates the mean width of prediction intervals on the validation set.
If `window` is not `NULL`, a matrix of the rolling means of interval
width is also returned. If `includemedian` is `TRUE`, the information of
the median interval width will be returned.

## Usage

``` r
avgWidth(
  object,
  level = 95,
  includemedian = FALSE,
  window = NULL,
  na.rm = FALSE
)
```

## Arguments

- object:

  An object of class `bb_cvforecast` or `cb_cvforecast`.

- level:

  Target confidence level for prediction intervals.

- includemedian:

  If `TRUE`, the median interval width will also be returned.

- window:

  If not `NULL`, the rolling mean (and rolling median if applicable)
  matrix for interval width will also be returned.

- na.rm:

  A logical indicating whether `NA` values should be stripped before the
  rolling mean and rolling median computation proceeds.

## Value

A list of class `width` with the following components:

- width:

  Forecast interval width as a multivariate time series, where the
  \\h\\th column holds the interval width for the forecast horizon
  \\h\\. The time index corresponds to the period for which the forecast
  is produced.

- mean:

  Mean interval width across the validation set.

- rollmean:

  If `window` is not `NULL`, a matrix of the rolling means of interval
  width will be returned.

- median:

  Median interval width across the validation set.

- rollmedian:

  If `window` is not `NULL`, a matrix of the rolling medians of interval
  width will be returned.
