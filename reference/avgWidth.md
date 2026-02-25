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

## Examples

``` r
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1055
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(n)) |>
  mutate(
    # Add x_lags
    x_lag = lag_matrix(x_lag_000, 5)) |>
  unpack(x_lag, names_sep = "_") |>
  mutate(
    # Response variable
    y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
    # Add an index to the data set
    inddd = seq(1, n)) |>
  drop_na() |>
  select(inddd, y, starts_with("x_lag")) |>
  # Make the data set a `tsibble`
  as_tsibble(index = inddd)

# Training set
sim_train <- sim_data[1:1000, ]
# Test set
sim_test <- sim_data[1001:1050, ]

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Model fitting
pprModel <- model_ppr(data = sim_train,
                      yvar = "y",
                      index.vars = index.vars)
#> [1] "model 1"
                      
# Conformal bootstrap prediction intervals (2-steps-ahead interval forecasts)
set.seed(12345)
pprModel_cb <- cb_cvforecast(object = pprModel,
                             data = sim_data,
                             yvar = "y",
                             predictor.vars = index.vars,
                             h = 2,
                             ncal = 30,
                             num.futures = 100,
                             window = 1000)
#> [1] "This is 1000"
#> [1] "This is 1001"
#> [1] "This is 1002"
#> [1] "This is 1003"
#> [1] "This is 1004"
#> [1] "This is 1005"
#> [1] "This is 1006"
#> [1] "This is 1007"
#> [1] "This is 1008"
#> [1] "This is 1009"
#> [1] "This is 1010"
#> [1] "This is 1011"
#> [1] "This is 1012"
#> [1] "This is 1013"
#> [1] "This is 1014"
#> [1] "This is 1015"
#> [1] "This is 1016"
#> [1] "This is 1017"
#> [1] "This is 1018"
#> [1] "This is 1019"
#> [1] "This is 1020"
#> [1] "This is 1021"
#> [1] "This is 1022"
#> [1] "This is 1023"
#> [1] "This is 1024"
#> [1] "This is 1025"
#> [1] "This is 1026"
#> [1] "This is 1027"
#> [1] "This is 1028"
#> [1] "This is 1029"
#> [1] "This is 1030"
#> [1] "This is 1031"
#> [1] "This is 1032"
#> [1] "This is 1033"
#> [1] "This is 1034"
#> [1] "This is 1035"
#> [1] "This is 1036"
#> [1] "This is 1037"
#> [1] "This is 1038"
#> [1] "This is 1039"
#> [1] "This is 1040"
#> [1] "This is 1041"
#> [1] "This is 1042"
#> [1] "This is 1043"
#> [1] "This is 1044"
#> [1] "This is 1045"
#> [1] "This is 1046"
#> [1] "This is 1047"
#> [1] "This is 1048"
#> [1] "This is 30"
#> [1] "This is 31"
#> [1] "This is 32"
#> [1] "This is 33"
#> [1] "This is 34"
#> [1] "This is 35"
#> [1] "This is 36"
#> [1] "This is 37"
#> [1] "This is 38"
#> [1] "This is 39"
#> [1] "This is 40"
#> [1] "This is 41"
#> [1] "This is 42"
#> [1] "This is 43"
#> [1] "This is 44"
#> [1] "This is 45"
#> [1] "This is 46"
#> [1] "This is 47"
#> [1] "This is 48"
                                 
# Mean width of generated 95% conformal bootstrap prediction intervals
width_data <- avgWidth(object = pprModel_cb)
width_data$mean
#>       h=1       h=2 
#> 0.4185444 0.3772966 
                                 
```
