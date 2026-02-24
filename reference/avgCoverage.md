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

## Examples

``` r
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1105
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
sim_test <- sim_data[1001:1100, ]

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Model fitting
pprModel <- model_ppr(data = sim_train,
                      yvar = "y",
                      index.vars = index.vars)
#> [1] "model 1"
                      
# Block bootstrap prediction intervals (3-steps-ahead interval forecasts)
set.seed(12345)
pprModel_bb <- bb_cvforecast(object = pprModel,
                             data = sim_data,
                             yvar = "y",
                             predictor.vars = index.vars,
                             h = 3,
                             num.futures = 50,
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
#> [1] "This is 1049"
#> [1] "This is 1050"
#> [1] "This is 1051"
#> [1] "This is 1052"
#> [1] "This is 1053"
#> [1] "This is 1054"
#> [1] "This is 1055"
#> [1] "This is 1056"
#> [1] "This is 1057"
#> [1] "This is 1058"
#> [1] "This is 1059"
#> [1] "This is 1060"
#> [1] "This is 1061"
#> [1] "This is 1062"
#> [1] "This is 1063"
#> [1] "This is 1064"
#> [1] "This is 1065"
#> [1] "This is 1066"
#> [1] "This is 1067"
#> [1] "This is 1068"
#> [1] "This is 1069"
#> [1] "This is 1070"
#> [1] "This is 1071"
#> [1] "This is 1072"
#> [1] "This is 1073"
#> [1] "This is 1074"
#> [1] "This is 1075"
#> [1] "This is 1076"
#> [1] "This is 1077"
#> [1] "This is 1078"
#> [1] "This is 1079"
#> [1] "This is 1080"
#> [1] "This is 1081"
#> [1] "This is 1082"
#> [1] "This is 1083"
#> [1] "This is 1084"
#> [1] "This is 1085"
#> [1] "This is 1086"
#> [1] "This is 1087"
#> [1] "This is 1088"
#> [1] "This is 1089"
#> [1] "This is 1090"
#> [1] "This is 1091"
#> [1] "This is 1092"
#> [1] "This is 1093"
#> [1] "This is 1094"
#> [1] "This is 1095"
#> [1] "This is 1096"
#> [1] "This is 1097"
                             
# Mean coverage of generated 95% block bootstrap prediction intervals
cov_data <- avgCoverage(object = pprModel_bb)
cov_data$mean
#>       h=1       h=2       h=3 
#> 0.8775510 0.9285714 0.9285714 
                                 
```
