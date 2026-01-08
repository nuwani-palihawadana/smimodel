# Forecasting using GAMs

Returns forecasts and other information for GAMs.

## Usage

``` r
# S3 method for class 'gamFit'
forecast(
  object,
  h = 1,
  level = c(80, 95),
  newdata,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `gamFit`. Usually the result of a call to
  [`model_gam`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gam.md).

- h:

  Forecast horizon.

- level:

  Confidence level for prediction intervals.

- newdata:

  The set of new data on for which the forecasts are required (i.e. test
  set; should be a `tsibble`).

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string.

- recursive:

  Whether to obtain recursive forecasts or not (default - `FALSE`).

- recursive_colRange:

  If `recursive = TRUE`, the range of column numbers in `newdata` to be
  filled with forecasts.

- ...:

  Other arguments not currently used.

## Value

An object of class `forecast`. Here, it is a list containing the
following elements:

- method:

  The name of the forecasting method as a character string.

- model:

  The fitted model.

- mean:

  Point forecasts as a time series.

- residuals:

  Residuals from the fitted model.

- fitted:

  Fitted values (one-step forecasts).

## Examples

``` r
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1015
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
sim_test <- sim_data[1001:1010, ]

# Predictors taken as non-linear variables
s.vars <- colnames(sim_data)[3:6]

# Predictors taken as linear variables
linear.vars <- colnames(sim_data)[7:8]

# Model fitting
gamModel <- model_gam(data = sim_train,
                      yvar = "y",
                      s.vars = s.vars,
                      linear.vars = linear.vars)
#> [1] "model 1"

forecast(gamModel, newdata = sim_test)
#>    Point Forecast
#> 1     1.245005200
#> 2     2.378001533
#> 3     1.659849053
#> 4    -0.228756505
#> 5     0.602818141
#> 6     0.434930625
#> 7     0.026903622
#> 8    -0.002923976
#> 9     0.463371054
#> 10   -0.161817694
```
