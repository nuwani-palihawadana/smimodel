# Forecasting using nonparametric additive models with backward elimination

Returns forecasts and other information for nonparametric additive
models with backward elimination.

## Usage

``` r
# S3 method for class 'backward'
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

  An object of class `backward`. Usually the result of a call to
  [`model_backward`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md).

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
n = 1215
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
# Validation set
sim_val <- sim_data[1001:1200, ]
# Test set
sim_test <- sim_data[1201:1210, ]

# Predictors taken as non-linear variables
s.vars <- colnames(sim_data)[3:8]

# Model fitting
backwardModel <- model_backward(data = sim_train,
                                val.data = sim_val,
                                yvar = "y",
                                s.vars = s.vars)
#> [1] "Model 1 fitted!"
forecast(backwardModel, newdata = sim_test)
#>    Point Forecast
#> 1       2.8363732
#> 2       1.5316512
#> 3      -0.3403322
#> 4       2.1579572
#> 5       1.1067827
#> 6       0.4067579
#> 7       3.2047193
#> 8       1.5858674
#> 9       0.5086024
#> 10      0.9975114
```
