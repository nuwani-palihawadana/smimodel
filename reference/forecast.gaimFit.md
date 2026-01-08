# Forecasting using GAIMs

Returns forecasts and other information for GAIMs.

## Usage

``` r
# S3 method for class 'gaimFit'
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

  An object of class `gaimFit`. Usually the result of a call to
  [`model_gaim`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gaim.md).

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

# Predictors taken as index variables
index.vars <- colnames(sim_data)[3:7]

# Assign group indices for each predictor
index.ind = c(rep(1, 3), rep(2, 2))

# Predictors taken as non-linear variables not entering indices
s.vars = "x_lag_005"

# Model fitting
gaimModel <- model_gaim(data = sim_train,
                        yvar = "y",
                        index.vars = index.vars,
                        index.ind = index.ind,
                        s.vars = s.vars)
#> [1] "model 1"
                        
forecast(gaimModel, newdata = sim_test)
#>    Point Forecast
#> 1      1.12501980
#> 2      2.35295038
#> 3      1.43434126
#> 4      0.03578922
#> 5      0.65434928
#> 6      0.49947660
#> 7     -0.01380044
#> 8      0.08583838
#> 9      0.43564545
#> 10    -0.01533303
```
