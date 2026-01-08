# Forecasting using SMI models

Returns forecasts and other information for SMI models.

## Usage

``` r
# S3 method for class 'smimodel'
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

  An object of class `smimodel`. Usually the result of a call to
  [`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md).

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
if (FALSE) { # \dontrun{
library(dplyr)
library(ROI)
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

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Model fitting
smimodel_ppr <- model_smimodel(data = sim_train,
                               yvar = "y",
                               index.vars = index.vars,
                               initialise = "ppr")

forecast(smimodel_ppr, newdata = sim_test)
} # }
```
