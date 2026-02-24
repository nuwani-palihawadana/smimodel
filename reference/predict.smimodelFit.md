# Obtaining forecasts on a test set from a `smimodelFit`

Gives forecasts on a test set.

## Usage

``` r
# S3 method for class 'smimodelFit'
predict(
  object,
  newdata,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL,
  ...
)
```

## Arguments

- object:

  A `smimodelFit` object.

- newdata:

  The set of new data on for which the forecasts are required (i.e. test
  set; should be a `tsibble`).

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string. (Since the nonlinear
  functions are estimated using splines, extrapolation is not desirable.
  Hence, if any predictor variable in `newdata` that is treated
  non-linearly in the estimated model, will be truncated to be in the
  in-sample range before obtaining predictions. If any variables are
  listed here will be excluded from such truncation.)

- recursive:

  Whether to obtain recursive forecasts or not (default - `FALSE`).

- recursive_colRange:

  If `recursive = TRUE`, the range of column numbers in `newdata` to be
  filled with forecasts. Recursive/autoregressive forecasting is
  required when the lags of the response variable itself are used as
  predictor variables into the model. Make sure such lagged variables
  are positioned together in increasing lag order (i.e.
  `lag_1, lag_2, ..., lag_m`, `lag_m =` maximum lag used) in `newdata`,
  with no break in the lagged variable sequence even if some of the
  intermediate lags are not used as predictors.

- ...:

  Other arguments not currently used.

## Value

A `tibble` with forecasts on test set.

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

predict(object = smimodel_ppr, newdata = sim_test)
} # }
```
