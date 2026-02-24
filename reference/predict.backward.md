# Obtaining forecasts on a test set from a fitted `backward`

Gives forecasts on a test set.

## Usage

``` r
# S3 method for class 'backward'
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

  A `backward` object.

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

A `tsibble` with forecasts on test set.

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
predict(object = backwardModel, newdata = sim_test)
#> # A tsibble: 10 x 10 [1]
#> # Key:       dummy_key [1]
#>    inddd      y x_lag_000 x_lag_001 x_lag_002 x_lag_003 x_lag_004 x_lag_005
#>    <int>  <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
#>  1  1206 2.93      0.989     0.487     0.718     0.489     0.887     0.859 
#>  2  1207 0.841     0.0648    0.989     0.487     0.718     0.489     0.887 
#>  3  1208 0.0602    0.158     0.0648    0.989     0.487     0.718     0.489 
#>  4  1209 1.89      0.785     0.158     0.0648    0.989     0.487     0.718 
#>  5  1210 1.01      0.542     0.785     0.158     0.0648    0.989     0.487 
#>  6  1211 0.451     0.417     0.542     0.785     0.158     0.0648    0.989 
#>  7  1212 3.39      0.999     0.417     0.542     0.785     0.158     0.0648
#>  8  1213 1.29      0.256     0.999     0.417     0.542     0.785     0.158 
#>  9  1214 0.387     0.508     0.256     0.999     0.417     0.542     0.785 
#> 10  1215 0.542     0.0790    0.508     0.256     0.999     0.417     0.542 
#> # ℹ 2 more variables: dummy_key <dbl>, .predict <dbl>
```
