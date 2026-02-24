# Obtaining forecasts on a test set from a fitted `gamFit`

Gives forecasts on a test set.

## Usage

``` r
# S3 method for class 'gamFit'
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

  A `gamFit` object.

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

predict(object = gamModel, newdata = sim_test)
#> # A tsibble: 10 x 10 [1]
#> # Key:       dummy_key [1]
#>    inddd       y x_lag_000 x_lag_001 x_lag_002 x_lag_003 x_lag_004 x_lag_005
#>    <int>   <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
#>  1  1006  1.13      0.478     0.848     0.853     0.160     0.594     0.274 
#>  2  1007  2.72      0.774     0.478     0.848     0.853     0.160     0.594 
#>  3  1008  1.54      0.295     0.774     0.478     0.848     0.853     0.160 
#>  4  1009  0.0723    0.0656    0.295     0.774     0.478     0.848     0.853 
#>  5  1010  0.401     0.441     0.0656    0.295     0.774     0.478     0.848 
#>  6  1011  0.524     0.462     0.441     0.0656    0.295     0.774     0.478 
#>  7  1012  0.0927    0.341     0.462     0.441     0.0656    0.295     0.774 
#>  8  1013  0.310     0.185     0.341     0.462     0.441     0.0656    0.295 
#>  9  1014  0.279     0.507     0.185     0.341     0.462     0.441     0.0656
#> 10  1015 -0.0138    0.0192    0.507     0.185     0.341     0.462     0.441 
#> # ℹ 2 more variables: dummy_key <dbl>, .predict <dbl>
```
