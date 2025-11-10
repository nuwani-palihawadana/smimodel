# Linear Regression models

A wrapper for [`lm`](https://rdrr.io/r/stats/lm.html) enabling multiple
linear models based on a grouping variable.

## Usage

``` r
model_lm(data, yvar, neighbour = 0, linear.vars, ...)
```

## Arguments

- data:

  Training data set on which models will be trained. Must be a data set
  of class `tsibble`.(Make sure there are no additional date or time
  related variables except for the `index` of the `tsibble`). If
  multiple models are fitted, the grouping variable should be the `key`
  of the `tsibble`. If a `key` is not specified, a dummy key with only
  one level will be created.

- yvar:

  Name of the response variable as a character string.

- neighbour:

  If multiple models are fitted: Number of neighbours of each key (i.e.
  grouping variable) to be considered in model fitting to handle
  smoothing over the key. Should be an `integer`. If `neighbour = x`,
  `x` number of keys before the key of interest and `x` number of keys
  after the key of interest are grouped together for model fitting. The
  default is `neighbour = 0` (i.e. no neighbours are considered for
  model fitting).

- linear.vars:

  A character vector of names of the predictor variables.

- ...:

  Other arguments not currently used.

## Value

An object of class `lmFit`. This is a `tibble` with two columns:

- key:

  The level of the grouping variable (i.e. key of the training data
  set).

- fit:

  Information of the fitted model corresponding to the `key`.

Each row of the column `fit` is an object of class `lm`. For details
refer [`stats::lm`](https://rdrr.io/r/stats/lm.html).

## See also

[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md),
[`model_backward`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md),
[`model_gaim`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gaim.md),
[`model_ppr`](https://nuwani-palihawadana.github.io/smimodel/reference/model_ppr.md),
[`model_gam`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gam.md)

## Examples

``` r
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1005
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

# Predictor variables
linear.vars <- colnames(sim_data)[3:8]

# Model fitting
lmModel <- model_lm(data = sim_data,
                    yvar = "y",
                    linear.vars = linear.vars)
#> [1] "model 1"
# Fitted model
lmModel$fit[[1]]
#> 
#> Call:
#> stats::lm(formula = as.formula(pre.formula), data = df_cat)
#> 
#> Coefficients:
#> (Intercept)    x_lag_000    x_lag_001    x_lag_002    x_lag_003    x_lag_004  
#>   -1.793753     2.770954     1.821332     0.048464     1.414372     0.038152  
#>   x_lag_005  
#>   -0.006478  
#> 
```
