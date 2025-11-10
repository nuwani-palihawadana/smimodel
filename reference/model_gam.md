# Generalised Additive Models

A wrapper for [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html)
enabling multiple GAMs based on a grouping variable.

## Usage

``` r
model_gam(
  data,
  yvar,
  family = gaussian(),
  neighbour = 0,
  s.vars,
  s.basedim = NULL,
  linear.vars = NULL,
  ...
)
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

- family:

  A description of the error distribution and link function to be used
  in the model (see [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`family`](https://rdrr.io/r/stats/family.html)).

- neighbour:

  If multiple models are fitted: Number of neighbours of each key (i.e.
  grouping variable) to be considered in model fitting to handle
  smoothing over the key. Should be an `integer`. If `neighbour = x`,
  `x` number of keys before the key of interest and `x` number of keys
  after the key of interest are grouped together for model fitting. The
  default is `neighbour = 0` (i.e. no neighbours are considered for
  model fitting).

- s.vars:

  A `character` vector of names of the predictor variables for which
  splines should be fitted (i.e. non-linear predictors).

- s.basedim:

  Dimension of the bases used to represent the smooth terms
  corresponding to `s.vars`. (For more information refer
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html).)

- linear.vars:

  A `character` vector of names of the predictor variables that should
  be included linearly into the model (i.e. linear predictors).

- ...:

  Other arguments not currently used.

## Value

An object of class `gamFit`. This is a `tibble` with two columns:

- key:

  The level of the grouping variable (i.e. key of the training data
  set).

- fit:

  Information of the fitted model corresponding to the `key`.

Each row of the column `fit` is an object of class `gam`. For details
refer [`mgcv::gamObject`](https://rdrr.io/pkg/mgcv/man/gamObject.html).

## See also

[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md),
[`model_backward`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md),
[`model_gaim`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gaim.md),
[`model_ppr`](https://nuwani-palihawadana.github.io/smimodel/reference/model_ppr.md),
[`model_lm`](https://nuwani-palihawadana.github.io/smimodel/reference/model_lm.md)

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

# Predictors taken as non-linear variables
s.vars <- colnames(sim_data)[3:6]

# Predictors taken as linear variables
linear.vars <- colnames(sim_data)[7:8]

# Model fitting
gamModel <- model_gam(data = sim_data,
                      yvar = "y",
                      s.vars = s.vars,
                      linear.vars = linear.vars)
#> [1] "model 1"

# Fitted model
gamModel$fit[[1]]
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> y ~ s(x_lag_000, bs = "cr") + s(x_lag_001, bs = "cr") + s(x_lag_002, 
#>     bs = "cr") + s(x_lag_003, bs = "cr") + x_lag_004 + x_lag_005
#> 
#> Estimated degrees of freedom:
#> 4.62 3.51 1.39 2.41  total = 14.93 
#> 
#> REML score: 397.9157     
```
