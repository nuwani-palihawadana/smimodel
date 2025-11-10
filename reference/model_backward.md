# Nonparametric Additive Model with Backward Elimination

Fits a nonparametric additive model, with simultaneous variable
selection through a backward elimination procedure as proposed by Fan
and Hyndman (2012).

## Usage

``` r
model_backward(
  data,
  val.data,
  yvar,
  neighbour = 0,
  family = gaussian(),
  s.vars = NULL,
  s.basedim = NULL,
  linear.vars = NULL,
  refit = TRUE,
  tol = 0.001,
  parallel = FALSE,
  workers = NULL,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL
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

- val.data:

  Validation data set. (The data set on which the model selection will
  be performed.) Must be a data set of class `tsibble`.

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

- family:

  A description of the error distribution and link function to be used
  in the model (see [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`family`](https://rdrr.io/r/stats/family.html)).

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

- refit:

  Whether to refit the model combining training and validation sets
  after model selection. If `FALSE`, the final model will be estimated
  only on the training set.

- tol:

  Tolerance for the ratio of relative change in validation set MSE, used
  in model selection.

- parallel:

  Whether to use parallel computing in model selection or not.

- workers:

  If `parallel = TRUE`, number of workers to use.

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string. (Since the nonlinear
  functions are estimated using splines, extrapolation is not desirable.
  Hence, if any predictor variable in `val.data` that is treated
  non-linearly in the estimated model, will be truncated to be in the
  in-sample range before obtaining predictions. If any variables are
  listed here will be excluded from such truncation.)

- recursive:

  Whether to obtain recursive forecasts or not (default - `FALSE`).

- recursive_colRange:

  If `recursive = TRUE`, the range of column numbers in `val.data` to be
  filled with forecasts. Recursive/autoregressive forecasting is
  required when the lags of the response variable itself are used as
  predictor variables into the model. Make sure such lagged variables
  are positioned together in increasing lag order (i.e.
  `lag_1, lag_2, ..., lag_m`, `lag_m =` maximum lag used) in `val.data`,
  with no break in the lagged variable sequence even if some of the
  intermediate lags are not used as predictors.

## Value

An object of class `backward`. This is a `tibble` with two columns:

- key:

  The level of the grouping variable (i.e. key of the training data
  set).

- fit:

  Information of the fitted model corresponding to the `key`.

Each row of the column `fit` is an object of class `gam`. For details
refer [`mgcv::gamObject`](https://rdrr.io/pkg/mgcv/man/gamObject.html).

## Details

This function fits a nonparametric additive model formulated through
Backward Elimination, as proposed by Fan and Hyndman (2012). The process
starts with all predictors included in an additive model, and predictors
are progressively omitted until the best model is obtained based on the
validation set. Once the best model is obtained, the final model is
re-fitted for the data set combining training and validation sets. For
more details see reference.

## References

Fan, S. & Hyndman, R.J. (2012). Short-Term Load Forecasting Based on a
Semi-Parametric Additive Model. *IEEE Transactions on Power Systems*,
27(1), 134-141. <http://doi.org/10.1109/TPWRS.2011.2162082>

## See also

[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md),
[`model_gaim`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gaim.md),
[`model_ppr`](https://nuwani-palihawadana.github.io/smimodel/reference/model_ppr.md),
[`model_gam`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gam.md),
[`model_lm`](https://nuwani-palihawadana.github.io/smimodel/reference/model_lm.md)

## Examples

``` r
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)
#> 
#> Attaching package: ‘tsibble’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, union

# Simulate data
n = 1205
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(n)) |>
  mutate(
    # Add x_lags
    x_lag = lag_matrix(x_lag_000, 5)) |>
  unpack(x_lag, names_sep = "_") |>
  mutate(
    # Response variable
    y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
    # Add an index to the data set
    inddd = seq(1, n)) |>
  drop_na() |>
  select(inddd, y1, starts_with("x_lag")) |>
  # Make the data set a `tsibble`
  as_tsibble(index = inddd)

# Training set
sim_train <- sim_data[1:1000, ]
# Validation set
sim_val <- sim_data[1001:1200, ]

# Predictors taken as non-linear variables
s.vars <- colnames(sim_data)[3:8]

# Model fitting
backwardModel <- model_backward(data = sim_train,
                                val.data = sim_val,
                                yvar = "y1",
                                s.vars = s.vars)
#> [1] "Model 1 fitted!"
# Fitted model
backwardModel$fit[[1]]
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> y1 ~ +s(x_lag_000, bs = "cr") + s(x_lag_001, bs = "cr") + s(x_lag_002, 
#>     bs = "cr") + s(x_lag_003, bs = "cr") + s(x_lag_005, bs = "cr")
#> 
#> Estimated degrees of freedom:
#> 4.80 3.51 1.00 3.28 1.00  total = 14.59 
#> 
#> REML score: 480.794     
```
