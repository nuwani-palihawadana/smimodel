# Projection Pursuit Regression (PPR) models

A wrapper for [`stats::ppr()`](https://rdrr.io/r/stats/ppr.html)
enabling multiple PPR models based on a grouping variable.

## Usage

``` r
model_ppr(data, yvar, neighbour = 0, index.vars, num_ind = 5, ...)
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

- index.vars:

  A `character` vector of names of the predictor variables for which
  indices should be estimated.

- num_ind:

  An `integer` that specifies the number of indices to be used in the
  model(s). (Corresponds to `nterms` in
  [`stats::ppr()`](https://rdrr.io/r/stats/ppr.html).)

- ...:

  Other arguments not currently used. (For more information on other
  arguments that can be passed, refer
  [`stats::ppr()`](https://rdrr.io/r/stats/ppr.html).)

## Value

An object of class `pprFit`. This is a `tibble` with two columns:

- key:

  The level of the grouping variable (i.e. key of the training data
  set).

- fit:

  Information of the fitted model corresponding to the `key`.

Each row of the column `fit` is an object of class
`c("ppr.form", "ppr")`. For details refer
[`stats::ppr()`](https://rdrr.io/r/stats/ppr.html).

## Details

A Projection Pursuit Regression (PPR) model (Friedman & Stuetzle (1981))
is given by \$\$y\_{i} = \sum\_{j=1}^{p}
{g\_{j}(\boldsymbol{\alpha}\_{j}^{T}\boldsymbol{x}\_{i})} +
\varepsilon\_{i}, \quad i = 1, \dots, n,\$\$ where \\y\_{i}\\ is the
response, \\\boldsymbol{x}\_{i}\\ is the \\q\\-dimensional predictor
vector, \\\boldsymbol{\alpha}\_{j} = ( \alpha\_{j1}, \dots, \alpha\_{jp}
)^{T}\\, \\j = 1, \dots, p\\ are \\q\\-dimensional projection vectors
(or vectors of "index coefficients"), \\g\_{j}\\'s are unknown nonlinear
functions, and \\\varepsilon\_{i}\\ is the random error.

## References

Friedman, J. H. & Stuetzle, W. (1981). Projection pursuit regression.
*Journal of the American Statistical Association*, 76, 817–823.
[doi:10.2307/2287576](https://www.tandfonline.com/doi/abs/10.1080/01621459.1981.10477729).

## See also

[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md),
[`model_backward`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md),
[`model_gaim`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gaim.md),
[`model_gam`](https://nuwani-palihawadana.github.io/smimodel/reference/model_gam.md),
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

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Model fitting
pprModel <- model_ppr(data = sim_data,
                      yvar = "y",
                      index.vars = index.vars)
#> [1] "model 1"

# Fitted model
pprModel$fit[[1]]
#> Call:
#> ppr(formula = as.formula(pre.formula), data = df_cat, nterms = num_ind)
#> 
#> Goodness of fit:
#>  5 terms 
#> 9.210028 
```
