# Groupwise Additive Index Models (GAIM)

A wrapper for
[`cgaim::cgaim()`](https://rdrr.io/pkg/cgaim/man/cgaim.html) enabling
multiple GAIM models based on a grouping variable. Currently does not
support Constrained GAIM (CGAIM)s.

## Usage

``` r
model_gaim(
  data,
  yvar,
  neighbour = 0,
  index.vars,
  index.ind,
  s.vars = NULL,
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

- index.ind:

  An `integer` vector that assigns group index for each predictor in
  `index.vars`.

- s.vars:

  A `character` vector of names of the predictor variables for which
  splines should be fitted individually (rather than considering as part
  of an index).

- linear.vars:

  A `character` vector of names of the predictor variables that should
  be included linearly into the model.

- ...:

  Other arguments not currently used. (Note that the arguments in
  [`cgaim::cgaim()`](https://rdrr.io/pkg/cgaim/man/cgaim.html) related
  to constrained GAIMs are currently not supported. Furthermore, the
  argument `subset` is also not supported due to a bug in
  [`cgaim::cgaim()`](https://rdrr.io/pkg/cgaim/man/cgaim.html).)

## Value

An object of class `gaimFit`. This is a `tibble` with two columns:

- key:

  The level of the grouping variable (i.e. key of the training data
  set).

- fit:

  Information of the fitted model corresponding to the `key`.

Each row of the column `fit` is an object of class `cgaim`. For details
refer [`cgaim::cgaim()`](https://rdrr.io/pkg/cgaim/man/cgaim.html).

## Details

Group-wise Additive Index Model (GAIM) can be written in the form
\$\$y\_{i} = \sum\_{j = 1}^{p}
g\_{j}(\boldsymbol{\alpha}\_{j}^{T}\boldsymbol{x}\_{ij}) +
\varepsilon\_{i}, \quad i = 1, \dots, n,\$\$ where \\y\_{i}\\ is the
univariate response, \\\boldsymbol{x}\_{ij} \in \mathbb{R}^{l{j}}\\, \\j
= 1, \dots, p\\ are pre-specified non-overlapping subsets of
\\\boldsymbol{x}\_{i}\\, and \\\boldsymbol{\alpha}\_j\\ are the
corresponding index coefficients, \\g\_{j}\\ is an unknown (possibly
nonlinear) component function, and \\\varepsilon\_{i}\\ is the random
error, which is independent of \\\boldsymbol{x}\_{i}\\.

## See also

[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md),
[`model_backward`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md),
[`model_ppr`](https://nuwani-palihawadana.github.io/smimodel/reference/model_ppr.md),
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
    y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
    # Add an index to the data set
    inddd = seq(1, n)) |>
  drop_na() |>
  select(inddd, y1, starts_with("x_lag")) |>
  # Make the data set a `tsibble`
  as_tsibble(index = inddd)

# Predictors taken as index variables
index.vars <- colnames(sim_data)[3:7]

# Assign group indices for each predictor
index.ind = c(rep(1, 3), rep(2, 2))

# Predictors taken as non-linear variables not entering indices
s.vars = "x_lag_005"

# Model fitting
gaimModel <- model_gaim(data = sim_data,
                        yvar = "y1",
                        index.vars = index.vars,
                        index.ind = index.ind,
                        s.vars = s.vars)
#> [1] "model 1"
# Fitted model
gaimModel$fit[[1]]
#> Formula:
#> y1 ~ g(x_lag_000, x_lag_001, x_lag_002) + g(x_lag_003, x_lag_004) + 
#>     s(x_lag_005)
#> 
#> Coefficients:
#> (Intercept)   x_lag_000   x_lag_003   x_lag_005 
#>  1.23026577  0.99071595  0.41458553  0.03298014 
#> 
#> Indices weights:
#> x_lag_000 
#>   x_lag_000   x_lag_001   x_lag_002 
#> 0.604670638 0.392787126 0.002542237 
#> x_lag_003 
#>   x_lag_003   x_lag_004 
#>  0.98152724 -0.01847276 
#> 
#> Residual sum of squares: 0.06248322
```
