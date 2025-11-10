# Sparse Multiple Index (SMI) Models

Fits nonparametric multiple index model(s), with simultaneous predictor
selection (hence "sparse") and predictor grouping. Possible to fit
multiple SMI models based on a grouping variable.

## Usage

``` r
model_smimodel(
  data,
  yvar,
  neighbour = 0,
  family = gaussian(),
  index.vars,
  initialise = c("ppr", "additive", "linear", "multiple", "userInput"),
  num_ind = 5,
  num_models = 5,
  seed = 123,
  index.ind = NULL,
  index.coefs = NULL,
  s.vars = NULL,
  linear.vars = NULL,
  lambda0 = 1,
  lambda2 = 1,
  M = 10,
  max.iter = 50,
  tol = 0.001,
  tolCoefs = 0.001,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE
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

- family:

  A description of the error distribution and link function to be used
  in the model (see [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`family`](https://rdrr.io/r/stats/family.html)).

- index.vars:

  A `character` vector of names of the predictor variables for which
  indices should be estimated.

- initialise:

  The model structure with which the estimation process should be
  initialised. The default is `"ppr"`, where the initial model is
  derived from projection pursuit regression. The other options are
  `"additive"` - nonparametric additive model, `"linear"` - linear
  regression model (i.e. a special case single-index model, where the
  initial values of the index coefficients are obtained through a linear
  regression), `"multiple"` - multiple models are fitted starting with
  different initial models (number of indices = `num_ind`; `num_models`
  random instances of the model (i.e. the predictor assignment to
  indices and initial index coefficients are generated randomly) are
  considered), and the final optimal model with the lowest loss is
  returned, and `"userInput"` - user specifies the initial model
  structure (i.e. the number of indices and the placement of index
  variables among indices) and the initial index coefficients through
  `index.ind` and `index.coefs` arguments respectively.

- num_ind:

  If `initialise = "ppr"` or `"multiple"`: an `integer` that specifies
  the number of indices to be used in the model(s). The default is
  `num_ind = 5`.

- num_models:

  If `initialise = "multiple"`: an `integer` that specifies the number
  of starting models to be checked. The default is `num_models = 5`.

- seed:

  If `initialise = "multiple"`: the seed to be set when generating
  random starting points.

- index.ind:

  If `initialise = "userInput"`: an `integer` vector that assigns group
  index for each predictor in `index.vars`.

- index.coefs:

  If `initialise = "userInput"`: a `numeric` vector of index
  coefficients.

- s.vars:

  A `character` vector of names of the predictor variables for which
  splines should be fitted individually (rather than considering as part
  of an index).

- linear.vars:

  A `character` vector of names of the predictor variables that should
  be included linearly into the model.

- lambda0:

  Penalty parameter for L0 penalty.

- lambda2:

  Penalty parameter for L2 penalty.

- M:

  Big-M value to be used in MIP.

- max.iter:

  Maximum number of MIP iterations performed to update index
  coefficients for a given model.

- tol:

  Tolerance for the objective function value (loss) of MIP.

- tolCoefs:

  Tolerance for coefficients.

- TimeLimit:

  A limit for the total time (in seconds) expended in a single MIP
  iteration.

- MIPGap:

  Relative MIP optimality gap.

- NonConvex:

  The strategy for handling non-convex quadratic objectives or
  non-convex quadratic constraints in Gurobi solver.

- verbose:

  The option to print detailed solver output.

## Value

An object of class `smimodel`. This is a `tibble` with two columns:

- key:

  The level of the grouping variable (i.e. key of the training data
  set).

- fit:

  Information of the fitted model corresponding to the `key`.

Each row of the column `fit` contains a list with two elements:

- initial:

  A list of information of the model initialisation. (For descriptions
  of the list elements see
  [`make_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)).

- best:

  A list of information of the final optimised model. (For descriptions
  of the list elements see
  [`make_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)).

## Details

Sparse Multiple Index (SMI) model is a semi-parametric model that can be
written as \$\$y\_{i} = \beta\_{0} + \sum\_{j =
1}^{p}g\_{j}(\boldsymbol{\alpha}\_{j}^{T}\boldsymbol{x}\_{ij}) +
\sum\_{k = 1}^{d}f\_{k}(w\_{ik}) +
\boldsymbol{\theta}^{T}\boldsymbol{u}\_{i} + \varepsilon\_{i}, \quad i =
1, \dots, n,\$\$ where \\y\_{i}\\ is the univariate response,
\\\beta\_{0}\\ is the model intercept, \\\boldsymbol{x}\_{ij} \in
\mathbb{R}^{l\_{j}}\\, \\j = 1, \dots, p\\ are \\p\\ subsets of
predictors entering indices, \\\boldsymbol{\alpha}\_{j}\\ is a vector of
index coefficients corresponding to the index \\h\_{ij} =
\boldsymbol{\alpha}\_{j}^{T}\boldsymbol{x}\_{ij}\\, and \\g\_{j}\\ is a
smooth nonlinear function (estimated by a penalised cubic regression
spline). The model also allows for predictors that do not enter any
indices, including covariates \\w\_{ik}\\ that relate to the response
through nonlinear functions \\f\_{k}\\, \\k = 1, \dots, d\\, and linear
covariates \\\boldsymbol{u}\_{i}\\.

In the model formulation related to this implementation, both the number
of indices \\p\\ and the predictor grouping among indices are assumed to
be unknown prior to model estimation. Suppose we observe
\\y_1,\dots,y_n\\, along with a set of potential predictors,
\\\boldsymbol{x}\_1,\dots,\boldsymbol{x}\_n\\, with each vector
\\\boldsymbol{x}\_i\\ containing \\q\\ predictors. This function
implements algorithmic variable selection for index variables (i.e.
predictors entering indices) of the SMI model by allowing for zero index
coefficients for predictors. Non-overlapping predictors among indices
are assumed (i.e. no predictor enters more than one index). For
algorithmic details see reference.

## References

Palihawadana, N.K., Hyndman, R.J. & Wang, X. (2024). Sparse Multiple
Index Models for High-Dimensional Nonparametric Forecasting.
<https://www.monash.edu/business/ebs/research/publications/ebs/2024/wp16-2024.pdf>.

## See also

[`greedy_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/greedy_smimodel.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
library(ROI)
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

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Model fitting
smimodel_ppr <- model_smimodel(data = sim_data,
                               yvar = "y1",
                               index.vars = index.vars,
                               initialise = "ppr")

# Best (optimised) fitted model
smimodel_ppr$fit[[1]]$best
} # }
```
