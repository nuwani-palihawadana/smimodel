# Updating index coefficients and non-linear functions iteratively

Iteratively updates index coefficients and non-linear functions using
mixed integer programming. (A helper function used within
[`update_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/update_smimodelFit.md);
users are not expected to directly call this function.)

## Usage

``` r
inner_update(
  x,
  data,
  yvar,
  family = gaussian(),
  index.vars,
  s.vars,
  linear.vars,
  num_ind,
  dgz,
  alpha_old,
  lambda0 = 1,
  lambda2 = 1,
  M = 10,
  max.iter = 50,
  tol = 0.001,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE
)
```

## Arguments

- x:

  Fitted `gam`.

- data:

  Training data set on which models will be trained. Should be a
  `tsibble`.

- yvar:

  Name of the response variable as a character string.

- family:

  A description of the error distribution and link function to be used
  in the model (see [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`family`](https://rdrr.io/r/stats/family.html)).

- index.vars:

  A `character` vector of names of the predictor variables for which
  indices should be estimated.

- s.vars:

  A `character` vector of names of the predictor variables for which
  splines should be fitted individually (rather than considering as part
  of an index).

- linear.vars:

  A `character` vector of names of the predictor variables that should
  be included linearly into the model.

- num_ind:

  Number of indices.

- dgz:

  The `tibble` of derivatives of the estimated smooths.

- alpha_old:

  Current vector of index coefficients.

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

  Tolerance for loss.

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

A list containing following elements:

- best_alpha:

  The vector of best index coefficient estimates.

- min_loss:

  Minimum value of the objective function(loss).

- index.ind:

  An `integer` vector that assigns group index for each predictor,
  corresponding to `best_alpha`.

- ind_pos:

  A list that indicates which predictors belong to which index,
  corresponding to `best_alpha`.

- X_new:

  A matrix of selected predictor variables, corresponding to
  `best_alpha`.
