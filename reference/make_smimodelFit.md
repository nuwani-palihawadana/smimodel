# Converting a fitted `gam` object to a `smimodelFit` object

Converts a given object of class `gam` to an object of class
`smimodelFit`.

## Usage

``` r
make_smimodelFit(
  x,
  data,
  yvar,
  neighbour,
  index.vars,
  index.ind,
  index.data,
  index.names,
  alpha,
  s.vars = NULL,
  linear.vars = NULL,
  lambda0 = NULL,
  lambda2 = NULL,
  M = NULL,
  max.iter = NULL,
  tol = NULL,
  tolCoefs = NULL,
  TimeLimit = NULL,
  MIPGap = NULL,
  NonConvex = NULL
)
```

## Arguments

- x:

  A fitted `gam` object.

- data:

  The original training data set.

- yvar:

  Name of the response variable as a character string.

- neighbour:

  `neighbour` argument passed from the outer function.

- index.vars:

  A `character` vector of names of the predictor variables for which
  indices are estimated.

- index.ind:

  An `integer` vector that assigns group index for each predictor in
  `index.vars`.

- index.data:

  A `tibble` including columns for the constructed indices.

- index.names:

  A `character` vector of names of the constructed indices.

- alpha:

  A vector of index coefficients.

- s.vars:

  A `character` vector of names of the predictor variables for which
  splines are fitted individually (rather than considering as part of an
  index).

- linear.vars:

  A `character` vector of names of the predictor variables that are
  included linearly in the model.

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

## Value

An object of class `smimodelFit`, which is a list that contains
following elements:

- alpha:

  A sparse matrix of index coefficients vectors. Each column of the
  matrix corresponds to the index coefficient vector of each index.

- derivatives:

  A `tibble` of derivatives of the estimated smooths.

- var_y:

  Name of the response variable.

- vars_index:

  A `character` vector of names of the predictor variables for which
  indices are estimated.

- vars_s:

  A `character` vector of names of the predictor variables for which
  splines are fitted individually.

- vars_linear:

  A `character` vector of names of the predictor variables that are
  included linearly in the model.

- neighbour:

  Number of neighbours of each key considered in model fitting.

- gam:

  Fitted `gam`.

- lambda0:

  L0 penalty parameter used for model fitting.

- lambda2:

  L2 penalty parameter used for model fitting.

- M:

  Big-M value used in MIP.

- max.iter:

  Maximum number of MIP iterations for a single round of index
  coefficients update.

- tol:

  Tolerance for the objective function value (loss) used in solving MIP.

- tolCoefs:

  Tolerance for coefficients used in updating index coefficients.

- TimeLimit:

  Limit for the total time (in seconds) expended in a single MIP
  iteration.

- MIPGap:

  Relative MIP optimality gap used.

- Nonconvex:

  The strategy used for handling non-convex quadratic objectives or
  non-convex quadratic constraints in Gurobi solver.
