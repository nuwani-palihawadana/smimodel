# Updating index coefficients using MIP

Updates index coefficients by solving a mixed integer program.

## Usage

``` r
update_alpha(
  Y,
  X,
  num_pred,
  num_ind,
  index.ind,
  dgz,
  alpha_old,
  lambda0 = 1,
  lambda2 = 1,
  M = 10,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE
)
```

## Arguments

- Y:

  Column matrix of response.

- X:

  Matrix of predictors (size adjusted to number of indices).

- num_pred:

  Number of predictors.

- num_ind:

  Number of indices.

- index.ind:

  An integer vector that assigns group index for each predictor.

- dgz:

  The `tibble` of derivatives of the estimated smooths from previous
  iteration.

- alpha_old:

  Vector of index coefficients from previous iteration.

- lambda0:

  Penalty parameter for L0 penalty.

- lambda2:

  Penalty parameter for L2 penalty.

- M:

  Big-M value to be used in MIP.

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

A vector of normalised index coefficients.
