# Updating a `smimodelFit`

Optimises and updates a given `smimodelFit`.

## Usage

``` r
update_smimodelFit(
  object,
  data,
  lambda0 = 1,
  lambda2 = 1,
  M = 10,
  max.iter = 50,
  tol = 0.001,
  tolCoefs = 0.001,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  A `smimodelFit` object.

- data:

  Training data set on which models will be trained. Must be a data set
  of class `tsibble`.(Make sure there are no additional date or time
  related variables except for the `index` of the `tsibble`).

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

- ...:

  Other arguments not currently used.

## Value

A list of optimised model information. For descriptions of the list
elements see
[`make_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)).
