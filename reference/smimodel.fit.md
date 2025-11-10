# SMI model estimation

Fits a single nonparametric multiple index model to the data. This is a
helper function designed to be called from user-facing wrapper
functions,
[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md)
and
[`greedy_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/greedy_smimodel.md).

## Usage

``` r
smimodel.fit(
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
  related variables except for the `index` of the `tsibble`).

- yvar:

  Name of the response variable as a character string.

- neighbour:

  `neighbour` argument passed from the outer function.

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

A list with two elements:

- initial:

  A list of information of the model initialisation. (For descriptions
  of the list elements see
  [`make_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)).

- best:

  A list of information of the final optimised model. (For descriptions
  of the list elements see
  [`make_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)).
