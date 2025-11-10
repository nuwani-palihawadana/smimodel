# SMI model with a given penalty parameter combination

Fits a nonparametric multiple index model to the data for a given
combination of the penalty parameters (lambda0, lambda2), and returns
the validation set mean squared error (MSE). (Used within
[`greedy.fit`](https://nuwani-palihawadana.github.io/smimodel/reference/greedy.fit.md);
users are not expected to use this function directly.)

## Usage

``` r
tune_smimodel(
  data,
  val.data,
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
  lambda.comb = c(1, 1),
  M = 10,
  max.iter = 50,
  tol = 0.001,
  tolCoefs = 0.001,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL
)
```

## Arguments

- data:

  Training data set on which models will be trained. Must be a data set
  of class `tsibble`.(Make sure there are no additional date or time
  related variables except for the `index` of the `tsibble`).

- val.data:

  Validation data set. (The data set on which the penalty parameter
  selection will be performed.) Must be a data set of class `tsibble`.
  (Once the penalty parameter selection is completed, the best model
  will be re-fitted for the combined data set `data + val.data`.)

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

- lambda.comb:

  A `numeric` vector (of length two) indicating the values for the two
  penalty parameters lambda0 and lambda2.

- M:

  Big-M value used in MIP.

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

A `numeric`.
