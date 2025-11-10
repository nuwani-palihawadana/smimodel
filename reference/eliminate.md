# Eliminate a variable and fit a nonparametric additive model

Eliminates a specified variable and fits a nonparametric additive model
with remaining variables, and returns validation set MSE. This is an
internal function of the package, and designed to be called from
[`model_backward`](https://nuwani-palihawadana.github.io/smimodel/reference/model_backward.md).

## Usage

``` r
eliminate(
  ind,
  train,
  val,
  yvar,
  family = gaussian(),
  s.vars = NULL,
  s.basedim = NULL,
  linear.vars = NULL,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colRange = NULL
)
```

## Arguments

- ind:

  An `integer` corresponding to the position of the predictor variable
  to be eliminated when fitting the model. (i.e. the function will
  combine `s.vars` and `linear.vars` in a single vector and eliminate
  the element corresponding to `ind`.)

- train:

  The data set on which the model(s) will be trained. Must be a data set
  of class `tsibble`.

- val:

  Validation data set. (The data set on which the model selection will
  be performed.) Must be a data set of class `tsibble`.

- yvar:

  Name of the response variable as a character string.

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

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string. (Since the nonlinear
  functions are estimated using splines, extrapolation is not desirable.
  Hence, if any predictor variable in `val` that is treated non-linearly
  in the estimated model, will be truncated to be in the in-sample range
  before obtaining predictions. If any variables are listed here will be
  excluded from such truncation.)

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
