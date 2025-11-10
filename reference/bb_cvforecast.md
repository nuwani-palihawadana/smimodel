# Single season block bootstrap prediction intervals through time series cross-validation forecasting

Compute prediction intervals by applying the single season block
bootstrap method to subsets of time series data using a rolling forecast
origin.

## Usage

``` r
bb_cvforecast(
  object,
  data,
  yvar,
  neighbour = 0,
  predictor.vars,
  h = 1,
  season.period = 1,
  m = 1,
  num.futures = 1000,
  level = c(80, 95),
  forward = TRUE,
  initial = 1,
  window = NULL,
  roll.length = 1,
  exclude.trunc = NULL,
  recursive = FALSE,
  recursive_colNames = NULL,
  na.rm = TRUE,
  ...
)
```

## Arguments

- object:

  Fitted model object of class `smimodel`, `backward`, `gaimFit` or
  `pprFit`.

- data:

  Data set. Must be a data set of class `tsibble`.(Make sure there are
  no additional date or time related variables except for the `index` of
  the `tsibble`). If multiple models are fitted, the grouping variable
  should be the `key` of the `tsibble`. If a `key` is not specified, a
  dummy key with only one level will be created.

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

- predictor.vars:

  A character vector of names of the predictor variables.

- h:

  Forecast horizon.

- season.period:

  Length of the seasonal period.

- m:

  Multiplier. (Block size = `NULLseason.period * m`)

- num.futures:

  Number of possible future sample paths to be generated.

- level:

  Confidence level for prediction intervals.

- forward:

  If `TRUE`, the final forecast origin for forecasting is \\y_T\\.
  Otherwise, the final forecast origin is \\y\_{T-1}\\.

- initial:

  Initial period of the time series where no cross-validation
  forecasting is performed.

- window:

  Length of the rolling window. If `NULL`, a rolling window will not be
  used.

- roll.length:

  Number of observations by which each rolling/expanding window should
  be rolled forward.

- exclude.trunc:

  The names of the predictor variables that should not be truncated for
  stable predictions as a character string. (Since the nonlinear
  functions are estimated using splines, extrapolation is not desirable.
  Hence, if any predictor variable is treated non-linearly in the
  estimated model, will be truncated to be in the in-sample range before
  obtaining predictions. If any variables are listed here will be
  excluded from such truncation.)

- recursive:

  Whether to obtain recursive forecasts or not (default - `FALSE`).

- recursive_colNames:

  If `recursive = TRUE`, a character vector giving the names of the
  columns in test data to be filled with forecasts.
  Recursive/autoregressive forecasting is required when the lags of the
  response variable itself are used as predictor variables into the
  model. Make sure such lagged variables are positioned together in
  increasing lag order (i.e. `lag_1, lag_2, ..., lag_m`, `lag_m =`
  maximum lag used) in `data`, with no break in the lagged variable
  sequence even if some of the intermediate lags are not used as
  predictors.

- na.rm:

  logical; if `TRUE` (default), any `NA` and `NaN`'s are removed from
  the sample before the quantiles are computed.

- ...:

  Other arguments not currently used.

## Value

An object of class `bb_cvforecast`, which is a list that contains
following elements:

- x:

  The original time series.

- method:

  A character string "bb_cvforecast".

- fit_times:

  The number of times the model is fitted in cross-validation.

- mean:

  Point forecasts as a multivariate time series, where the \\h^{th}\\
  column holds the point forecasts for forecast horizon \\h\\. The time
  index corresponds to the period for which the forecast is produced.

- res:

  The matrix of in-sample residuals produced in cross-validation. The
  number of rows corresponds to `fit_times`, and the row names
  corresponds the time index of the forecast origin of the corresponding
  cross-validation iteration.

- model_fit:

  Models fitted in cross-validation.

- level:

  The confidence values associated with the prediction intervals.

- lower:

  A list containing lower bounds for prediction intervals for each
  level. Each element within the list will be a multivariate time series
  with the same dimensional characteristics as `mean`.

- upper:

  A list containing upper bounds for prediction intervals for each
  level. Each element within the list will be a multivariate time series
  with the same dimensional characteristics as `mean`.

- possible_futures:

  A list of matrices containing future sample paths generated at each
  cross-validation step.

## See also

[`cb_cvforecast`](https://nuwani-palihawadana.github.io/smimodel/reference/cb_cvforecast.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
library(ROI)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1105
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(n)) |>
  mutate(
    # Add x_lags
    x_lag = lag_matrix(x_lag_000, 5)) |>
  unpack(x_lag, names_sep = "_") |>
  mutate(
    # Response variable
    y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 +
    (0.35*x_lag_002 + 0.7*x_lag_005)^2 + rnorm(n, sd = 0.1),
    # Add an index to the data set
    inddd = seq(1, n)) |>
  drop_na() |>
  select(inddd, y, starts_with("x_lag")) |>
  # Make the data set a `tsibble`
  as_tsibble(index = inddd)

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Training set
sim_train <- sim_data[1:1000, ]
# Test set
sim_test <- sim_data[1001:1100, ]

# Model fitting
smimodel_ppr <- model_smimodel(data = sim_train,
                               yvar = "y",
                               index.vars = index.vars,
                               initialise = "ppr")

# Block bootstrap prediction intervals (3-steps-ahead interval forecasts)
set.seed(12345)
smimodel_ppr_bb <- bb_cvforecast(object = smimodel_ppr,
                                 data = sim_data,
                                 yvar = "y",
                                 predictor.vars = index.vars,
                                 h = 3,
                                 num.futures = 50,
                                 window = 1000)
} # }
```
