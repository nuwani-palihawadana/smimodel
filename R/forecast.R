#' Forecasting using SMI models
#'
#' Returns forecasts and other information for SMI models.
#'
#'
#' @param object An object of class \code{smimodel}. Usually the result of a call
#'   to \code{\link{model_smimodel}}.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{forecast}. Here, it is a list containing the
#' following elements: \item{method}{The name of the forecasting method as a
#' character string.} \item{model}{The fitted model.} \item{mean}{Point
#' forecasts as a time series.} \item{residuals}{Residuals from the fitted
#' model.} \item{fitted}{Fitted values (one-step forecasts).}
#'
#' @method forecast smimodel
#' 
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1015
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'   
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1010, ]
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' smimodel_ppr <- model_smimodel(data = sim_train,
#'                                yvar = "y",
#'                                index.vars = index.vars,
#'                                initialise = "ppr")
#'
#' forecast(smimodel_ppr, newdata = sim_test)
#' }
#'
#' @export
forecast.smimodel <- function(object, h = 1, level = c(80, 95), newdata, 
                             exclude.trunc = NULL,
                             recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Sparse Multiple Index Model"
  pred <- predict(object = object, newdata = newdata, 
                  exclude.trunc = exclude.trunc, recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  fitted <- augment(object)$.fitted
  residuals <- augment(object)$.resid
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}


#' Forecasting using GAIMs
#'
#' Returns forecasts and other information for GAIMs.
#'
#'
#' @param object An object of class \code{gaimFit}. Usually the result of a call
#'   to \code{\link{model_gaim}}.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{forecast}. Here, it is a list containing the
#' following elements: \item{method}{The name of the forecasting method as a
#' character string.} \item{model}{The fitted model.} \item{mean}{Point
#' forecasts as a time series.} \item{residuals}{Residuals from the fitted
#' model.} \item{fitted}{Fitted values (one-step forecasts).}
#'
#' @method forecast gaimFit
#' 
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1015
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'   
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1010, ]
#'
#' # Predictors taken as index variables
#' index.vars <- colnames(sim_data)[3:7]
#'
#' # Assign group indices for each predictor
#' index.ind = c(rep(1, 3), rep(2, 2))
#'
#' # Predictors taken as non-linear variables not entering indices
#' s.vars = "x_lag_005"
#'
#' # Model fitting
#' gaimModel <- model_gaim(data = sim_train,
#'                         yvar = "y",
#'                         index.vars = index.vars,
#'                         index.ind = index.ind,
#'                         s.vars = s.vars)
#'                         
#' forecast(gaimModel, newdata = sim_test)
#'
#' @export
forecast.gaimFit <- function(object, h = 1, level = c(80, 95), newdata, 
                             exclude.trunc = NULL,
                             recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Groupwise Additive Index Model"
  pred <- predict(object = object, newdata = newdata, 
                  exclude.trunc = exclude.trunc, recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  fitted <- augment(object)$.fitted
  residuals <- augment(object)$.resid
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}


#' Forecasting using PPR models
#'
#' Returns forecasts and other information for PPR models.
#'
#'
#' @param object An object of class \code{pprFit}. Usually the result of a call
#'   to \code{\link{model_ppr}}.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{forecast}. Here, it is a list containing the
#' following elements: \item{method}{The name of the forecasting method as a
#' character string.} \item{model}{The fitted model.} \item{mean}{Point
#' forecasts as a time series.} \item{residuals}{Residuals from the fitted
#' model.} \item{fitted}{Fitted values (one-step forecasts).}
#'
#' @method forecast pprFit
#' 
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1015
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'   
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1010, ]
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' pprModel <- model_ppr(data = sim_train,
#'                       yvar = "y",
#'                       index.vars = index.vars)
#'
#' forecast(pprModel, newdata = sim_test)
#'
#' @export
forecast.pprFit <- function(object, h = 1, level = c(80, 95), newdata, 
                            exclude.trunc = NULL,
                            recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Projection Pursuit Regression Model"
  pred <- predict(object = object, newdata = newdata, 
                  exclude.trunc = exclude.trunc,
                  recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  fitted <- augment(object)$.fitted
  residuals <- augment(object)$.resid
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}


#' Forecasting using nonparametric additive models with backward elimination
#'
#' Returns forecasts and other information for nonparametric additive models
#' with backward elimination.
#'
#'
#' @param object An object of class \code{backward}. Usually the result of a
#'   call to \code{\link{model_backward}}.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{forecast}. Here, it is a list containing the
#' following elements: \item{method}{The name of the forecasting method as a
#' character string.} \item{model}{The fitted model.} \item{mean}{Point
#' forecasts as a time series.} \item{residuals}{Residuals from the fitted
#' model.} \item{fitted}{Fitted values (one-step forecasts).}
#'
#' @method forecast backward
#' 
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1215
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Validation set
#' sim_val <- sim_data[1001:1200, ]
#' # Test set
#' sim_test <- sim_data[1201:1210, ]
#'
#' # Predictors taken as non-linear variables
#' s.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' backwardModel <- model_backward(data = sim_train,
#'                                 val.data = sim_val,
#'                                 yvar = "y",
#'                                 s.vars = s.vars)
#' forecast(backwardModel, newdata = sim_test)
#'
#' @export
forecast.backward <- function(object, h = 1, level = c(80, 95), newdata, 
                              exclude.trunc = NULL,
                              recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Nonparametric Additive Model with Backward Elimination"
  pred <- predict(object = object, newdata = newdata, 
                  exclude.trunc = exclude.trunc,
                  recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  fitted <- augment(object)$.fitted
  residuals <- augment(object)$.resid
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}


#' Forecasting using GAMs
#'
#' Returns forecasts and other information for GAMs.
#'
#' @param object An object of class \code{gamFit}. Usually the result of a call
#'   to \code{\link{model_gam}}.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{forecast}. Here, it is a list containing the
#' following elements: \item{method}{The name of the forecasting method as a
#' character string.} \item{model}{The fitted model.} \item{mean}{Point
#' forecasts as a time series.} \item{residuals}{Residuals from the fitted
#' model.} \item{fitted}{Fitted values (one-step forecasts).}
#'
#' @method forecast gamFit
#' 
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1015
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'   
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1010, ]
#'
#' # Predictors taken as non-linear variables
#' s.vars <- colnames(sim_data)[3:6]
#'
#' # Predictors taken as linear variables
#' linear.vars <- colnames(sim_data)[7:8]
#'
#' # Model fitting
#' gamModel <- model_gam(data = sim_train,
#'                       yvar = "y",
#'                       s.vars = s.vars,
#'                       linear.vars = linear.vars)
#'
#' forecast(gamModel, newdata = sim_test)
#'
#' @export
forecast.gamFit <- function(object, h = 1, level = c(80, 95), newdata, 
                            exclude.trunc = NULL,
                            recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Generalised Additive Model"
  pred <- predict(object = object, newdata = newdata, 
                  exclude.trunc = exclude.trunc,
                  recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  fitted <- augment(object)$.fitted
  residuals <- augment(object)$.resid
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}