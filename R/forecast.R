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
#' @export
forecast.gaimFit <- function(object, h = 1, level = c(80, 95), newdata, 
                             recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Groupwise Additive Index Model"
  pred <- predict(object = object, newdata = newdata, recursive = recursive,
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
#' @export
forecast.pprFit <- function(object, h = 1, level = c(80, 95), newdata, 
                             recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Projection Pursuit Regression Model"
  pred <- predict(object = object, newdata = newdata, recursive = recursive,
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
#' @export
forecast.backward <- function(object, h = 1, level = c(80, 95), newdata, 
                            recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Nonparametric Additive Model with Backward Elimination"
  pred <- predict(object = object, newdata = newdata, recursive = recursive,
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
#' @export
forecast.gamFit <- function(object, h = 1, level = c(80, 95), newdata, 
                            recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Generalised Additive Model"
  pred <- predict(object = object, newdata = newdata, recursive = recursive,
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