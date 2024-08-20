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
  # # NULL prediction interval bounds for now
  # nint <- length(level)
  # lower <- matrix(NA, ncol = nint, nrow = h)
  # upper <- lower
  # for (i in 1:nint) {
  #   lower[, i] <- NA
  #   upper[, i] <- NA
  # }
  # colnames(lower) <- colnames(upper) <- paste(level, "%", sep = "")
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        # lower = lower,
                        # upper = upper,
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}


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