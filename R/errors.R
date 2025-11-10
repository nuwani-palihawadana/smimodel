#' @rdname point_measures
#' @export
MAE <- function(residuals, na.rm = TRUE, ...){
  mean(abs(residuals), na.rm = na.rm)
}

#' @rdname point_measures
#' @export
MSE <- function(residuals, na.rm = TRUE, ...){
  mean(residuals^2, na.rm = na.rm)
}


#' Point estimate accuracy measures
#'
#' @param residuals A vector of residuals from either the validation or test
#'   data.
#' @param na.rm If \code{TRUE}, remove the missing values before calculating the
#'   accuracy measure.
#' @param ... Additional arguments for each measure.
#'
#' @return For the individual functions (`MAE`, `MSE`), returns a single numeric
#' scalar giving the requested accuracy measure.
#'
#' For the exported object `point_measures`, returns a named list of functions
#' that can be supplied to higher-level accuracy routines.
#' 
#' @examples
#' set.seed(123)
#' ytrain <- rnorm(100)
#' ytest  <- rnorm(30)
#' yhat   <- ytest + rnorm(30, sd = 0.3)
#' resid   <- ytest - yhat
#'
#' MAE(resid)
#' MSE(resid)
#'
#' @name point_measures
#' @export
point_measures <- list(MAE = MAE, MSE = MSE)
