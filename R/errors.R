#' @rdname point_measures
#' @export
MSE <- function(residuals, na.rm = TRUE, ...){
  mean(residuals^2, na.rm = na.rm)
}


#' @rdname point_measures
#' @export
MAE <- function(residuals, na.rm = TRUE, ...){
  mean(abs(residuals), na.rm = na.rm)
}


#' Point estimate accuracy measures
#'
#' @param residuals A vector of residuals from either the validation or test
#'   data.
#' @param na.rm If \code{TRUE}, remove the missing values before calculating the
#'   accuracy measure.
#' @param ... Additional arguments for each measure.
#'
#' @export
point_measures <- list(MAE = MAE, MSE = MSE)
