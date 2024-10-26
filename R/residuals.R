#' Extract residuals from a fitted \code{smimodel}
#'
#' Generates residuals from a fitted \code{smimodel} object.
#'
#' @param object A \code{smimodel} object.
#' @param ... Other arguments not currently used.
#' @return A \code{numeric} vector of residuals.
#'
#' @method residuals smimodel
#'
#' @export
residuals.smimodel <- function(object, ...){
  resids <- augment(object)$.resid
  return(resids)
}