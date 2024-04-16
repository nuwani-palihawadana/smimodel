#' Extract residuals from a fitted `smimodel`
#'
#' Generates residuals from a fitted `smimodel` object.
#'
#' @param object A `smimodel` object.
#' @param ... Other arguments not currently used.
#'
#' @method residuals smimodel
#'
#' @export
residuals.smimodel <- function(object, ...){
  resids <- augment(object)$.resid
  return(resids)
}