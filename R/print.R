#' Printing a `smimodel` object
#'
#' The default print method for a `smimodel` object.
#'
#' @param x A model object of class `smimodel` as produced by `model_smimodel()`
#'   function.
#' @param ... Other arguments not currently used.
#'
#' @method print smimodel
#'
#' @export
print.smimodel <- function(x, ...) {
  cat("Fitted SMI Model(s):\n")
  NextMethod()
}



#' Printing a `smimodelFit` object
#'
#' The default print method for a `smimodelFit` object.
#'
#' @param x A model object of class `smimodelFit` as produced by
#'   `new_smimodelFit()` or `update_smimodelFit()` functions.
#' @param ... Other arguments not currently used.
#'
#' @method print smimodelFit
#'
#' @export
print.smimodelFit <- function(x, ...) {
  cat("SMI Model Fit:\n")
  cat("Index coefficients:\n")
  print(x$alpha)
  cat("\n")
  cat("Response variable:\n")
  print(x$var_y)
  cat("\n")
  cat("Index variables:\n")
  print(x$vars_index)
  cat("\n")
  cat("Spline variables (non-index):\n")
  print(x$vars_s)
  cat("\n")
  cat("Linear variables:\n")
  print(x$vars_linear)
  cat("\n")
  cat("GAM Fit:")
  print(x$gam)
  cat("\n")
}



#' Printing a `backward` object
#'
#' The default print method for a `backward` object.
#'
#' @param x A model object of class `backward` as produced by `model_backward()`
#'   function.
#' @param ... Other arguments not currently used.
#'
#' @method print backward
#'
#' @export
print.backward <- function(x, ...) {
  cat("Fitted Model(s):\n")
  NextMethod()
}