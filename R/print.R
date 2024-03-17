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



#' Printing a `pprFit` object
#'
#' The default print method for a `pprFit` object.
#'
#' @param x A model object of class `pprFit` as produced by `model_ppr()`
#'   function.
#' @param ... Other arguments not currently used.
#'
#' @method print pprFit
#'
#' @export
print.pprFit <- function(x, ...) {
  cat("Fitted Model(s):\n")
  NextMethod()
}



#' Printing a `gaimFit` object
#'
#' The default print method for a `gaimFit` object.
#'
#' @param x A model object of class `gaimFit` as produced by `model_gaim()`
#'   function.
#' @param ... Other arguments not currently used.
#'
#' @method print gaimFit
#'
#' @export
print.gaimFit <- function(x, ...) {
  cat("Fitted Model(s):\n")
  NextMethod()
}