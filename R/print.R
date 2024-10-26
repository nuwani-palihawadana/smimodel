#' Printing a \code{smimodel} object
#'
#' The default print method for a \code{smimodel} object.
#'
#' @param x An object of class \code{smimodel}.
#' @param ... Other arguments not currently used.
#'
#' @method print smimodel
#'
#' @export
print.smimodel <- function(x, ...) {
  cat("Fitted SMI Model(s):\n")
  NextMethod()
}


#' Printing a \code{smimodelFit} object
#'
#' The default print method for a \code{smimodelFit} object.
#'
#' @param x An object of class \code{smimodelFit}.
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


#' Printing a \code{backward} object
#'
#' The default print method for a \code{backward} object.
#'
#' @param x A model object of class \code{backward} as produced by
#'   \code{model_backward()} function.
#' @param ... Other arguments not currently used.
#'
#' @method print backward
#'
#' @export
print.backward <- function(x, ...) {
  cat("Fitted Model(s):\n")
  NextMethod()
}


#' Printing a \code{pprFit} object
#'
#' The default print method for a \code{pprFit} object.
#'
#' @param x A model object of class \code{pprFit} as produced by
#'   \code{model_ppr()} function.
#' @param ... Other arguments not currently used.
#'
#' @method print pprFit
#'
#' @export
print.pprFit <- function(x, ...) {
  cat("Fitted Model(s):\n")
  NextMethod()
}


#' Printing a \code{gaimFit} object
#'
#' The default print method for a \code{gaimFit} object.
#'
#' @param x A model object of class \code{gaimFit} as produced by
#'   \code{model_gaim()} function.
#' @param ... Other arguments not currently used.
#'
#' @method print gaimFit
#'
#' @export
print.gaimFit <- function(x, ...) {
  cat("Fitted Model(s):\n")
  NextMethod()
}