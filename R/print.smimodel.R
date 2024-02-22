#' Printing a `smimodel` object
#'
#' The default print method for a `smimodel` object.
#'
#' @param x A model object of class `smimodel` as produced by
#'   `new_smimodel()` or `update_smimodel()` functions.
#' @param ... Other arguments not currently used.
#'
#' @method print smimodel
#'
#' @export
print.smimodel <- function(x, ...) {
  cat("Fitted smimodel:\n")
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
  cat("Fitted gam:")
  print(x$gam)
  cat("\n")
}