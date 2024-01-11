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
  cat("Estimated index coefficients:\n")
  list_index <- x[1:(length(x)-4)]
  for(i in 1:length(list_index)){
    cat("index", i, "\n")
    print(x[[i]]$coefficients)
    cat("\n")
  }
  cat("Response variable:\n")
  print(x$var_y)
  cat("\n")
  cat("Index variables:\n")
  print(x$vars_index)
  cat("\n")
  cat("Linear variables:\n")
  print(x$vars_linear)
  cat("\n")
  cat("Fitted gam:")
  print(x$gam)
  cat("\n")
}