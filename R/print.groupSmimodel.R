#' Printing a `groupSmimodel` object
#'
#' The default print method for a `groupSmimodel` object.
#'
#' @param x A model object of class `groupSmimodel` as produced by
#'   `new_groupSmimodel()` or `update_groupSmimodel()` functions.
#' @param ... Other arguments not currently used.
#'
#' @method print groupSmimodel
#'
#' @export
print.groupSmimodel <- function(x, ...) {
  cat("Fitted groupSmimodel:\n")
  NextMethod()
}