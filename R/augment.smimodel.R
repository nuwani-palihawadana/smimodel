#' Augment function for class `smimodel`
#'
#' Generates residuals and fitted values of a fitted `smimodel` object.
#'
#' @param x A `smimodel` object.
#' @param ... Other arguments not currently used.
#'
#' @importFrom generics augment
#'
#' @method augment smimodel
#'
#' @export
augment.smimodel <- function(x, ...) {
  smimodel.resid <- x$gam$residuals
  smimodel.fitted <- x$gam$fitted.values
  df <- tibble::tibble(
    .resid = smimodel.resid,
    .fitted = smimodel.fitted
  )
  return(df)
}
#' @export
generics::augment