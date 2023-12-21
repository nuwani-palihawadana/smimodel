#' Augment function for class `smimodel`
#'
#' Generates residuals and fitted values of a fitted `smimodel` object.
#'
#' @param x A `smimodel` object.
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param ... Other arguments not currently used.
#'
#' @importFrom generics augment
#'
#' @method augment smimodel
#'
#' @export
augment.smimodel <- function(x, data, ...) {
  fitted_gam <- make_gam(x = x, data = data)
  smimodel.resid <- fitted_gam$residuals
  smimodel.fitted <- fitted_gam$fitted.values
  df <- tibble::tibble(
    .resid = smimodel.resid,
    .fitted = smimodel.fitted
  )
  return(df)
}
#' @export
generics::augment