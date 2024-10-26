#' Plot estimated smooths from a fitted \code{smimodel}
#'
#' Plots the graphs of fitted spline(s). If a set of multiple models are fitted,
#' plots graphs of fitted spline(s) of a specified model (in argument
#' \code{model}) out of the set of multiple models fitted.
#'
#' @param object A \code{smimodel} object.
#' @param model An \code{integer} to indicate the smooths of which model (out of
#'   the set of multiple models fitted) to be plotted.
#' @param ... Other arguments not currently used.
#'
#' @return Plot(s) of fitted spline(s).
#'
#' @method autoplot smimodel
#'
#' @export
autoplot.smimodel <- function(object, model = 1, ...){
  draw(object$fit[[model]]$best$gam)
}