#' Plot estimated smooths from a fitted `smimodel`
#'
#' Plots the graphs of fitted spline(s) for a specified predictor variable, by
#' fixing all the other predictor values. If a set of multiple models are
#' fitted, plots individual graphs of fitted values of all the models in a
#' single plot except for tensor smooths fitted in `model_distlag()`. out of the
#' set of multiple models fitted)
#'
#' @param object A `smimodel` object returned from `model_smimodel()`.
#' @param model An `integer` to indicate the smooths of which model (out of the
#'   set of multiple models fitted) to be plotted.
#' @param ... Other arguments not currently used.
#' 
#' @importFrom ggplot2 autoplot
#' @importFrom gratia draw
#' 
#' @return Plot(s) of fitted spline(s).
#'
#' @method autoplot smimodel
#'
#' @export
autoplot.smimodel <- function(object, model = 1, ...){
  draw(object$fit[[model]]$best$gam)
}
#' @export
ggplot2::autoplot