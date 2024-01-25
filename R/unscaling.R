#' Unscale a fitted `smimodel`
#'
#' Transforms back the index coefficients to suit original-scale index variables
#' if the same were standardised when estimating the `smimodel` (happens in
#' `initialise = "ppr"` in `smimodel()`). Users are not expected to directly use
#' this function; usually called within `smimodel()`.
#'
#' @param object A `smimodel` object.
#' @param scaledInfo The list returned from a call of the function `scaling()`.
#'   (Relates to the argument `data` in the corresponding call of `smimodel()`.)
#'
#' @export
unscaling <- function(object, scaledInfo){
  scaledInfo <- scaledInfo$scaled_info
  list_index <- object$alpha[ , 2:NCOL(object$alpha)]
  for(b in 1:length(list_index)){
    temp_coef <- list_index[ , b]/scaledInfo
    object$alpha[ , (b+1)] <- temp_coef
    # temp_intercept <- ((list_index[[b]]$coefficients*scaledInfo$scaled_means*(-1))/scaledInfo$scaled_scales)
    # object[[b]]$intercept <- sum(temp_intercept)
  }
  return(object)
}