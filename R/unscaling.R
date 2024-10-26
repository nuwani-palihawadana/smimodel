#' Unscale a fitted \code{smimodel}
#'
#' Transforms back the index coefficients to suit original-scale index variables
#' if the same were scaled when estimating the \code{smimodel} (happens in
#' \code{initialise = "ppr"} in \code{\link{model_smimodel}} or
#' \code{\link{greedy_smimodel}}). Users are not expected to directly use this
#' function; usually called within \code{\link{smimodel.fit}}.
#'
#' @param object A \code{smimodel} object.
#' @param scaledInfo The list returned from a call of the function
#'   \code{\link{scaling}}.

unscaling <- function(object, scaledInfo){
  scaledInfo <- scaledInfo$scaled_info
  list_index <- object$alpha
  for(b in 1:NCOL(list_index)){
    temp_coef <- list_index[ , b]/scaledInfo
    object$alpha[ , b] <- temp_coef
  }
  return(object)
}