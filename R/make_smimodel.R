#' Converting a fitted `gam` object to a `smimodel` object
#'
#' Converts a given object of class `gam` to an object of class `smimodel`.
#'
#' @param x A fitted `gam` object.
#' @param yvar Name of the response variable as a character string.
#' @param index.vars A character vector of names of the predictor variables for
#'   which are estimated.
#' @param index.ind An integer vector that assigns group index for each
#'   predictor in `index.vars`.
#' @param index.data A `tibble` including columns for the constructed indices.
#' @param index.names A character vector of names of the constructed indices.
#' @param alpha A vector of index coefficients.
#' @param linear.vars A character vector of names of the predictor variables
#'   that are included linearly in the model.
#'
#' @importFrom gratia derivatives
#'
#' @export
make_smimodel <- function(x, yvar, index.vars, index.ind, index.data,
                          index.names, alpha, linear.vars = NULL){
  # Constructing a new index coefficient vector to have all predictors in each
  # index, and structuring output
  ind_pos <- split(1:length(index.ind), index.ind)
  names(alpha) <- NULL
  alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
  new_index_info <- allpred_index(num_pred = length(index.vars), 
                                  num_ind = length(ind_pos), 
                                  ind_pos = ind_pos, 
                                  alpha = alpha)
  new_alpha <- split(new_index_info$alpha_init_new, new_index_info$index)
  names(new_alpha) <- index.names
  var_names <- tibble::tibble(vars_index = index.vars)
  alpha <- dplyr::bind_cols(var_names, dplyr::bind_cols(new_alpha))
  # # Constructing the class `smimodel`
  # smimodel <- vector(mode = "list", length = (length(new_alpha)+4))
  # generating derivatives, and structuring the output
  if(!is.null(index.data)){
    if(length(x$smooth) == 0){
      # Derivatives of the fitted smooths
      dgz <- vector(length = length(index.names), mode = "list")
      for (i in seq_along(index.names)) {
        dgz[[i]] <- rep(1, NROW(index.data))
      }
    }else{
      # Derivatives of the fitted smooths
      dgz <- vector(length = length(index.names), mode = "list")
      for (i in seq_along(index.names)) {
        temp <- gratia::derivatives(x, type = "central",
                                    data = index.data,
                                    term = paste0("s(", index.names[i], ")"))
        dgz[[i]] <- temp$derivative
      }
    }
    names(dgz) <- paste0("d", seq_along(index.names))
    derivs <- dplyr::bind_cols(dgz)
    # for(i in 1:length(new_alpha)){
    #   smimodel[[i]] <- list("coefficients" = new_alpha[[i]],
    #                         "derivatives" = dgz[[i]])
    # }
  }else if(is.null(index.data)){
    derivs <- NULL
    # for(i in 1:length(new_alpha)){
    #   smimodel[[i]] <- list("coefficients" = new_alpha[[i]],
    #                         "derivatives" = NULL)
    # }
  }
  # smimodel[[(length(new_alpha)+1)]] <- yvar
  # smimodel[[(length(new_alpha)+2)]] <- index.vars
  # if(!is.null(linear.vars)){
  #   smimodel[[(length(new_alpha)+3)]] <- linear.vars
  # }
  # smimodel[[(length(new_alpha)+4)]] <- x
  smimodel <- list("alpha" = alpha, "derivatives" = derivs, "var_y" = yvar, 
                   "vars_index" = index.vars, "vars_linear" = linear.vars,
                   "gam" = x)
  # names(smimodel) <- c(index.names, "var_y", "vars_index", "vars_linear", "gam")
  class(smimodel) <- c("smimodel", "list")
  return(smimodel)
}