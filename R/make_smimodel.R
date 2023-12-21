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
  # Derivatives of the fitted smooths
  dgz <- vector(length = length(index.names), mode = "list")
  dgz_names <- character(length = length(index.names))
  for (i in seq_along(index.names)) {
    temp <- gratia::derivatives(x, type = "central", data = index.data, 
                        term = paste0("s(", paste0(index.names[i]), ")"))
    dgz[[i]] <- temp$derivative
    dgz_names[i] <- paste0("d", i)
  }
  names(dgz) <- dgz_names
  # Constructing a new index coefficient vector to have all predictors in each index
  ind_pos <- split(1:length(index.ind), index.ind)
  new_index_info <- allpred_index(num_pred = length(index.vars), 
                                  num_ind = length(ind_pos), 
                                  ind_pos = ind_pos, 
                                  alpha = alpha)
  new_alpha <- split(new_index_info$alpha_init_new, new_index_info$index)
  # Constructing the class `smimodel`
  smimodel <- vector(mode = "list", length = (length(new_alpha)+3))
  for(i in 1:length(new_alpha)){
    smimodel[[i]] <- list("coefficients" = new_alpha[[i]],
                          "derivatives" = dgz[[i]])
  }
  smimodel[[(length(new_alpha)+1)]] <- yvar
  smimodel[[(length(new_alpha)+2)]] <- index.vars
  if(!is.null(linear.vars)){
    smimodel[[(length(new_alpha)+3)]] <- linear.vars
  }
  names(smimodel) <- c(index.names, "var_y","vars_index", "vars_linear")
  class(smimodel) <- c("smimodel", "list")
  return(smimodel)
}