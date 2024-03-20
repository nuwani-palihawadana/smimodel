#' Converting a fitted `gam` object to a `smimodelFit` object
#'
#' Converts a given object of class `gam` to an object of class `smimodelFit`.
#'
#' @param x A fitted `gam` object.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting.
#' @param index.vars A character vector of names of the predictor variables for
#'   which are estimated.
#' @param index.ind An integer vector that assigns group index for each
#'   predictor in `index.vars`.
#' @param index.data A `tibble` including columns for the constructed indices.
#' @param index.names A character vector of names of the constructed indices.
#' @param alpha A vector of index coefficients.
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted individually (rather than considering as a
#'   part of an index considered in `index.vars`).
#' @param linear.vars A character vector of names of the predictor variables
#'   that are included linearly in the model.
#'
#' @importFrom gratia derivatives
#' @importFrom methods as

make_smimodelFit <- function(x, yvar, neighbour, index.vars, index.ind, index.data,
                          index.names, alpha, s.vars = NULL, linear.vars = NULL){
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
  alpha <- as(as.matrix(dplyr::bind_cols(new_alpha)), "sparseMatrix")
  rownames(alpha) <- index.vars
  # Constructing the class `smimodelFit`
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
  }else if(is.null(index.data)){
    derivs <- NULL
  }
  smimodel <- list("alpha" = alpha, "derivatives" = derivs, "var_y" = yvar, 
                   "vars_index" = index.vars, "vars_s" = s.vars,
                   "vars_linear" = linear.vars, 
                   "neighbour" = neighbour, "gam" = x)
  class(smimodel) <- c("smimodelFit", "list")
  return(smimodel)
}