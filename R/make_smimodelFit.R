#' Converting a fitted \code{gam} object to a \code{smimodelFit} object
#'
#' Converts a given object of class \code{gam} to an object of class
#' \code{smimodelFit}.
#'
#' @param x A fitted \code{gam} object.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an \code{integer}. If \code{neighbour =
#'   x}, \code{x} number of keys before the key of interest and \code{x} number
#'   of keys after the key of interest are grouped together for model fitting.
#'   The default is \code{neighbour = 0} (i.e. no neighbours are considered for
#'   model fitting).
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices are estimated.
#' @param index.ind An \code{integer} vector that assigns group index for each
#'   predictor in \code{index.vars}.
#' @param index.data A \code{tibble} including columns for the constructed
#'   indices.
#' @param index.names A \code{character} vector of names of the constructed
#'   indices.
#' @param alpha A vector of index coefficients.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines are fitted individually (rather than considering as part
#'   of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that are included linearly in the model.
#' @return An object of class \code{smimodelFit}, which is a list that contains
#' following elements: \item{alpha}{A sparse matrix of index coefficients vectors.
#' Each column of the matrix corresponds to the index coefficient vector of each
#' index.} \item{derivatives}{A \code{tibble} of derivatives of the estimated
#' smooths.} \item{var_y}{Name of the response variable.}
#' \item{vars_index}{A \code{character} vector of names of the predictor
#'   variables for which indices are estimated.}
#'   \item{vars_s}{A \code{character} vector of names of the predictor variables
#'   for which splines are fitted individually.}
#'   \item{vars_linear}{A \code{character} vector of names of the predictor
#'   variables that are included linearly in the model.}
#'   \item{neighbour}{Number of neighbours of each key considered in model
#'   fitting.} \item{gam}{Fitted \code{gam}.}

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
                                    select = paste0("s(", index.names[i], ")"))
        dgz[[i]] <- temp$.derivative
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