#' Sparse Multiple Index (SMI) models - Alternative function (trial for changes
#' in algorithm)
#'
#' Fits a nonparametric multiple index model to data with simultaneous variable
#' selection (hence "sparse").
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param yvar Name of the response variable as a character string.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is "additive", where the initial model
#'   will be a nonparametric additive model. The other options are "linear" -
#'   linear regression model (i.e. a special case single-index model, where the
#'   initial values of the index coefficients are obtained through a linear
#'   regression), and "userInput" - user specifies the initial model structure
#'   (i.e. the number of indices and the placement of index variables among
#'   indices) and the initial index coefficients through `index.ind` and
#'   `index.coefs` arguments respectively.
#' @param index.ind If `initialise = "userInput"`: an integer vector that
#'   assigns group index for each predictor in `index.vars`.
#' @param index.coefs If `initialise = "userInput"`: a numeric vector of index
#'   coefficients.
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param verbose The option to print detailed solver output.
#'
#' @export
alt_smimodel <- function(data, yvar, index.vars, 
                         initialise = c("additive", "linear", "userInput"),
                         index.ind = NULL, index.coefs = NULL, 
                         linear.vars = NULL, lambda0 = 1, lambda2 = 1, 
                         M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                         TimeLimit = Inf, verbose = FALSE){
  initialise <- match.arg(initialise)
  # Constructing the initial `smimodel`
  init_smimodel <- new_smimodel(data = data, yvar = yvar, 
                                index.vars = index.vars, 
                                initialise = initialise,
                                index.ind = index.ind, 
                                index.coefs = index.coefs, 
                                linear.vars = linear.vars)
  # Optimising the initial `smimodel`
  opt_smimodel <- alt_update_smimodel(object = init_smimodel, data = data, 
                                  lambda0 = lambda0, lambda2 = lambda2, 
                                  M = M, max.iter = max.iter, 
                                  tol = tol, tolCoefs = tolCoefs,
                                  TimeLimit = TimeLimit, verbose = verbose)
  output <- list("initial" = init_smimodel, "best" = opt_smimodel)
  return(output)
}