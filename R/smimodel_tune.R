#' Sparse Multiple Index (SMI) model with a given penalty parameter combination
#'
#' Fits a nonparametric multiple index model to the data for a givcen
#' combination of the penalty parameters (lambda0, lambda2), and returns the
#' in-sample mean squared error (MSE). (Used within `greedy()`; users are not
#' expected to use this function directly.)
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param yvar Name of the response variable as a character string.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is "ppr", where the initial model is
#'   derived from projection pursuit regression. The other options are
#'   "additive" - nonparametric additive model, "linear" - linear regression
#'   model (i.e. a special case single-index model, where the initial values of
#'   the index coefficients are obtained through a linear regression),
#'   "multiple" - multiple models are fitted starting with different initial
#'   models (default `num_ind` (number of indices) = 5; five random instances of
#'   the model (i.e. the predictor assignment to indices and initial index
#'   coefficients are generated randomly) are considered), and the final optimal
#'   model with the lowest loss is returned (user can change the number of
#'   indices considered using `num_ind` argument), and "userInput" - user
#'   specifies the initial model structure (i.e. the number of indices and the
#'   placement of index variables among indices) and the initial index
#'   coefficients through `index.ind` and `index.coefs` arguments respectively.
#' @param num_ind If `initialise = "ppr" or "multiple"`: an integer that
#'   specifies the number of indices to be used in the model(s).
#' @param num_models If `initialise = "multiple"`: an integer that specifies the
#'   number of starting points to be checked.
#' @param seed If `initialise = "multiple"`: the seed to be set when generating
#'   random starting points.
#' @param index.ind If `initialise = "userInput"`: an integer vector that
#'   assigns group index for each predictor in `index.vars`.
#' @param index.coefs If `initialise = "userInput"`: a numeric vector of index
#'   coefficients.
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#' @param lambda.comb A numeric vector (of length two) indicating the values for
#'   the two penalty parameters lambda0 and lambda2.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param verbose The option to print detailed solver output.
#' 
#' @importFrom fabletools MSE
#'
#' @export
smimodel_tune <- function(data, yvar, index.vars, 
                          initialise = c("ppr", "additive", "linear", "multiple", "userInput"),
                          num_ind = 5, num_models = 5, seed = 123, index.ind = NULL, 
                          index.coefs = NULL, linear.vars = NULL, 
                          lambda.comb = c(1, 1), 
                          M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                          TimeLimit = Inf, MIPGap = 1e-4, verbose = FALSE){
  # Estimating smimodel
  smimodel <- smimodel(data = data, yvar = yvar, 
                       index.vars = index.vars, 
                       initialise = initialise, 
                       num_ind = num_ind, num_models = num_models, 
                       seed = seed,
                       index.ind = index.ind, 
                       index.coefs = index.coefs,
                       linear.vars = linear.vars,
                       lambda0 = lambda.comb[[1]], lambda2 = lambda.comb[[2]], 
                       M = M, max.iter = max.iter, 
                       tol = tol, tolCoefs = tolCoefs,
                       TimeLimit = TimeLimit, MIPGap = MIPGap,
                       verbose = verbose)
  # In-sample MSE
  smimodel_mse <- fabletools::MSE(.resid = smimodel$best$gam$residuals)
  return(smimodel_mse)
}