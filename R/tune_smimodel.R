#' Sparse Multiple Index (SMI) model with a given penalty parameter combination
#'
#' Fits a nonparametric multiple index model to the data for a given combination
#' of the penalty parameters (lambda0, lambda2), and returns the validation set
#' mean squared error (MSE). (Used within `greedy.fit()`; users are not expected
#' to use this function directly.)
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param val.data Validation data set. (The data set on which the penalty
#'   parameter selection will be performed.) Must be a data set of class
#'   `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
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
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted individually (rather than considering as a
#'   part of an index considered in `index.vars`).
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
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `val.data` to be filled with forecasts.
#'
#' @importFrom fabletools MSE

tune_smimodel <- function(data, val.data, yvar, neighbour = 0,
                          family = gaussian(), index.vars, 
                          initialise = c("ppr", "additive", "linear", 
                                         "multiple", "userInput"),
                          num_ind = 5, num_models = 5, seed = 123, 
                          index.ind = NULL, index.coefs = NULL, 
                          s.vars = NULL, linear.vars = NULL, lambda.comb = c(1, 1), 
                          M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                          TimeLimit = Inf, MIPGap = 1e-4, 
                          NonConvex = -1, verbose = FALSE, recursive = FALSE,
                          recursive_colRange = NULL){
  # Estimating smimodel
  smimodel <- smimodel.fit(data = data, yvar = yvar, 
                           neighbour = neighbour,
                           family = family,
                           index.vars = index.vars, 
                           initialise = initialise, 
                           num_ind = num_ind, num_models = num_models, 
                           seed = seed,
                           index.ind = index.ind, 
                           index.coefs = index.coefs,
                           s.vars = s.vars,
                           linear.vars = linear.vars,
                           lambda0 = lambda.comb[[1]], lambda2 = lambda.comb[[2]], 
                           M = M, max.iter = max.iter, 
                           tol = tol, tolCoefs = tolCoefs,
                           TimeLimit = TimeLimit, MIPGap = MIPGap,
                           NonConvex = NonConvex, verbose = verbose)
  # Validation set MSE
  # Predictions
  pred <- predict(object = smimodel, newdata = val.data, recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  # MSE
  smimodel_mse = MSE(.resid = (as.numeric(as.matrix(val.data[,{{yvar}}], ncol = 1)) - pred))
  # # In-sample MSE
  # smimodel_mse <- fabletools::MSE(.resid = smimodel$best$gam$residuals)
  return(smimodel_mse)
}