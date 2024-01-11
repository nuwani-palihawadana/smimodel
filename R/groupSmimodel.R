#' Sparse Multiple Index (SMI) models based on a grouping variable
#'
#' Fits nonparametric multiple index models, with simultaneous variable
#' selection for each group based on a grouping variable of interest.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class `tsibble`.(Make sure there are no additional
#'   date/time/date-time/yearmonth/POSIXct/POSIXt variables except for the
#'   `index` of the `tsibble`). If multiple models are fitted, the grouping
#'   variable should be the key of the `tsibble`.
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
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
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
#' @importFrom dplyr arrange filter mutate rename
#' @importFrom tsibble as_tsibble is_tsibble index key
#' @importFrom vctrs vec_as_names
#'
#' @export
groupSmimodel <- function(data, yvar, index.vars, 
                          initialise = c("additive", "linear", "userInput"),
                          index.ind = NULL, index.coefs = NULL, 
                          neighbour = 0, linear.vars = NULL, 
                          lambda0 = 1, lambda2 = 1, 
                          M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                          TimeLimit = Inf, verbose = FALSE){
  initialise <- match.arg(initialise)
  # Constructing the initial `groupSmimodel`
  init_groupSmimodel <- new_groupSmimodel(data = data, yvar = yvar, 
                                          index.vars = index.vars, 
                                          initialise = initialise, 
                                          index.ind = index.ind, 
                                          index.coefs = index.coefs,
                                          neighbour = 0,
                                          linear.vars = linear.vars)
  # Optimising the initial `groupSmimodel`
  opt_groupSmimodel <- update_groupSmimodel(object = init_groupSmimodel, 
                                            data = data, neighbour = neighbour, 
                                            lambda0 = lambda0, lambda2 = lambda2, 
                                            M = M, max.iter = max.iter, 
                                            tol = tol, tolCoefs = tolCoefs,
                                            TimeLimit = TimeLimit, verbose = verbose)
  output <- list("initial" = init_groupSmimodel, "best" = opt_groupSmimodel)
  return(output)
}