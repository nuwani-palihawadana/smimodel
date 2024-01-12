#' Sparse Multiple Index (SMI) models
#'
#' Fits a nonparametric multiple index model to the data, with simultaneous
#' variable selection (hence "sparse").
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
#' @importFrom stats runif
#'
#' @export
smimodel_multistart <- function(data, yvar, index.vars, 
                                initialise = c("additive", "linear", "multiple", "userInput"),
                                index.ind = NULL, index.coefs = NULL, linear.vars = NULL, 
                                lambda0 = 1, lambda2 = 1, 
                                M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                                TimeLimit = Inf, verbose = FALSE){
  initialise <- match.arg(initialise)
  if(initialise == "multiple"){
    Y_data <- as.matrix(data[ , yvar])
    num_pred <- length(index.vars)
    # Different number of indices to be considered (The single-index model 
    # (i.e. numInds = 1) is initialised using the linear model.)
    numInds <- c(1, 2, 3, 5)
    smimodels_initial <- vector(mode = "list", length = 4)
    smimodels_optimised <- vector(mode = "list", length = 4)
    smimodels_loss <- vector(mode = "list", length = 4)
    # Starting point 1: Linear model (Single-index)
    print("Multiple starting points: 1")
    smimodels_initial[[1]] <- new_smimodel(data = data, yvar = yvar, 
                                           index.vars = index.vars, 
                                           initialise = "linear", 
                                           index.ind = index.ind, 
                                           index.coefs = index.coefs, 
                                           linear.vars = linear.vars)
    smimodels_optimised[[1]] <- update_smimodel(object = smimodels_initial[[1]], 
                                                data = data, 
                                                lambda0 = lambda0, 
                                                lambda2 = lambda2, 
                                                M = M, max.iter = max.iter, 
                                                tol = tol, tolCoefs = tolCoefs,
                                                TimeLimit = TimeLimit, 
                                                verbose = verbose)
    # Preparing alpha - index coefficients vector
    list_index <- smimodels_optimised[[1]][1:(length(smimodels_optimised[[1]])-4)]
    num_ind <- length(list_index)
    alpha <- vector(mode = "list", length = num_ind)
    for(i in 1:num_ind){
      alpha[[i]] <- list_index[[i]]$coefficients
    }
    alpha <- unlist(alpha)
    # Calculating loss
    smimodels_loss[[1]] <- LossFunction(Y = Y_data, 
                                        Yhat = smimodels_optimised[[1]]$gam$fitted.values, 
                                        alpha = alpha, 
                                        lambda0 = lambda0, lambda2 = lambda2)
    # Other starting points
    for(j in 2:4){
      print(paste0("Multiple starting points: ", j))
      set.seed(123)
      temp <- runif(numInds[j], min = 0.01, max = 0.99)
      prob <- temp/norm(matrix(temp, ncol = 1))
      indexInd <- sample(1:numInds[j], num_pred, replace = TRUE, prob = prob)
      indexCoefs <- runif(num_pred)
      smimodels_initial[[j]] <- new_smimodel(data = data, yvar = yvar, 
                                             index.vars = index.vars, 
                                             initialise = "userInput", 
                                             index.ind = indexInd, 
                                             index.coefs = indexCoefs, 
                                             linear.vars = linear.vars)
      smimodels_optimised[[j]] <- update_smimodel(object = smimodels_initial[[j]], 
                                                  data = data, 
                                                  lambda0 = lambda0, 
                                                  lambda2 = lambda2, 
                                                  M = M, max.iter = max.iter, 
                                                  tol = tol, tolCoefs = tolCoefs,
                                                  TimeLimit = TimeLimit, 
                                                  verbose = verbose)
      # Preparing alpha - index coefficients vector
      list_index <- smimodels_optimised[[j]][1:(length(smimodels_optimised[[j]])-4)]
      num_ind <- length(list_index)
      alpha <- vector(mode = "list", length = num_ind)
      for(k in 1:num_ind){
        alpha[[k]] <- list_index[[k]]$coefficients
      }
      alpha <- unlist(alpha)
      # Calculating loss
      smimodels_loss[[j]] <- LossFunction(Y = Y_data, 
                                          Yhat = smimodels_optimised[[j]]$gam$fitted.values, 
                                          alpha = alpha, 
                                          lambda0 = lambda0, lambda2 = lambda2)
    }
    min_loss <- which(unlist(smimodels_loss) == min(unlist(smimodels_loss)))
    init_smimodel <- smimodels_initial[[min_loss]]
    opt_smimodel <- smimodels_optimised[[min_loss]]
  }else{
    # Constructing the initial `smimodel`
    init_smimodel <- new_smimodel(data = data, yvar = yvar, 
                                  index.vars = index.vars, 
                                  initialise = initialise, 
                                  index.ind = index.ind, 
                                  index.coefs = index.coefs, 
                                  linear.vars = linear.vars)
    # Optimising the initial `smimodel`
    opt_smimodel <- update_smimodel(object = init_smimodel, data = data, 
                                    lambda0 = lambda0, lambda2 = lambda2, 
                                    M = M, max.iter = max.iter, 
                                    tol = tol, tolCoefs = tolCoefs,
                                    TimeLimit = TimeLimit, verbose = verbose)
  }
  output <- list("initial" = init_smimodel, "best" = opt_smimodel)
  return(output)
}