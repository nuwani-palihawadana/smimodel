#' Sparse Multiple Index (SMI) model estimation
#'
#' Fits a single nonparametric multiple index model to the data. This is a
#' helper function designed to be called from `smimodel()`.
#'
#' @param data Training data set on which models will be trained. Should be a
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
#'   models (number of indices = `num_ind`; `num_models` random instances of the
#'   model (i.e. the predictor assignment to indices and initial index
#'   coefficients are generated randomly) are considered), and the final optimal
#'   model with the lowest loss is returned, and "userInput" - user specifies
#'   the initial model structure (i.e. the number of indices and the placement
#'   of index variables among indices) and the initial index coefficients
#'   through `index.ind` and `index.coefs` arguments respectively.
#' @param num_ind If `initialise = "ppr" or "multiple"`: an integer that
#'   specifies the number of indices to be used in the model(s). The default is
#'   5.
#' @param num_models If `initialise = "multiple"`: an integer that specifies the
#'   number of starting models to be checked. The default is 5.
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
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
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

smimodel.fit <- function(data, yvar, neighbour = 0, 
                         family = gaussian(), index.vars, 
                         initialise = c("ppr", "additive", "linear", 
                                        "multiple", "userInput"),
                         num_ind = 5, num_models = 5, seed = 123, index.ind = NULL, 
                         index.coefs = NULL, s.vars = NULL, linear.vars = NULL, 
                         lambda0 = 1, lambda2 = 1, 
                         M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                         TimeLimit = Inf, MIPGap = 1e-4, 
                         NonConvex = -1, verbose = FALSE){
  stopifnot(tsibble::is_tsibble(data))
  data_index <- index(data)
  data_key <- key(data)[[1]]
  data1 <- data |>
    tibble::as_tibble() |>
    dplyr::arrange({{data_index}})
  initialise <- match.arg(initialise)
  if(initialise == "ppr"){
    scaled <- scaling(data = data1, index.vars = index.vars)
    scaledInfo <- scaled$scaled_info
    data1 <- scaled$scaled_data
    pre.formula <- lapply(index.vars, function(var) paste0(var)) |> 
      paste(collapse = "+") 
    pre.formula <- paste(yvar, "~", pre.formula)
    # Fitting PPR
    ppr_fit <- stats::ppr(as.formula(pre.formula), data = data1, 
                          nterms = num_ind)
    # Sparsifying indices
    ppr_coefs <- ppr_fit$alpha
    threshold <- max(abs(ppr_fit$alpha))*0.1
    zero_ind <- which(abs(ppr_fit$alpha) < threshold)
    ppr_coefs[zero_ind] <- 0 
    for(i in 1:NROW(ppr_coefs)){
      maxCoef <- which.max(abs(ppr_coefs[i, ]))
      ppr_coefs[i, ][-maxCoef] <- 0 
    }
    index.ind <- vector(mode = "list", length = NCOL(ppr_coefs))
    index.coefs <- vector(mode = "list", length = NCOL(ppr_coefs))
    for(a in 1:NCOL(ppr_coefs)){
      index.ind[[a]] <- rep(a, NROW(ppr_coefs))
      index.coefs[[a]] <- ppr_coefs[ , a]
    }
    index.ind <- unlist(index.ind)
    index.coefs <- unlist(index.coefs)
    names(index.coefs) <- NULL
    data1 <- data1 |>
      tsibble::as_tsibble(index = {{data_index}}, key = {{data_key}})
    # Constructing the initial `smimodel`
    init_smimodel <- new_smimodelFit(data = data1, yvar = yvar, 
                                     neighbour = neighbour,
                                     family = family,
                                     index.vars = index.vars, 
                                     initialise = "userInput",  
                                     index.ind = index.ind, 
                                     index.coefs = index.coefs, 
                                     s.vars = s.vars,
                                     linear.vars = linear.vars)
    # Optimising the initial `smimodel`
    opt_smimodel_temp <- update_smimodelFit(object = init_smimodel, data = data1, 
                                            lambda0 = lambda0, lambda2 = lambda2, 
                                            M = M, max.iter = max.iter, 
                                            tol = tol, tolCoefs = tolCoefs,
                                            TimeLimit = TimeLimit, 
                                            MIPGap = MIPGap, NonConvex = NonConvex,
                                            verbose = verbose)
    opt_smimodel <- unscaling(object = opt_smimodel_temp, scaledInfo = scaled)
  }else if(initialise == "multiple"){
    Y_data <- as.matrix(data[ , yvar])
    num_pred <- length(index.vars)
    smimodels_initial <- vector(mode = "list", length = num_models)
    smimodels_optimised <- vector(mode = "list", length = num_models)
    smimodels_loss <- vector(mode = "list", length = num_models)
    # Multiple starting points
    permutes <- gtools::permutations(num_ind, num_ind)
    set.seed(seed)
    for(j in 1:num_models){
      print(paste0("Multiple starting points: ", j))
      permute_ind <- sample(1:dim(permutes)[1], 1)
      rest <- sample(1:num_ind, (num_pred - num_ind), replace = TRUE)
      indexInd <- c(permutes[permute_ind, ], rest)
      indexCoefs <- runif(num_pred)
      smimodels_initial[[j]] <- new_smimodelFit(data = data, yvar = yvar, 
                                                neighbour = neighbour,
                                                family = family,
                                                index.vars = index.vars, 
                                                initialise = "userInput", 
                                                index.ind = indexInd, 
                                                index.coefs = indexCoefs, 
                                                s.vars = s.vars,
                                                linear.vars = linear.vars)
      smimodels_optimised[[j]] <- update_smimodelFit(object = smimodels_initial[[j]], 
                                                     data = data, 
                                                     lambda0 = lambda0, 
                                                     lambda2 = lambda2, 
                                                     M = M, max.iter = max.iter, 
                                                     tol = tol, tolCoefs = tolCoefs,
                                                     TimeLimit = TimeLimit, 
                                                     MIPGap = MIPGap,
                                                     NonConvex = NonConvex,
                                                     verbose = verbose)
      # Preparing alpha - index coefficients vector
      list_index <- smimodels_optimised[[j]]$alpha
      numInd <- NCOL(list_index)
      alpha <- vector(mode = "list", length = numInd)
      for(k in 1:numInd){
        alpha[[k]] <- list_index[ , k]
      }
      alpha <- unlist(alpha)
      names(alpha) <- NULL
      # Calculating loss
      smimodels_loss[[j]] <- loss(Y = Y_data, 
                                  Yhat = smimodels_optimised[[j]]$gam$fitted.values, 
                                  alpha = alpha, 
                                  lambda0 = lambda0, lambda2 = lambda2)
    }
    min_loss <- which(unlist(smimodels_loss) == min(unlist(smimodels_loss)))
    init_smimodel <- smimodels_initial[[min_loss]]
    opt_smimodel <- smimodels_optimised[[min_loss]]
  }else{
    # Constructing the initial `smimodel`
    init_smimodel <- new_smimodelFit(data = data, yvar = yvar, 
                                     neighbour = neighbour,
                                     family = family,
                                     index.vars = index.vars, 
                                     initialise = initialise, 
                                     index.ind = index.ind, 
                                     index.coefs = index.coefs, 
                                     s.vars = s.vars,
                                     linear.vars = linear.vars)
    # Optimising the initial `smimodel`
    opt_smimodel <- update_smimodelFit(object = init_smimodel, data = data,
                                       lambda0 = lambda0, lambda2 = lambda2, 
                                       M = M, max.iter = max.iter, 
                                       tol = tol, tolCoefs = tolCoefs,
                                       TimeLimit = TimeLimit, 
                                       MIPGap = MIPGap, NonConvex = NonConvex,
                                       verbose = verbose)
  }
  output <- list("initial" = init_smimodel, "best" = opt_smimodel)
  return(output)
}
utils::globalVariables(".")