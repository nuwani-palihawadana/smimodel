#' Sparse Multiple Index (SMI) Models
#'
#' Fits nonparametric multiple index model(s), with simultaneous predictor
#' selection (hence "sparse") and predictor grouping. Possible to fit multiple
#' SMI models based on a grouping variable.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}). If
#'   multiple models are fitted, the grouping variable should be the \code{key}
#'   of the \code{tsibble}. If a \code{key} is not specified, a dummy key with
#'   only one level will be created.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an \code{integer}. If \code{neighbour =
#'   x}, \code{x} number of keys before the key of interest and \code{x} number
#'   of keys after the key of interest are grouped together for model fitting.
#'   The default is \code{neighbour = 0} (i.e. no neighbours are considered for
#'   model fitting).
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is \code{"ppr"}, where the initial model
#'   is derived from projection pursuit regression. The other options are
#'   \code{"additive"} - nonparametric additive model, \code{"linear"} - linear
#'   regression model (i.e. a special case single-index model, where the initial
#'   values of the index coefficients are obtained through a linear regression),
#'   \code{"multiple"} - multiple models are fitted starting with different
#'   initial models (number of indices = \code{num_ind}; \code{num_models}
#'   random instances of the model (i.e. the predictor assignment to indices and
#'   initial index coefficients are generated randomly) are considered), and the
#'   final optimal model with the lowest loss is returned, and
#'   \code{"userInput"} - user specifies the initial model structure (i.e. the
#'   number of indices and the placement of index variables among indices) and
#'   the initial index coefficients through \code{index.ind} and
#'   \code{index.coefs} arguments respectively.
#' @param num_ind If \code{initialise = "ppr"} or \code{"multiple"}: an
#'   \code{integer} that specifies the number of indices to be used in the
#'   model(s). The default is \code{num_ind = 5}.
#' @param num_models If \code{initialise = "multiple"}: an \code{integer} that
#'   specifies the number of starting models to be checked. The default is
#'   \code{num_models = 5}.
#' @param seed If \code{initialise = "multiple"}: the seed to be set when
#'   generating random starting points.
#' @param index.ind If \code{initialise = "userInput"}: an \code{integer} vector
#'   that assigns group index for each predictor in \code{index.vars}.
#' @param index.coefs If \code{initialise = "userInput"}: a \code{numeric}
#'   vector of index coefficients.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted individually (rather than considering as
#'   part of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value to be used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for the objective function value (loss) of MIP.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @return  An object of class \code{smimodel}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.}
#'   Each row of the column \code{fit} contains a list with two elements:
#'   \item{initial}{A list of information of the model initialisation. (For
#'   descriptions of the list elements see \code{\link{make_smimodelFit}}).}
#'   \item{best}{A list of information of the final optimised model. (For
#'   descriptions of the list elements see \code{\link{make_smimodelFit}}).}
#'
#' @details Sparse Multiple Index (SMI) model is a semi-parametric model that
#'   can be written as \deqn{y_{i} = \beta_{0} +
#' \sum_{j = 1}^{p}g_{j}(\boldsymbol{\alpha}_{j}^{T}\boldsymbol{x}_{ij}) +
#' \sum_{k = 1}^{d}f_{k}(w_{ik}) + \boldsymbol{\theta}^{T}\boldsymbol{u}_{i} +
#' \varepsilon_{i}, \quad i = 1, \dots, n,} where \eqn{y_{i}} is the univariate
#' response, \eqn{\beta_{0}} is the model intercept, \eqn{\boldsymbol{x}_{ij} \in
#' \mathbb{R}^{l_{j}}}, \eqn{j = 1, \dots, p} are \eqn{p} subsets of predictors
#'   entering indices, \eqn{\boldsymbol{\alpha}_{j}} is a vector of index
#'   coefficients corresponding to the index \eqn{h_{ij} =
#'   \boldsymbol{\alpha}_{j}^{T}\boldsymbol{x}_{ij}}, and \eqn{g_{j}} is a
#'   smooth nonlinear function (estimated by a penalised cubic regression
#'   spline). The model also allows for predictors that do not enter any
#'   indices, including covariates \eqn{w_{ik}} that relate to the response
#'   through nonlinear functions \eqn{f_{k}}, \eqn{k = 1, \dots, d}, and linear
#'   covariates \eqn{\boldsymbol{u}_{i}}.
#'
#'   In the model formulation related to this implementation, both the number of
#'   indices \eqn{p} and the predictor grouping among indices are assumed to be
#'   unknown prior to model estimation. Suppose we observe \eqn{y_1,\dots,y_n},
#'   along with a set of potential predictors,
#'   \eqn{\boldsymbol{x}_1,\dots,\boldsymbol{x}_n}, with each vector
#'   \eqn{\boldsymbol{x}_i} containing \eqn{q} predictors. This function
#'   implements algorithmic variable selection for index variables (i.e.
#'   predictors entering indices) of the SMI model by allowing for zero index
#'   coefficients for predictors. Non-overlapping predictors among indices are
#'   assumed (i.e. no predictor enters more than one index). For algorithmic
#'   details see reference.
#'
#' @references Palihawadana, N.K., Hyndman, R.J. & Wang, X. (2024). Sparse
#'   Multiple Index Models for High-Dimensional Nonparametric Forecasting.
#'   \url{https://www.monash.edu/business/ebs/research/publications/ebs/2024/wp16-2024.pdf}.
#'
#' @seealso \code{\link{greedy_smimodel}}
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y1, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' smimodel_ppr <- model_smimodel(data = sim_data,
#'                                yvar = "y1",
#'                                index.vars = index.vars,
#'                                initialise = "ppr")
#'
#' # Best (optimised) fitted model
#' smimodel_ppr$fit[[1]]$best
#' }
#'
#' @export
model_smimodel <- function(data, yvar, neighbour = 0, family = gaussian(), 
                           index.vars, 
                           initialise = c("ppr", "additive", "linear", 
                                          "multiple", "userInput"),
                           num_ind = 5, num_models = 5, seed = 123, 
                           index.ind = NULL, index.coefs = NULL, 
                           s.vars = NULL, linear.vars = NULL, 
                           lambda0 = 1, lambda2 = 1, 
                           M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                           TimeLimit = Inf, MIPGap = 1e-4, 
                           NonConvex = -1, verbose = FALSE){
  
  # Message for gurobi installation
  message("Do you have Gurobi solver installed? 
  Make sure you have an active installation of Gurobi solver (https://www.gurobi.com/) 
  in your local machine before using this function. 
  Refer the section 'Other Required Software' in the README for installation help.")
  
  # Check for `tsibble`
  stopifnot(tsibble::is_tsibble(data))
  
  initialise <- match.arg(initialise)
  data1 <- data
  data_index <- index(data1)
  data_key <- key(data1)
  if (length(key(data1)) == 0) {
    data1 <- data1 |>
      dplyr::mutate(dummy_key = rep(1, NROW(data1))) |>
      tsibble::as_tsibble(index = data_index, key = dummy_key)
    data_key <- key(data1)
  }
  key11 <- key(data1)[[1]]
  key_unique <- unique(as.character(sort(dplyr::pull((data1[, {{ key11 }}])[, 1]))))
  key_num <- seq_along(key_unique)
  ref <- data.frame(key_unique, key_num)
  data1 <- data1 |>
    dplyr::mutate(
      num_key = as.numeric(factor(as.character({{ key11 }}), levels = key_unique))
    )
  smimodels_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    df_cat <- df_cat |>
      drop_na()
    smimodels_list[[i]] <- smimodel.fit(data = df_cat, yvar = yvar, 
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
                                        lambda0 = lambda0, lambda2 = lambda2, 
                                        M = M, max.iter = max.iter, 
                                        tol = tol, tolCoefs = tolCoefs,
                                        TimeLimit = TimeLimit, MIPGap = MIPGap,
                                        NonConvex = NonConvex, verbose = verbose)
  }
  data_list <- list(key_unique, smimodels_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("smimodel", "tbl_df", "tbl", "data.frame")
  return(models)
}


#' SMI model estimation
#'
#' Fits a single nonparametric multiple index model to the data. This is a
#' helper function designed to be called from user-facing wrapper functions,
#' \code{\link{model_smimodel}} and \code{\link{greedy_smimodel}}.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}).
#' @param yvar Name of the response variable as a character string.
#' @param neighbour `neighbour` argument passed from the outer function.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is \code{"ppr"}, where the initial model
#'   is derived from projection pursuit regression. The other options are
#'   \code{"additive"} - nonparametric additive model, \code{"linear"} - linear
#'   regression model (i.e. a special case single-index model, where the initial
#'   values of the index coefficients are obtained through a linear regression),
#'   \code{"multiple"} - multiple models are fitted starting with different
#'   initial models (number of indices = \code{num_ind}; \code{num_models}
#'   random instances of the model (i.e. the predictor assignment to indices and
#'   initial index coefficients are generated randomly) are considered), and the
#'   final optimal model with the lowest loss is returned, and
#'   \code{"userInput"} - user specifies the initial model structure (i.e. the
#'   number of indices and the placement of index variables among indices) and
#'   the initial index coefficients through \code{index.ind} and
#'   \code{index.coefs} arguments respectively.
#' @param num_ind If \code{initialise = "ppr"} or \code{"multiple"}: an
#'   \code{integer} that specifies the number of indices to be used in the
#'   model(s). The default is \code{num_ind = 5}.
#' @param num_models If \code{initialise = "multiple"}: an \code{integer} that
#'   specifies the number of starting models to be checked. The default is
#'   \code{num_models = 5}.
#' @param seed If \code{initialise = "multiple"}: the seed to be set when
#'   generating random starting points.
#' @param index.ind If \code{initialise = "userInput"}: an \code{integer} vector
#'   that assigns group index for each predictor in \code{index.vars}.
#' @param index.coefs If \code{initialise = "userInput"}: a \code{numeric}
#'   vector of index coefficients.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted individually (rather than considering as
#'   part of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value to be used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for the objective function value (loss) of MIP.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @return A list with two elements: \item{initial}{A list of information of the
#'   model initialisation. (For descriptions of the list elements see
#'   \code{\link{make_smimodelFit}}).} \item{best}{A list of information of the
#'   final optimised model. (For descriptions of the list elements see
#'   \code{\link{make_smimodelFit}}).}
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


#' Constructor function for the class \code{smimodelFit}
#'
#' Constructs an object of class \code{smimodelFit} using the information passed
#' to arguments.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}).
#' @param yvar Name of the response variable as a character string.
#' @param neighbour `neighbour` argument passed from the outer function.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is "additive", where the initial model
#'   will be a nonparametric additive model. The other options are "linear" -
#'   linear regression model (i.e. a special case single-index model, where the
#'   initial values of the index coefficients are obtained through a linear
#'   regression), and "userInput" - user specifies the initial model structure
#'   (i.e. the number of indices and the placement of index variables among
#'   indices) and the initial index coefficients through \code{index.ind} and
#'   \code{index.coefs} arguments respectively.
#' @param index.ind If \code{initialise = "userInput"}: an \code{integer} vector
#'   that assigns group index for each predictor in \code{index.vars}.
#' @param index.coefs If \code{initialise = "userInput"}: a \code{numeric}
#'   vector of index coefficients.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted individually (rather than considering as
#'   part of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model.
#' @return  A list of initial model information. For descriptions of the list
#'   elements see \code{\link{make_smimodelFit}}).
new_smimodelFit <- function(data, yvar, neighbour = 0, 
                            family = gaussian(), index.vars, 
                            initialise = c("additive", "linear", "userInput"), 
                            index.ind = NULL, index.coefs = NULL, 
                            s.vars = NULL, linear.vars = NULL){
  stopifnot(tsibble::is_tsibble(data))
  data_index <- index(data)
  data_key <- key(data)[[1]]
  initialise <- match.arg(initialise)
  data <- data |>
    drop_na()
  Y_data <- as.matrix(data[ , yvar])
  X_index <- as.matrix(data[ , index.vars])
  if(initialise == "linear"){
    # Initialise the model with a linear model
    index.ind <- rep(1, length(index.vars))
    ind_pos <- split(seq_along(index.ind), index.ind)
    pre.formula <- lapply(index.vars, function(var) paste0(var)) |>
      paste(collapse = "+") 
    pre.formula <- paste(yvar, "~", pre.formula)
    if (!is.null(s.vars)){
      svars.formula <- lapply(s.vars, function(var) paste0(var)) |>
        paste(collapse = "+")
      pre.formula <- paste(pre.formula, "+", svars.formula)
    }
    if (!is.null(linear.vars)){
      linear.formula <- lapply(linear.vars, function(var) paste0(var)) |>
        paste(collapse = "+")
      pre.formula <- paste(pre.formula, "+", linear.formula)
    }
    pre.formula <- paste(pre.formula, "-", 1)
    fun1 <- mgcv::gam(as.formula(pre.formula), data = data, family = family,
                      method = "REML")
    add <- data |>
      drop_na() |>
      select({{ data_index }}, {{ data_key }})
    fun1$model <- bind_cols(add, fun1$model)
    fun1$model <- as_tsibble(fun1$model,
                             index = data_index,
                             key = all_of(data_key))
    # Index coefficients
    alpha <- index.coefs <- fun1$coefficients[1:length(index.vars)]
    alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    # Calculating indices
    ind <- vector(length = length(ind_pos), mode = "list")
    for(i in 1:length(ind)){
      if(length(ind_pos[[i]]) == 1){
        temp <- i
      }else{
        temp <- unlist(lapply(1:length(ind_pos[[i]]), function(x) paste0(i, x))) 
      }
      ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                               as.matrix(alpha[match(temp, names(alpha))], ncol = 1))
    }
    dat_names <- names(ind) <- paste0("index", 1:length(ind))
    dat <- tibble::as_tibble(ind)
    dat_new <- dplyr::bind_cols(data, dat)
  }else{
    if(initialise == "additive"){
      # Initialise the model with a nonparametric additive model
      index.ind <- seq_along(index.vars)
      # Index positions
      ind_pos <- split(seq_along(index.ind), index.ind)
      # Index coefficients
      alpha <- index.coefs <- rep(1, length(index.vars))
      alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    }else if(initialise == "userInput"){
      if (is.null(index.ind) | is.null(index.coefs)) stop("index.ind and/or index.coefs are/is not provided.")
      # Number of (index) predictors
      num_pred <- length(index.vars)
      # Index positions
      ind_pos <- split(seq_along(index.ind), index.ind)
      # Number of indices
      num_ind <- length(ind_pos)
      # Index coefficients
      alpha <- unlist(tapply(index.coefs, index.ind, normalise_alpha))
      # Constructing a new index coefficient vector to have all predictors in each index
      newIndex <- allpred_index(num_pred = num_pred,
                                num_ind = num_ind,
                                ind_pos = ind_pos,
                                alpha = alpha)
      alpha <- newIndex$alpha_init_new
      index.ind <- newIndex$index
      ind_pos <- newIndex$index_positions
      # Adjusting X (matrix of predictors) to fit number of indices
      X_index <- do.call(cbind, replicate(num_ind, X_index, simplify = FALSE))
      # Checking for all zero indices
      ind_rm_id <- numeric()
      ind_rm_pos <- numeric()
      for(i in 1:num_ind){
        if(all(alpha[ind_pos[[i]]] == 0)){
          ind_rm_id <- c(ind_rm_id, i)
          ind_rm_pos <- c(ind_rm_pos, ind_pos[[i]])
          message(paste0('Initial model', ': All coefficients of index', i,
                         ' are zero. Removing index', i, 
                         '. However, the variables in the removed index are considered in subsequent model searches.')) 
        }
      }
      if(length(ind_rm_id) != 0){
        X_index <- as.matrix(X_index[ , -ind_rm_pos])
        alpha <- alpha[-ind_rm_pos]
        num_ind  <- num_ind - length(ind_rm_id)
        index.ind <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          index.ind[[i]] <- rep(i, num_pred)
        }
        index.ind <- unlist(index.ind)
        ind_pos <- split(1:length(index.ind), index.ind)
      }
      names(alpha) <- NULL
      alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    }
    # Calculating indices
    ind <- vector(length = length(ind_pos), mode = "list")
    for(i in 1:length(ind)){
      if(length(ind_pos[[i]]) == 1){
        temp <- i
      }else{
        temp <- unlist(lapply(1:length(ind_pos[[i]]), function(x) paste0(i, x))) 
      }
      ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                               as.matrix(alpha[match(temp, names(alpha))], ncol = 1))
    }
    dat_names <- names(ind) <- paste0("index", 1:length(ind))
    dat <- tibble::as_tibble(ind)
    # Constructing the formula
    pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ', bs="cr")')) |>
      paste(collapse = "+") 
    pre.formula <- paste(yvar, "~", pre.formula)
    if (!is.null(s.vars)){
      svars.formula <- lapply(s.vars, function(var) paste0("s(", var, ', bs="cr")')) |>
        paste(collapse = "+") 
      pre.formula <- paste(pre.formula, "+", svars.formula)
    }
    if (!is.null(linear.vars)){
      linear.formula <- lapply(linear.vars, function(var) paste0(var)) |>
        paste(collapse = "+")
      pre.formula <- paste(pre.formula, "+", linear.formula)
    }
    # Model fitting
    dat_new <- dplyr::bind_cols(data, dat)
    dat_new <- dat_new |>
      tsibble::as_tsibble(index = {{data_index}}, key = {{data_key}})
    fun1 <- mgcv::gam(as.formula(pre.formula), data = dat_new, family = family,
                      method = "REML")
    add <- dat_new |>
      drop_na() |>
      select({{ data_index }}, {{ data_key }})
    fun1$model <- bind_cols(add, fun1$model)
    fun1$model <- as_tsibble(fun1$model,
                             index = data_index,
                             key = all_of(data_key))
  }
  smimodel <- make_smimodelFit(x = fun1, yvar = yvar, 
                               neighbour = neighbour,
                               index.vars = index.vars, 
                               index.ind = index.ind, index.data = dat_new,
                               index.names = dat_names, alpha = alpha,
                               s.vars = s.vars, linear.vars = linear.vars)
  return(smimodel)
}


#' Converting a fitted \code{gam} object to a \code{smimodelFit} object
#'
#' Converts a given object of class \code{gam} to an object of class
#' \code{smimodelFit}.
#'
#' @param x A fitted \code{gam} object.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour `neighbour` argument passed from the outer function.
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
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value to be used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for the objective function value (loss) of MIP.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
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
#'   fitting.} \item{gam}{Fitted \code{gam}.} \item{lambda0}{L0 penalty
#'   parameter used for model fitting.} \item{lambda2}{L2 penalty
#'   parameter used for model fitting.} \item{M}{Big-M value used in MIP.}
#'   \item{max.iter}{Maximum number of MIP iterations for a single round of
#'   index coefficients update.} \item{tol}{Tolerance for the objective function
#'    value (loss) used in solving MIP.} \item{tolCoefs}{Tolerance for
#'    coefficients used in updating index coefficients.} \item{TimeLimit}{Limit
#'    for the total time (in seconds) expended in a single MIP iteration.}
#'    \item{MIPGap}{Relative MIP optimality gap used.} \item{Nonconvex}{The
#'    strategy used for handling non-convex quadratic objectives or non-convex
#'    quadratic constraints in Gurobi solver.}
make_smimodelFit <- function(x, yvar, neighbour, index.vars, index.ind, index.data,
                             index.names, alpha, s.vars = NULL, linear.vars = NULL,
                             lambda0 = NULL, lambda2 = NULL, 
                             M = NULL, max.iter = NULL, 
                             tol = NULL, tolCoefs = NULL,
                             TimeLimit = NULL, 
                             MIPGap = NULL, NonConvex = NULL){
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
                   "neighbour" = neighbour, "gam" = x,
                   "lambda0" = lambda0, "lambda2" = lambda2,
                   "M" = M, "max.iter" = max.iter,
                   "tol" = tol, "tolCoefs" = tolCoefs,
                   "TimeLimit" = TimeLimit, "MIPGap" = MIPGap,
                   "Nonconvex" = NonConvex)
  class(smimodel) <- c("smimodelFit", "list")
  return(smimodel)
}


#' Updating a \code{smimodelFit}
#'
#' Optimises and updates a given \code{smimodelFit}.
#'
#' @param object A \code{smimodelFit} object.
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}).
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value to be used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for the objective function value (loss) of MIP.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @param ... Other arguments not currently used.
#' @return  A list of optimised model information. For descriptions of the list
#'   elements see \code{\link{make_smimodelFit}}).
update_smimodelFit <- function(object, data, lambda0 = 1, lambda2 = 1, 
                               M = 10, max.iter = 50, 
                               tol = 0.001, tolCoefs = 0.001,
                               TimeLimit = Inf, MIPGap = 1e-4, 
                               NonConvex = -1, verbose = FALSE, ...){
  if (!tsibble::is_tsibble(data)) stop("data is not a tsibble.")
  data_index <- index(data)
  data_key <- key(data)[[1]]
  data <- data |> drop_na()
  # Preparing inputs to `inner_update()`
  list_index <- object$alpha
  num_ind <- NCOL(list_index)
  alpha <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    alpha[[i]] <- list_index[ , i]
  }
  alpha <- unlist(alpha)
  names(alpha) <- NULL
  dgz <- as.matrix(object$derivatives)
  # Optimising the model
  best_alpha1 <- inner_update(x = object$gam, data = data, yvar = object$var_y,
                              family = object$gam$family$family,
                              index.vars = object$vars_index, 
                              s.vars = object$vars_s,
                              linear.vars = object$vars_linear, 
                              num_ind = num_ind, dgz = dgz, 
                              alpha_old = alpha, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, max.iter = max.iter, 
                              tol = tol, TimeLimit = TimeLimit, 
                              MIPGap = MIPGap, NonConvex = NonConvex, 
                              verbose = verbose)
  if(all(best_alpha1$best_alpha == 0)){
    # Constructing the formula and model fitting
    if (!is.null(object$vars_s) & !is.null(object$vars_linear)){
      pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) |>
        paste(collapse = "+") 
      pre.formula <- paste(object$var_y, "~", pre.formula)
      linear.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
        paste(collapse = "+")
      pre.formula <- paste(pre.formula, "+", linear.formula)
    }else if(!is.null(object$vars_s)){
      pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) |>
        paste(collapse = "+")
      pre.formula <- paste(object$var_y, "~", pre.formula)
    }else if(!is.null(object$vars_linear)){
      pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
        paste(collapse = "+")
      pre.formula <- paste(object$var_y, "~", pre.formula)
    }else{
      pre.formula <- paste(object$var_y, "~", 1)
    }
    fun1_final <- mgcv::gam(as.formula(pre.formula), data = data, 
                            family = object$gam$family$family, method = "REML")
    add <- data |>
      drop_na() |>
      select({{ data_index }}, {{ data_key }})
    fun1_final$model <- bind_cols(add, fun1_final$model)
    fun1_final$model <- as_tsibble(fun1_final$model,
                                   index = data_index,
                                   key = all_of(data_key))
    index.names <- paste0("index", 1:length(best_alpha1$ind_pos))
    print("Final model fitted!")
    final_smimodel <- make_smimodelFit(x = fun1_final, yvar = object$var_y, 
                                       neighbour = object$neighbour,
                                       index.vars = object$vars_index, 
                                       index.ind = best_alpha1$index.ind, 
                                       index.data = NULL, index.names = index.names,
                                       alpha = best_alpha1$best_alpha, 
                                       s.vars = object$vars_s,
                                       linear.vars = object$vars_linear,
                                       lambda0 = lambda0, lambda2 = lambda2, 
                                       M = M, max.iter = max.iter, 
                                       tol = tol, tolCoefs = tolCoefs,
                                       TimeLimit = TimeLimit, 
                                       MIPGap = MIPGap, NonConvex = NonConvex)
  }else{
    # Checking models with higher number of indices
    alpha_current <- best_alpha1$best_alpha
    loss_current <- best_alpha1$min_loss
    index_current <- best_alpha1$index.ind
    ind_pos_current <- best_alpha1$ind_pos
    X_new_current <- best_alpha1$X_new
    j <- length(ind_pos_current)+1
    num_pred <- length(object$vars_index)
    while(j <= num_pred){
      # Checking whether a new index can be added without variable repetition;
      # if not, terminating the loop.
      index_list <- split(alpha_current, index_current)
      index_mat <- do.call(rbind, index_list)
      drop_pred_ind <- which(colSums(index_mat) == 0)
      if(length(drop_pred_ind) == 0){
        break
      }else{
        # Initialising the new index to be added to the current model
        drop_pred_name <- object$vars_index[drop_pred_ind]
        X_init <- as.matrix(data[, drop_pred_name])
        index.ind <- rep(1, length(drop_pred_name))
        # Calculating indices
        ind <- vector(length = length(ind_pos_current), mode = "list")
        for(i in 1:length(ind)){
          ind[[i]] <- as.numeric(X_new_current[, ind_pos_current[[i]]] %*% 
                                   as.matrix(alpha_current[ind_pos_current[[i]]], ncol = 1))
        }
        dat_names <- names(ind) <- paste0("index", 1:length(ind))
        dat <- as_tibble(ind)
        # Nonlinear function update 
        # Constructing the formula
        pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) |>
          paste(collapse = "+") 
        pre.formula <- paste(object$var_y, "~", pre.formula)
        if (!is.null(object$vars_s)){
          svars.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ',bs="cr")')) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", svars.formula)
        }
        if (!is.null(object$vars_linear)){
          linear.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", linear.formula)
        }
        # Model fitting
        dat <- dplyr::bind_cols(data, dat)
        fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, 
                          family = object$gam$family$family, method = "REML")
        Yhat <- as.matrix(fun1$fitted.values, ncol = 1, nrow = length(fun1$fitted.values))
        # Residuals (R) matrix
        R <- as.matrix(data[ , object$var_y] - Yhat)
        coefs_init <- init_alpha(Y = R, X = X_init, index.ind = index.ind,
                                 init.type = "reg")$alpha_init
        new_ind <- numeric(length = num_pred)
        new_ind[drop_pred_ind] <- coefs_init
        alpha <- c(alpha_current, new_ind)
        X_index <- data[ , object$vars_index]
        # Adjusting X (matrix of predictors) to fit number of indices
        X_new <- as.matrix(do.call(cbind, replicate(j, X_index, simplify = FALSE)))
        index <- vector(mode = "list", length = j)
        for(i in 1:j){
          index[[i]] <- rep(i, num_pred)
        }
        index <- unlist(index)
        ind_pos <- split(1:length(index), index)
        # Calculating indices
        ind <- vector(length = length(ind_pos), mode = "list")
        for(i in 1:length(ind)){
          ind[[i]] <- as.numeric(X_new[, ind_pos[[i]]] %*% 
                                   as.matrix(alpha[ind_pos[[i]]], ncol = 1))
        }
        dat_names <- names(ind) <- paste0("index", 1:length(ind))
        dat <- as_tibble(ind)
        # Nonlinear function update 
        # Constructing the formula
        pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) |>
          paste(collapse = "+") 
        pre.formula <- paste(object$var_y, "~", pre.formula)
        if (!is.null(object$vars_s)){
          svars.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ',bs="cr")')) |>
            paste(collapse = "+") 
          pre.formula <- paste(pre.formula, "+", svars.formula)
        }
        if (!is.null(object$vars_linear)){
          linear.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", linear.formula)
        }
        # Model fitting
        dat <- dplyr::bind_cols(data, dat)
        gam2 <- mgcv::gam(as.formula(pre.formula), data = dat, 
                          family = object$gam$family$family, method = "REML")
        # Derivatives of the fitted smooths
        dgz <- vector(length = length(dat_names), mode = "list")
        for (i in seq_along(dat_names)) {
          temp <- gratia::derivatives(gam2, type = "central", data = dat, 
                                      select = paste0("s(", paste0(dat_names[i]), ")"))
          dgz[[i]] <- temp$.derivative
        }
        names(dgz) <- paste0("d", seq_along(dat_names))
        dgz <- as.matrix(as_tibble(dgz))
        # Optimising the new model
        best_alpha2 <- inner_update(x = gam2, data = data, yvar = object$var_y,
                                    family = object$gam$family$family,
                                    index.vars = object$vars_index, 
                                    s.vars = object$vars_s,
                                    linear.vars = object$vars_linear, 
                                    num_ind = j, dgz = dgz, 
                                    alpha_old = alpha, lambda0 = lambda0, 
                                    lambda2 = lambda2, M = M, max.iter = max.iter,
                                    tol = tol, TimeLimit = TimeLimit, 
                                    MIPGap = MIPGap, NonConvex = NonConvex,
                                    verbose = verbose)
        if(all(best_alpha2$best_alpha == 0)){
          # Constructing the formula and model fitting
          if (!is.null(object$vars_s) & !is.null(object$vars_linear)){
            pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) |>
              paste(collapse = "+")
            pre.formula <- paste(object$var_y, "~", pre.formula)
            linear.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
              paste(collapse = "+") 
            pre.formula <- paste(pre.formula, "+", linear.formula)
          }else if(!is.null(object$vars_s)){
            pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) |>
              paste(collapse = "+") 
            pre.formula <- paste(object$var_y, "~", pre.formula)
          }else if(!is.null(object$vars_linear)){
            pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
              paste(collapse = "+") 
            pre.formula <- paste(object$var_y, "~", pre.formula)
          }else{
            pre.formula <- paste(object$var_y, "~", 1)
          }
          fun_null <- mgcv::gam(as.formula(pre.formula), data = data, 
                                family = object$gam$family$family, 
                                method = "REML")
          add <- data |>
            drop_na() |>
            select({{ data_index }}, {{ data_key }})
          fun_null$model <- bind_cols(add, fun_null$model)
          fun_null$model <- as_tsibble(fun_null$model,
                                       index = data_index,
                                       key = all_of(data_key))
          Y <- as.matrix(data[ , object$var_y])
          Yhat <- as.matrix(fun_null$fitted.values, 
                            ncol = 1, nrow = length(fun_null$fitted.values))
          loss_null <- loss(Y = Y, Yhat = Yhat, alpha = best_alpha2$best_alpha,
                            lambda0 = lambda0, lambda2 = lambda2)
          if(loss_null >= loss_current){
            break
          }else{
            index.names <- paste0("index", 1:length(best_alpha2$ind_pos))
            alpha_current <- best_alpha2$best_alpha
            print("Final model fitted!")
            final_smimodel <- make_smimodelFit(x = fun_null, 
                                               yvar = object$var_y, 
                                               neighbour = object$neighbour,
                                               index.vars = object$vars_index, 
                                               index.ind = best_alpha2$index.ind, 
                                               index.data = NULL, 
                                               index.names = index.names,
                                               alpha = best_alpha2$best_alpha, 
                                               s.vars = object$vars_s,
                                               linear.vars = object$vars_linear,
                                               lambda0 = lambda0, lambda2 = lambda2, 
                                               M = M, max.iter = max.iter, 
                                               tol = tol, tolCoefs = tolCoefs,
                                               TimeLimit = TimeLimit, 
                                               MIPGap = MIPGap, NonConvex = NonConvex)
            break
          }
        }else{
          # Termination/continuation checks
          if((length(ind_pos_current) == length(best_alpha2$ind_pos))){
            comparison <- abs(alpha_current - best_alpha2$best_alpha) <= tolCoefs
            if(all(comparison)){
              if(best_alpha2$min_loss >= loss_current){
                break
              }else{
                alpha_current <- best_alpha2$best_alpha
                loss_current <- best_alpha2$min_loss
                index_current <- best_alpha2$index.ind
                ind_pos_current <- best_alpha2$ind_pos
                X_new_current <- best_alpha2$X_new
                break
              }
            }
            if(best_alpha2$min_loss >= loss_current){
              break
            }else{
              alpha_current <- best_alpha2$best_alpha
              loss_current <- best_alpha2$min_loss
              index_current <- best_alpha2$index.ind
              ind_pos_current <- best_alpha2$ind_pos
              X_new_current <- best_alpha2$X_new
            }
          }else{
            if(best_alpha2$min_loss >= loss_current){
              break
            }else{
              alpha_current <- best_alpha2$best_alpha
              loss_current <- best_alpha2$min_loss
              index_current <- best_alpha2$index.ind
              ind_pos_current <- best_alpha2$ind_pos
              X_new_current <- best_alpha2$X_new
            }
          }
          j <- length(ind_pos_current)+1
        }
      }
    }
    if(any(alpha_current != 0)){
      # Calculating indices
      ind <- vector(length = length(ind_pos_current), mode = "list")
      for(i in 1:length(ind)){
        ind[[i]] <- as.numeric(X_new_current[, ind_pos_current[[i]]] %*% 
                                 as.matrix(alpha_current[ind_pos_current[[i]]], ncol = 1))
      }
      dat_names <- names(ind) <- paste0("index", 1:length(ind))
      dat <- tibble::as_tibble(ind)
      ## Nonlinear function update 
      # Constructing the formula
      pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) |>
        paste(collapse = "+") 
      pre.formula <- paste(object$var_y, "~", pre.formula)
      if (!is.null(object$vars_s)){
        svars.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ',bs="cr")')) |>
          paste(collapse = "+") 
        pre.formula <- paste(pre.formula, "+", svars.formula)
      }
      if (!is.null(object$vars_linear)){
        linear.formula <- lapply(object$vars_linear, function(var) paste0(var)) |>
          paste(collapse = "+")
        pre.formula <- paste(pre.formula, "+", linear.formula)
      }
      dat <- dplyr::bind_cols(data, dat)
      # Model fitting
      fun1_final <- mgcv::gam(as.formula(pre.formula), data = dat, 
                              family = object$gam$family$family, method = "REML")
      add <- data |>
        drop_na() |>
        select({{ data_index }}, {{ data_key }})
      fun1_final$model <- bind_cols(add, fun1_final$model)
      fun1_final$model <- as_tsibble(fun1_final$model,
                                     index = data_index,
                                     key = all_of(data_key))
      print("Final model fitted!")
      final_smimodel <- make_smimodelFit(x = fun1_final, 
                                         yvar = object$var_y,
                                         neighbour = object$neighbour,
                                         index.vars = object$vars_index, 
                                         index.ind = index_current, 
                                         index.data = dat, 
                                         index.names = dat_names,
                                         alpha = alpha_current, 
                                         s.vars = object$vars_s,
                                         linear.vars = object$vars_linear,
                                         lambda0 = lambda0, lambda2 = lambda2, 
                                         M = M, max.iter = max.iter, 
                                         tol = tol, tolCoefs = tolCoefs,
                                         TimeLimit = TimeLimit, 
                                         MIPGap = MIPGap, NonConvex = NonConvex)
    }
  }
  return(final_smimodel)
}


#' Updating index coefficients and non-linear functions iteratively
#'
#' Iteratively updates index coefficients and non-linear functions using mixed
#' integer programming. (A helper function used within
#' \code{\link{update_smimodelFit}}; users are not expected to directly call
#' this function.)
#'
#' @param x Fitted \code{gam}.
#' @param data Training data set on which models will be trained. Should be a
#'   \code{tsibble}.
#' @param yvar Name of the response variable as a character string.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices should be estimated.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted individually (rather than considering as
#'   part of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model.
#' @param num_ind Number of indices.
#' @param dgz The \code{tibble} of derivatives of the estimated smooths.
#' @param alpha_old Current vector of index coefficients.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value to be used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @return A list containing following elements: \item{best_alpha}{The vector of
#'   best index coefficient estimates.} \item{min_loss}{Minimum value of the
#'   objective function(loss).}
#' \item{index.ind}{An \code{integer} vector that assigns group index for each
#' predictor, corresponding to \code{best_alpha}.}
#' \item{ind_pos}{A list that indicates which predictors belong to which index,
#' corresponding to \code{best_alpha}.}
#' \item{X_new}{A matrix of selected predictor variables, corresponding to
#' \code{best_alpha}.}
inner_update <- function(x, data, yvar, family = gaussian(), index.vars, 
                         s.vars, linear.vars, num_ind, dgz, alpha_old, 
                         lambda0 = 1, lambda2 = 1, M = 10, max.iter = 50, 
                         tol = 0.001, TimeLimit = Inf,
                         MIPGap = 1e-4, NonConvex = -1, verbose = FALSE){
  data <- data |>
    drop_na()
  data.Y <- as.matrix(data[ , yvar])
  X_index <- as.matrix(data[ , index.vars])
  num_pred <- NCOL(X_index)
  index <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    index[[i]] <- rep(i, num_pred)
  }
  index <- unlist(index)
  ind_pos <- split(1:length(index), index)
  Yhat <- as.matrix(x$fitted.values, ncol = 1)
  # Residuals (R) matrix
  R <- as.matrix(data.Y - Yhat)
  # Loss
  l2 <- loss(Y = as.matrix(data.Y), Yhat = Yhat, 
             alpha = alpha_old, 
             lambda0 = lambda0, lambda2 = lambda2)
  # Setting up a counter for increases in loss (If the loss increases for 
  # 3 consecutive iterations, the algorithm will terminate.)
  increase_count <- 0
  # Adjusting X (matrix of predictors) to fit number of indices
  X_new <- do.call(cbind, replicate(num_ind, X_index, simplify = FALSE))
  # Initial minimum loss and best index coefficient estimates
  best_l2 <- l2
  best_alpha <- alpha_old
  best_index <- index
  best_ind_pos <- ind_pos
  best_X_new <- X_new
  # Iteratively update index coefficients
  maxIt <- 1
  while(maxIt <= max.iter){
    alpha_new <- update_alpha(Y = R, X = X_new, num_pred = num_pred, 
                              num_ind = num_ind, index.ind = index, dgz = dgz, 
                              alpha_old = alpha_old, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, TimeLimit = TimeLimit,
                              MIPGap = MIPGap, NonConvex = NonConvex, 
                              verbose = verbose)
    if(all(alpha_new == 0)){
      best_l2 <- NULL
      best_alpha <- alpha_new
      best_index <- index
      best_ind_pos <- ind_pos
      best_X_new <- NULL
      print("Null indices are produced!")
      break
    }else{
      # Checking for all zero indices
      ind_rm_id <- numeric()
      ind_rm_pos <- numeric()
      for(i in 1:num_ind){
        if(all(alpha_new[ind_pos[[i]]] == 0)){
          ind_rm_id <- c(ind_rm_id, i)
          ind_rm_pos <- c(ind_rm_pos, ind_pos[[i]])
          message(paste0('Iteration ', maxIt, ': All coefficients of index', i,
                         ' are zero. Removing index', i, ' from subsequent iterations.')) 
        }
      }
      if(length(ind_rm_id) != 0){
        X_new <- as.matrix(X_new[ , -ind_rm_pos])
        alpha_new <- alpha_new[-ind_rm_pos]
        num_ind  <- num_ind - length(ind_rm_id)
        index <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          index[[i]] <- rep(i, num_pred)
        }
        index <- unlist(index)
        ind_pos <- split(1:length(index), index)
      }
      # Calculating indices
      ind <- vector(length = length(ind_pos), mode = "list")
      for(i in 1:length(ind)){
        ind[[i]] <- as.numeric(X_new[, ind_pos[[i]]] %*% 
                                 as.matrix(alpha_new[ind_pos[[i]]], ncol = 1))
      }
      dat_names <- names(ind) <- paste0("index", 1:length(ind))
      dat <- as_tibble(ind)
      # Nonlinear function update 
      # Constructing the formula
      pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ', bs="cr")')) |>
        paste(collapse = "+") 
      pre.formula <- paste(yvar, "~", pre.formula)
      if (!is.null(s.vars)){
        svars.formula <- lapply(s.vars, function(var) paste0("s(", var, ', bs="cr")')) |>
          paste(collapse = "+") 
        pre.formula <- paste(pre.formula, "+", svars.formula)
      }
      if (!is.null(linear.vars)){
        linear.formula <- lapply(linear.vars, function(var) paste0(var)) |>
          paste(collapse = "+")
        pre.formula <- paste(pre.formula, "+", linear.formula)
      }
      # Model fitting
      dat <- dplyr::bind_cols(data, dat)
      fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, family = family, 
                        method = "REML")
      Yhat <- as.matrix(fun1$fitted.values, ncol = 1)
      # Derivatives of the fitted smooths
      dgz <- vector(length = length(dat_names), mode = "list")
      for (i in seq_along(dat_names)) {
        temp <- gratia::derivatives(fun1, type = "central", data = dat, 
                                    select = paste0("s(", dat_names[i], ")"))
        dgz[[i]] <- temp$.derivative
      }
      names(dgz) <- paste0("d", seq_along(dat_names))
      dgz <- as.matrix(as_tibble(dgz))
      # Residuals (R) matrix
      R <- as.matrix(data.Y - Yhat)
      # Loss
      l2_new <- loss(Y = as.matrix(data.Y), Yhat = Yhat, 
                     alpha = alpha_new, 
                     lambda0 = lambda0, lambda2 = lambda2)
      eps <- (l2 - l2_new)/l2
      # Update loss and estimates
      alpha_old <- alpha_new
      l2 <- l2_new
      # Update minimum loss and best estimates
      if(l2_new < best_l2){
        best_l2 <- l2_new
        best_alpha <- alpha_new
        best_index <- index
        best_ind_pos <- ind_pos
        best_X_new <- X_new
      }
      # Check tolerance conditions
      if(eps <= 0){
        increase_count <- increase_count + 1
      }else{
        increase_count <- 0
      }
      if ((eps > 0) & (eps < tol)) { 
        print("Tolerance for loss reached!")
        break
      }else if(increase_count >= 3){
        print("Loss increased for 3 consecutive iterations!")
        break
      }
      maxIt <- maxIt + 1
      if (maxIt > max.iter) { 
        print("Maximum iterations reached!")
      }
    }
  }
  output <- list("best_alpha" = best_alpha, "min_loss" = best_l2, 
                 "index.ind" = best_index, "ind_pos" = best_ind_pos, 
                 "X_new" = best_X_new)
  return(output)
}


#' Updating index coefficients using MIP
#'
#' Updates index coefficients by solving a mixed integer program.
#'
#' @param Y Column matrix of response.
#' @param X Matrix of predictors (size adjusted to number of indices).
#' @param num_pred Number of predictors.
#' @param num_ind Number of indices.
#' @param index.ind An integer vector that assigns group index for each
#'   predictor.
#' @param dgz The \code{tibble} of derivatives of the estimated smooths from
#'   previous iteration.
#' @param alpha_old Vector of index coefficients from previous iteration.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value to be used in MIP.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @return A vector of normalised index coefficients.
update_alpha <- function(Y, X, num_pred, num_ind, index.ind, dgz, alpha_old, 
                         lambda0 = 1, lambda2 = 1, M = 10, TimeLimit = Inf,
                         MIPGap = 1e-4, NonConvex = -1, verbose = FALSE){
  dgz <- as.matrix(dgz[, index.ind])
  V <- X * dgz
  Vt <- t(V)
  VtV <- Vt %*% V
  VtR <- Vt %*% Y
  VtR_t <- t(VtR)
  VtV_alpha_old_t <- t(VtV %*% alpha_old)
  # Converting the sparse matrix to a normal matrix in R because `ROI::Q_objective()`
  # does not accept sparse matrices.
  Q <- as.matrix(2 * Matrix::bdiag(VtV + lambda2*diag(num_pred*num_ind),
                                   diag(0, num_pred*num_ind)))
  L <- t(as.matrix(c((-2*VtV_alpha_old_t - 2*VtR_t),
                     rep(lambda0, num_ind*num_pred)), 
                   ncol = 2*num_ind*num_pred, nrow = 1))
  # Objective function
  obj <- ROI::Q_objective(Q = Q, L = L)
  # Constraints
  C1 <- Matrix::bdiag(replicate(num_ind, diag(num_pred), simplify = FALSE))
  C2 <- Matrix::bdiag(replicate(num_ind, diag(-M, num_pred), simplify = FALSE))
  cons1 <- as.matrix(cbind(C1, C2))
  C3 <- Matrix::bdiag(replicate(num_ind, diag(-1, num_pred), simplify = FALSE))
  C4 <- Matrix::bdiag(replicate(num_ind, diag(-M, num_pred), simplify = FALSE))
  cons2 <- as.matrix(cbind(C3, C4))
  if(num_ind > 1){
    C5 <- do.call(cbind, replicate(num_ind, diag(0, num_pred), simplify = FALSE))
    C6 <- do.call(cbind, replicate(num_ind, diag(num_pred), simplify = FALSE))
    cons3 <- as.matrix(cbind(C5, C6))
    cons <- ROI::L_constraint(
      L = rbind(cons1, cons2, cons3),
      dir = c(rep("<=", (2*num_ind*num_pred)), rep("<=", num_pred)),
      rhs = c(rep(0, (2*num_ind*num_pred)), rep(1, num_pred))
    )
  }else{
    cons <- ROI::L_constraint(
      L = rbind(cons1, cons2),
      dir = rep("<=", (2*num_ind*num_pred)),
      rhs = rep(0, (2*num_ind*num_pred))
    )
  }
  # Optimization problem
  init <- ROI::OP(obj, cons, types = c(rep("C", num_ind*num_pred), rep("B", num_ind*num_pred)),
                  bounds = ROI::V_bound(li = 1:(num_ind*num_pred), lb = rep(-M, num_ind*num_pred),
                                        ui = 1:(num_ind*num_pred), ub = rep(M, num_ind*num_pred),
                                        nobj = 2*num_ind*num_pred), # lower default bound is 0
                  maximum = FALSE)
  # Solving
  sol <- ROI::ROI_solve(init, solver = "gurobi", 
                        TimeLimit = TimeLimit, 
                        MIPGap = MIPGap, 
                        NonConvex = NonConvex, 
                        verbose = verbose)
  # Binary variables check
  coefs <- sol$solution[1:(num_ind*num_pred)]
  indicators <- sol$solution[((num_ind*num_pred)+1):(num_ind*num_pred*2)]
  zero_pos <- which(indicators < 1e-5)
  coefs[zero_pos] <- 0
  # Index coefficients
  alpha <- unlist(tapply(coefs, index.ind, normalise_alpha))
  return(alpha)
}


#' Initialising index coefficients
#'
#' Initialises index coefficient vector through linear regression or penalised
#' linear regression.
#'
#' @param Y Column matrix of response.
#' @param X Matrix of predictors entering indices.
#' @param index.ind An \code{integer} vector that assigns group index for each
#'   predictor.
#' @param init.type Type of initialisation for index coefficients.
#'   (\code{"penalisedReg"} - Penalised linear regression; \code{"reg"} - Linear
#'   regression)
#' @param lambda0 If \code{init.type = "penalisedReg"}, penalty parameter for L0
#'   penalty.
#' @param lambda2 If \code{init.type = "penalisedReg"}, penalty parameter for L2
#'   penalty.
#' @param M If \code{init.type = "penalisedReg"}, the big-M value to be used in
#'   the MIP.
#' @return A list containing the following components:
#'   \item{alpha_init}{Normalised vector of index coefficients.}
#' \item{alpha_nonNormalised}{Non-normalised (i.e. prior to normalising)
#' vector of index coefficients.}
init_alpha <- function(Y, X, index.ind, init.type = "penalisedReg", 
                       lambda0 = 1, lambda2 = 1, M = 10){
  p = NCOL(X)
  Xt = t(X)
  XtX = Xt %*% X
  XtY = Xt %*% Y
  if(init.type == "reg"){
    Q <- 2 * XtX
    L <- t(as.matrix(-2 * t(XtY), nrow = 1, ncol = p))
    # Objective function
    obj = ROI::Q_objective(Q = Q, L = L)
    # Optimization problem
    init <- ROI::OP(obj, types = rep("C", p),
                    bounds = V_bound(li = 1:p, lb = rep(-Inf, p)),
                    maximum = FALSE)
    # Solving
    sol <- ROI::ROI_solve(init, solver = "gurobi")
    # Initial alpha values
    alpha_raw <- sol$solution[1:p]
    alpha_0 <- unlist(tapply(sol$solution[1:p], index.ind, normalise_alpha))
  }else if(init.type == "penalisedReg"){
    Q <- as.matrix(2 * Matrix::bdiag(XtX + (lambda2*diag(p)), diag(0, p))) 
    # Converting the sparse matrix to a normal matrix in R because `ROI::Q_objective()`
    # does not accept sparse matrices.
    L <- t(as.matrix(c(-2 * t(XtY), rep(lambda0, p)), nrow = 1, ncol = 2*p))
    # Objective function
    obj = ROI::Q_objective(Q = Q, L = L)
    # Constraints
    cons1 <- cbind(diag(p), - M * diag(p))
    cons2 <- cbind(- diag(p), - M * diag(p))
    cons = ROI::L_constraint(
      L = rbind(cons1, cons2),
      dir = rep("<=", (2*p)),
      rhs = rep(0, 2*p)
    )
    # Optimization problem
    init <- ROI::OP(obj, cons, types = c(rep("C", p), rep("B", p)), 
                    bounds = ROI::V_bound(li = 1:p, lb = rep(- M, p),
                                          ui = 1:p, ub = rep(M, p),
                                          nobj = 2*p), # lower default bound is 0
                    maximum = FALSE)
    # Solving
    sol <- ROI::ROI_solve(init, solver = "gurobi")
    # Binary variables check
    coefs <- sol$solution[1:p]
    indicators <- sol$solution[(p+1):(p*2)]
    zero_pos <- which(indicators < 1e-5)
    coefs[zero_pos] <- 0
    # Initial alpha values
    alpha_raw <- coefs
    alpha_0 <- unlist(tapply(coefs, index.ind, normalise_alpha))
  }
  output <- list("alpha_init" = alpha_0, "alpha_nonNormalised" = alpha_raw)
  return(output)
}


#' Calculating the loss of the MIP used to estimate a SMI model
#'
#' Calculates the value of the objective function (loss function) of the mixed
#' integer program used to estimate a SMI model.
#'
#' @param Y Column matrix of response.
#' @param Yhat Predicted value of the response.
#' @param alpha Vector of index coefficients.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @return A \code{numeric}.
loss <- function(Y, Yhat, alpha, lambda0, lambda2){
  sqdError <- sum((Y - Yhat)^2)
  L0Penalty <- lambda0*sum(alpha != 0)
  nonzeroAlpha <- alpha[alpha != 0]
  L2Penalty <- lambda2*sum(nonzeroAlpha^2)
  return((sqdError + L0Penalty + L2Penalty))
}


#' Scaling index coefficient vectors to have unit norm
#'
#' Scales a coefficient vector of a particular index to have unit norm.
#'
#' @param alpha A vector of index coefficients.
#' @return A \code{numeric} vector.
normalise_alpha <- function (alpha) {
  anorm <- norm(matrix(alpha, ncol = 1))
  if (!(is.na(anorm) | all(alpha == 0))) 
    alpha <- alpha/anorm
  return(alpha)
}


#' Scale data
#'
#' Scales the columns of the \code{data} corresponding to \code{index.vars}.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   \code{tibble}.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @return A list containing the following components: \item{scaled_data}{The
#'   scaled data set of class \code{tibble}.} \item{scaled_info}{A named
#'   \code{numeric} vector of standard deviations of \code{index.vars} that were
#'   used to scale the corresponding columns of \code{data}.}
scaling <- function(data, index.vars){
  scaleData <- scale(data[ , index.vars], center = FALSE, 
                     scale = apply(data[ , index.vars], 2, sd, na.rm = TRUE))
  scaleInfo <- attributes(scaleData)$`scaled:scale`
  data <- data |>
    dplyr::select(-{{index.vars}})
  data <- dplyr::bind_cols(data, scaleData)
  output <- list("scaled_data" = data, "scaled_info" = scaleInfo)
  return(output)
}


#' Unscale a fitted \code{smimodel}
#'
#' Transforms back the index coefficients to suit original-scale index variables
#' if the same were scaled when estimating the \code{smimodel} (happens in
#' \code{initialise = "ppr"} in \code{\link{model_smimodel}} or
#' \code{\link{greedy_smimodel}}). Users are not expected to directly use this
#' function; usually called within \code{\link{smimodel.fit}}.
#'
#' @param object A \code{smimodel} object.
#' @param scaledInfo The list returned from a call of the function
#'   \code{\link{scaling}}.
#' @return A \code{smimodel} object.
unscaling <- function(object, scaledInfo){
  scaledInfo <- scaledInfo$scaled_info
  list_index <- object$alpha
  for(b in 1:NCOL(list_index)){
    temp_coef <- list_index[ , b]/scaledInfo
    object$alpha[ , b] <- temp_coef
  }
  return(object)
}


#' Constructing index coefficient vectors with all predictors in each index
#'
#' Constructs vectors of coefficients for each index including a coefficient for
#' all the predictors that are entering indices. i.e. if a coefficient is not
#' provided for a particular predictor in a particular index, the function will
#' replace the missing coefficient with a zero.
#'
#' @param num_pred Number of predictors.
#' @param num_ind Number of indices.
#' @param ind_pos A list of length = \code{num_ind} that indicates which
#'   predictors belong to which index.
#' @param alpha A vector of index coefficients.
#' @return A list containing the following components: \item{alpha_init_new}{A
#'   \code{numeric} vector of index coefficients.}
#'   \item{index}{An \code{integer} vector that assigns group indices for each
#'   predictor.} \item{index_positions}{A list of length = \code{num_ind} that
#'   indicates which predictors belong to which index.}
allpred_index <- function(num_pred, num_ind, ind_pos, alpha){
  init_list <- vector(mode = "list", length = num_ind)
  index <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    if(length(ind_pos[[i]]) == 1){
      temp <- i
    }else{
      temp <- unlist(lapply(1:length(ind_pos[[i]]), function(x) paste0(i, x))) 
    }
    if(length(alpha) == (num_ind*num_pred)){
      init_list[[i]] <- alpha[match(temp, names(alpha))]
    }else{
      init_list[[i]] <- numeric(length = num_pred)
      init_list[[i]][ind_pos[[i]]] <- alpha[match(temp, names(alpha))]
    }
    index[[i]] <- rep(i, num_pred)
  }
  alpha_init <- unlist(init_list)
  index <- unlist(index)
  ind_pos <- split(1:length(index), index)
  output <- list("alpha_init_new" = alpha_init,
                 "index" = index,
                 "index_positions" = ind_pos)
  return(output)
}


#' Splitting predictors into multiple indices
#'
#' Splits a given number of predictors into a given number of indices.
#'
#' @param num_pred Number of predictors.
#' @param num_ind Number of indices.
#' @return A list containing the following components: \item{index}{An
#'   \code{integer} vector that assigns group indices for each predictor.}
#'   \item{index_positions}{A list of length = \code{num_ind} that indicates
#'   which predictors belong to which index.}
split_index <- function(num_pred, num_ind){
  split_list <- numeric()
  rem_pred <- num_pred
  rem_ind <- num_ind
  count = 1
  while(count <= num_ind){
    split_num <- ceiling(rem_pred/rem_ind)
    rem_pred <- rem_pred - split_num
    rem_ind <- rem_ind - 1
    split1 <- rep(count, split_num)
    split_list <- c(split_list, split1)
    count <- count + 1
  }
  split_pos <- split(1:length(split_list), split_list) 
  output <- list("index" = split_list, "index_positions" = split_pos)
  return(output)
}