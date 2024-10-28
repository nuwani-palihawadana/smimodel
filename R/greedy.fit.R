#' Greedy search for tuning penalty parameters
#'
#' Function to perform a greedy search over a given grid of penalty parameter
#' combinations (lambda0, lambda2), and fits a single SMI model with the best
#' (lowest validation set MSE) penalty parameter combination. This is a helper
#' function designed to be called from \code{\link{greedy_smimodel}}.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}). If
#'   multiple models are fitted, the grouping variable should be the \code{key}
#'   of the \code{tsibble}. If a \code{key} is not specified, a dummy key with
#'   only one level will be created.
#' @param val.data Validation data set. (The data set on which the penalty
#'   parameter selection will be performed.) Must be a data set of class
#'   \code{tsibble}. (Once the penalty parameter selection is completed, the
#'   best model will be re-fitted for the combined data set \code{data +
#'   val.data}.)
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
#' @param lambda0_seq A \code{numeric} vector of candidate values for lambda0
#'   (penalty parameter for L0 penalty).
#' @param lambda2_seq A \code{numeric} vector of candidate values for lambda2
#'   (penalty parameter for L2 penalty).
#' @param lambda_step Step size between two adjacent values in
#'   \code{lambda0_seq} and \code{lambda2_seq}.
#' @param lambda0_start_seq A subset from \code{lambda0_seq} as candidate
#'   starting points for the greedy search.
#' @param lambda2_start_seq A subset from \code{lambda2_seq} as candidate
#'   starting points for the greedy search.
#' @param refit Whether to refit the model combining training and validation
#'   sets after parameter tuning. If \code{FALSE}, the final model will be
#'   estimated only on the training set.
#' @param M Big-M value used in MIP.
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
#' @param parallel The option to use parallel processing in fitting SMI models
#'   for different penalty parameter combinations.
#' @param workers If \code{parallel = TRUE}: Number of cores to use.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{val.data} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{val.data}, with no break in the lagged variable sequence
#'   even if some of the intermediate lags are not used as predictors.
#' @return  A list that contains three elements: \item{initial}{A list of
#'   information of the model initialisation. (For descriptions of the list
#'   elements see \code{\link{make_smimodelFit}}).} \item{best}{A list of
#'   information of the final optimised model. (For descriptions of the list
#'   elements see \code{\link{make_smimodelFit}}).} \item{best_lambdas}{Selected
#'   penalty parameter combination.}

greedy.fit <- function(data, val.data, yvar, neighbour = 0, 
                       family = gaussian(), index.vars, 
                       initialise = c("ppr", "additive", "linear", 
                                      "multiple", "userInput"),
                       num_ind = 5, num_models = 5, seed = 123, index.ind = NULL, 
                       index.coefs = NULL, s.vars = NULL, linear.vars = NULL, 
                       lambda0_seq, lambda2_seq, lambda_step,
                       lambda0_start_seq, lambda2_start_seq, refit = TRUE,
                       M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                       TimeLimit = Inf, MIPGap = 1e-4, NonConvex = -1, 
                       verbose = FALSE, parallel = FALSE, workers = NULL,
                       recursive = FALSE, recursive_colRange = NULL){
  # Full grid
  grid1 <- expand.grid(lambda0_seq, lambda2_seq)
  # Data frame for storing all combinations searched
  all_comb <- data.frame()
  # Current minimum MSE
  current_MSE <- Inf
  # Current best lambdas
  current_lambdas <- numeric(length = 2)
  # Select map function
  if(parallel){
    future::plan("future::multisession", workers = workers)
    map_f <- furrr::future_map
  } else {
    map_f <- purrr::map
  }
  # Starting point options
  lambda_comb <- expand.grid(lambda0_start_seq, lambda2_start_seq)
  # Model fitting for each combination of lambdas
  MSE_list <- seq(1, NROW(lambda_comb), by = 1) |>
    map_f(~ tune_smimodel(data = data, val.data = val.data, yvar = yvar, 
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
                          lambda.comb = as.numeric(lambda_comb[., ]),
                          M = M, max.iter = max.iter, 
                          tol = tol, tolCoefs = tolCoefs,
                          TimeLimit = TimeLimit, MIPGap = MIPGap,
                          NonConvex = NonConvex, verbose = verbose,
                          recursive = recursive,
                          recursive_colRange = recursive_colRange))
  #if (parallel) future:::ClusterRegistry("stop")
  # Selecting best starting point
  min_lambda_pos <- which.min(unlist(MSE_list))
  min_MSE <- min(unlist(MSE_list))
  min_lambdas <- as.numeric(lambda_comb[min_lambda_pos, ])
  print("First round completed; starting point selected!")
  # Updating searched combinations store
  all_comb <- bind_rows(all_comb, lambda_comb)
  # Greedy search
  while(min_MSE < current_MSE){
    current_MSE <- min_MSE
    current_lambdas <- min_lambdas 
    # Constructing new search space
    lambda0_seq_new <- c((current_lambdas[1] - lambda_step), current_lambdas[1], 
                         (current_lambdas[1] + lambda_step))
    lambda2_seq_new <- c((current_lambdas[2] - lambda_step), current_lambdas[2], 
                         (current_lambdas[2] + lambda_step))
    lambda_comb_new <- expand.grid(lambda0_seq_new, lambda2_seq_new)
    lambda_exist1 <- do.call(paste0, lambda_comb_new) %in% do.call(paste0, grid1)
    lambda_comb_new1 <- lambda_comb_new[lambda_exist1 == TRUE, ]
    lambda_exist2 <- do.call(paste0, lambda_comb_new1) %in% do.call(paste0, all_comb)
    lambda_comb <- lambda_comb_new1[lambda_exist2 == FALSE, ]
    if(NROW(lambda_comb) == 0){
      break 
    }else{
      MSE_list <- seq(1, NROW(lambda_comb), by = 1) |>
        map_f(~ tune_smimodel(data = data, val.data = val.data, yvar = yvar, 
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
                              lambda.comb = as.numeric(lambda_comb[., ]),
                              M = M, max.iter = max.iter, 
                              tol = tol, tolCoefs = tolCoefs,
                              TimeLimit = TimeLimit, MIPGap = MIPGap,
                              NonConvex = NonConvex, verbose = verbose,
                              recursive = recursive,
                              recursive_colRange = recursive_colRange))
      # Selecting best starting point
      min_lambda_pos <- which.min(unlist(MSE_list))
      min_MSE <- min(unlist(MSE_list))
      min_lambdas <- as.numeric(lambda_comb[min_lambda_pos, ])
      print("Another round completed!")
      # Updating searched combinations store
      all_comb <- bind_rows(all_comb, lambda_comb)
    }
  }
  # Final smimodel: re-fit the best model for the combined data set training + validation
  # Data
  if(refit == TRUE){
    combinedData <- dplyr::bind_rows(data, val.data)
    final_smimodel_list <- smimodel.fit(data = combinedData, yvar = yvar, 
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
                                        lambda0 = current_lambdas[1], 
                                        lambda2 = current_lambdas[2],
                                        M = M, max.iter = max.iter, 
                                        tol = tol, tolCoefs = tolCoefs,
                                        TimeLimit = TimeLimit, MIPGap = MIPGap,
                                        NonConvex = NonConvex, verbose = verbose)
  }else{
    final_smimodel_list <- smimodel.fit(data = data, yvar = yvar, 
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
                                        lambda0 = current_lambdas[1], 
                                        lambda2 = current_lambdas[2],
                                        M = M, max.iter = max.iter, 
                                        tol = tol, tolCoefs = tolCoefs,
                                        TimeLimit = TimeLimit, MIPGap = MIPGap,
                                        NonConvex = NonConvex, verbose = verbose)
  }
  print("Final model fitted!")
  output <- list("initial" = final_smimodel_list$initial, 
                 "best" = final_smimodel_list$best,
                 "best_lambdas" = current_lambdas)
  return(output)
}
