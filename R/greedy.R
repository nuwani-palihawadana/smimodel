#' Greedy search for tuning penalty parameters
#'
#' Performs a greedy search over a given grid of penalty parameter combinations
#' (lambda0, lambda2), and fits a SMI model with the best (lowest in-sample MSE)
#' penalty parameter combination.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param yvar Name of the response variable as a character string.
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
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#' @param lambda0_seq A numeric vector of candidate values for lambda0 (penalty
#'   parameter for L0 penalty).
#' @param lambda2_seq A numeric vector of candidate values for lambda2 (penalty
#'   parameter for L2 penalty).
#' @param lambda_step Step size between two adjacent values in `lambda0_seq` and
#'   `lambda2_seq`.
#' @param lambda0_start_seq A subset from `lambda0_seq` as candidate starting
#'   points for the greedy search.
#' @param lambda2_start_seq A subset from `lambda2_seq` as candidate starting
#'   points for the greedy search.
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
#' @param parallel The option to use parallel processing in fitting `smimodel`s
#'   for different penalty parameter combinations.
#' @param workers If `parallel = TRUE`: Number of cores to use.
#'
#' @importFrom future plan
#' @importFrom furrr future_map
#' @importFrom purrr map
#' @importFrom stats gaussian
#'
#' @export
greedy <- function(data, yvar, family = gaussian(), index.vars, 
                   initialise = c("ppr", "additive", "linear", 
                                  "multiple", "userInput"),
                   num_ind = 5, num_models = 5, seed = 123, index.ind = NULL, 
                   index.coefs = NULL, linear.vars = NULL, 
                   lambda0_seq, lambda2_seq, lambda_step,
                   lambda0_start_seq, lambda2_start_seq, 
                   M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                   TimeLimit = Inf, MIPGap = 1e-4, NonConvex = -1, 
                   verbose = FALSE, parallel = FALSE, workers = NULL){
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
  MSE_list <- seq(1, NROW(lambda_comb), by = 1) %>%
    map_f(~ smimodel_tune(data = data, yvar = yvar, 
                          family = family,
                          index.vars = index.vars, 
                          initialise = initialise, 
                          num_ind = num_ind, num_models = num_models, 
                          seed = seed,
                          index.ind = index.ind, 
                          index.coefs = index.coefs,
                          linear.vars = linear.vars,
                          lambda.comb = as.numeric(lambda_comb[., ]),
                          M = M, max.iter = max.iter, 
                          tol = tol, tolCoefs = tolCoefs,
                          TimeLimit = TimeLimit, MIPGap = MIPGap,
                          NonConvex = NonConvex, verbose = verbose))
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
      MSE_list <- seq(1, NROW(lambda_comb), by = 1) %>%
        map_f(~ smimodel_tune(data = data, yvar = yvar, 
                              family = family,
                              index.vars = index.vars, 
                              initialise = initialise, 
                              num_ind = num_ind, num_models = num_models, 
                              seed = seed,
                              index.ind = index.ind, 
                              index.coefs = index.coefs,
                              linear.vars = linear.vars,
                              lambda.comb = as.numeric(lambda_comb[., ]),
                              M = M, max.iter = max.iter, 
                              tol = tol, tolCoefs = tolCoefs,
                              TimeLimit = TimeLimit, MIPGap = MIPGap,
                              NonConvex = NonConvex, verbose = verbose))
      # Selecting best starting point
      min_lambda_pos <- which.min(unlist(MSE_list))
      min_MSE <- min(unlist(MSE_list))
      min_lambdas <- as.numeric(lambda_comb[min_lambda_pos, ])
      print("Another round completed!")
      # Updating searched combinations store
      all_comb <- bind_rows(all_comb, lambda_comb)
    }
  }
  # Final smimodel
  final_smimodel_list <- smimodel(data = data, yvar = yvar, 
                                  family = family,
                                  index.vars = index.vars, 
                                  initialise = initialise, 
                                  num_ind = num_ind, num_models = num_models, 
                                  seed = seed,
                                  index.ind = index.ind, 
                                  index.coefs = index.coefs,
                                  linear.vars = linear.vars,
                                  lambda0 = current_lambdas[1], 
                                  lambda2 = current_lambdas[2],
                                  M = M, max.iter = max.iter, 
                                  tol = tol, tolCoefs = tolCoefs,
                                  TimeLimit = TimeLimit, MIPGap = MIPGap,
                                  NonConvex = NonConvex, verbose = verbose)
  print("Final model fitted!")
  output <- list("initial" = final_smimodel_list$initial, 
                 "best" = final_smimodel_list$best,
                 "best_lambdas" = current_lambdas)
  return(output)
}
