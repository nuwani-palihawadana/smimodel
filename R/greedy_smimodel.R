#' SMI model estimation through a greedy search for penalty parameters
#'
#' Performs a greedy search over a given grid of penalty parameter combinations
#' (lambda0, lambda2), and fits SMI model(s) with best (lowest validation set
#' MSE) penalty parameter combination(s). If a grouping variable is used,
#' penalty parameters are tuned separately for each individual model.
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
#' @param nlambda The number of values for lambda0 (penalty parameter for L0
#'   penalty) - default is 100.
#' @param lambda.min.ratio Smallest value for lambda0, as a fraction of
#'   lambda0.max (data derived).
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
#' @return  An object of class \code{smimodel}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.}
#'   Each row of the column \code{fit} contains a list with three elements:
#'   \item{initial}{A list of information of the model initialisation. (For
#'   descriptions of the list elements see \code{\link{make_smimodelFit}}).}
#'   \item{best}{A list of information of the final optimised model. (For
#'   descriptions of the list elements see \code{\link{make_smimodelFit}}).}
#'   \item{best_lambdas}{Selected penalty parameter combination.} The number of
#'   rows of the \code{tibble} equals to the number of levels in the grouping
#'   variable.
#'
#' @references Palihawadana, N.K., Hyndman, R.J. & Wang, X. (2024). Sparse
#'   Multiple Index Models for High-Dimensional Nonparametric Forecasting.
#'   \url{https://www.monash.edu/business/ebs/research/publications/ebs/2024/wp16-2024.pdf}.
#'
#' @seealso \code{\link{model_smimodel}}
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
#' n = 1205
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
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Validation set
#' sim_val <- sim_data[1001:1200, ]
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' smi_greedy <- greedy_smimodel(data = sim_train,
#'                               val.data = sim_val,
#'                               yvar = "y1",
#'                               index.vars = index.vars,
#'                               initialise = "ppr",
#'                               lambda.min.ratio = 0.1)
#'
#' # Best (optimised) fitted model
#' smi_greedy$fit[[1]]$best
#'
#' # Selected penalty parameter combination
#' smi_greedy$fit[[1]]$best_lambdas
#' }
#'
#' @export
greedy_smimodel <- function(data, val.data, yvar, neighbour = 0,
                            family = gaussian(), index.vars,
                            initialise = c("ppr", "additive", "linear",
                                           "multiple", "userInput"),
                            num_ind = 5, num_models = 5, seed = 123,
                            index.ind = NULL, index.coefs = NULL,
                            s.vars = NULL, linear.vars = NULL,
                            nlambda = 100, lambda.min.ratio = 0.0001,
                            refit = TRUE, M = 10, max.iter = 50,
                            tol = 0.001, tolCoefs = 0.001,
                            TimeLimit = Inf, MIPGap = 1e-4, NonConvex = -1,
                            verbose = FALSE, parallel = FALSE, workers = NULL,
                            recursive = FALSE, recursive_colRange = NULL){

  # Message for gurobi installation
  message("Do you have Gurobi solver installed?
  Make sure you have an active installation of Gurobi solver (https://www.gurobi.com/)
  in your local machine before using this function.
  Refer the section 'Other Required Software' in the README for installation help.")

  # Check for `tsibble`
  stopifnot(tsibble::is_tsibble(data))
  stopifnot(tsibble::is_tsibble(val.data))

  initialise <- match.arg(initialise)
  # Data Preparation
  data1 <- data
  data2 <- val.data
  data_index <- index(data1)
  data_key <- key(data1)
  if (length(key(data1)) == 0) {
    data1 <- data1 |>
      mutate(dummy_key = rep(1, NROW(data1))) |>
      as_tsibble(index = data_index, key = dummy_key)
    data_key <- key(data1)
    data2 <- data2 |>
      mutate(dummy_key = rep(1, NROW(data2))) |>
      as_tsibble(index = data_index, key = dummy_key)
  }
  key11 <- key(data1)[[1]]
  key_unique <- unique(as.character(sort(dplyr::pull((data1[, {{ key11 }}])[, 1]))))
  key_num <- seq_along(key_unique)
  ref <- data.frame(key_unique, key_num)
  data1 <- data1 |>
    dplyr::mutate(
      num_key = as.numeric(factor(as.character({{ key11 }}), levels = key_unique))
    )
  data2 <- data2 |>
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
    df_cat_val <- data2 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    smimodels_list[[i]] <- greedy.fit(data = df_cat, val.data = df_cat_val,
                                      yvar = yvar,
                                      neighbour = neighbour,
                                      family = family, index.vars = index.vars,
                                      initialise = initialise,
                                      num_ind = num_ind, num_models = num_models,
                                      seed = seed, index.ind = index.ind,
                                      index.coefs = index.coefs, s.vars = s.vars,
                                      linear.vars = linear.vars,
                                      nlambda = nlambda,
                                      lambda.min.ratio = lambda.min.ratio,
                                      refit = refit,
                                      M = M, max.iter = max.iter,
                                      tol = tol, tolCoefs = tolCoefs,
                                      TimeLimit = TimeLimit, MIPGap = MIPGap,
                                      NonConvex = NonConvex, verbose = verbose,
                                      parallel = parallel, workers = workers,
                                      recursive = recursive,
                                      recursive_colRange = recursive_colRange)
  }
  data_list <- list(key_unique, smimodels_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("smimodel", "tbl_df", "tbl", "data.frame")
  return(models)
}
utils::globalVariables(c("dummy_key", "num_key"))


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
#' @param nlambda The number of values for lambda0 (penalty parameter for L0
#'   penalty) - default is 100.
#' @param lambda.min.ratio Smallest value for lambda0, as a fraction of
#'   lambda0.max (data derived).
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
                       nlambda = 100, lambda.min.ratio = 0.0001,
                       refit = TRUE, M = 10, max.iter = 50,
                       tol = 0.001, tolCoefs = 0.001,
                       TimeLimit = Inf, MIPGap = 1e-4, NonConvex = -1,
                       verbose = FALSE, parallel = FALSE, workers = NULL,
                       recursive = FALSE, recursive_colRange = NULL){

  ## Calculating lambda0.max based on the scale of the first term in the
  ## SMI modelling objective function
  # Fit an additive model (gam) as a benchmark
  bench <- model_gam(data = data,
                     yvar = yvar,
                     family = family,
                     neighbour = neighbour,
                     s.vars = c(index.vars, s.vars),
                     linear.vars = linear.vars)
  # Residuals
  bench_resid <- augment(bench)$.resid
  # lambda0_seq
  lambda0_max <- sum(bench_resid^2)
  lambda0_min <- lambda0_max * lambda.min.ratio
  lambda0_seq <- c(0, exp(seq(log(lambda0_min), log(lambda0_max), length.out = nlambda)))

  # lambda2_seq - a sequence of power of tens
  # lambda2.max is taken as the power of ten that matches the scale of lambda0.max
  max_power10 <- round(log10(abs(lambda0_max)))
  lambda2_seq <- c(0, 10^seq(-2, max_power10, by = 1))

  l0_len <- length(lambda0_seq)
  l2_len <- length(lambda2_seq)

  # Data frame for storing all combinations searched
  all_comb <- data.frame()
  # Vector for storing validation set MSEs of all combinations searched
  all_comb_mse <- numeric()
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
  start_l0 <- c(lambda0_seq[1], lambda0_seq[ceiling(l0_len / 2)], lambda0_seq[l0_len])
  start_l2 <- c(lambda2_seq[1], lambda2_seq[ceiling(l2_len / 2)], lambda2_seq[l2_len])
  lambda_comb <- expand.grid(start_l0, start_l2)

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
  # Updating searched combinations MSE
  all_comb_mse <- c(all_comb_mse, unlist(MSE_list))

  # Greedy search
  while(min_MSE < current_MSE){
    current_MSE <- min_MSE
    current_lambdas <- min_lambdas
    # Constructing new search space
    # lambda0
    l0_indx <- which(lambda0_seq == current_lambdas[1])
    if(l0_indx < l0_len){
      lambda0_seq_new <- c(lambda0_seq[l0_indx - 1], current_lambdas[1],
                           lambda0_seq[l0_indx + 1])
    }else{
      lambda0_seq_new <- c(lambda0_seq[l0_indx - 1], current_lambdas[1])
    }
    # lambda2
    l2_indx <- which(lambda2_seq == current_lambdas[2])
    if(l2_indx < l2_len){
      lambda2_seq_new <- c(lambda2_seq[l2_indx - 1], current_lambdas[2],
                           lambda2_seq[l2_indx + 1])
    }else{
      lambda2_seq_new <- c(lambda2_seq[l2_indx - 1], current_lambdas[2])
    }
    lambda_comb_new <- expand.grid(lambda0_seq_new, lambda2_seq_new)
    # Already searched combinations
    lambda_exist <- do.call(paste0, lambda_comb_new) %in% do.call(paste0, all_comb)
    lambda_comb <- lambda_comb_new[lambda_exist == FALSE, ]
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
      # Updating searched combinations MSE
      all_comb_mse <- c(all_comb_mse, unlist(MSE_list))
    }
  }
  # tibble with the information on all searched lambda combinations
  all_comb_info <- as_tibble(all_comb) |> 
    rename(lambda0 = Var1,
           lambda2 = Var2) |>
    mutate(val_MSE = all_comb_mse)
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
                 "best_lambdas" = current_lambdas,
                 "lambda0_seq" = lambda0_seq,
                 "lambda2_seq" = lambda2_seq,
                 "searched" = all_comb_info)
  return(output)
}



#' SMI model with a given penalty parameter combination
#'
#' Fits a nonparametric multiple index model to the data for a given combination
#' of the penalty parameters (lambda0, lambda2), and returns the validation set
#' mean squared error (MSE). (Used within \code{\link{greedy.fit}}; users are
#' not expected to use this function directly.)
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
#' @param lambda.comb A \code{numeric} vector (of length two) indicating the
#'   values for the two penalty parameters lambda0 and lambda2.
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
#' @return A \code{numeric}.
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

  # # Predictions on validation set
  pred <- predict(object = smimodel$best, newdata = val.data, recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  # Validation set MSE
  # Convert to a tibble
  index_val <- index(val.data)
  val.data <- val.data |>
    as_tibble() |>
    arrange({{index_val}})
  smimodel_mse <- MSE(residuals = (as.numeric(as.matrix(val.data[,{{yvar}}], ncol = 1)) - pred))
  return(smimodel_mse)
}
