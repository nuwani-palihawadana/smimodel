#' SMI model estimation through a greedy search for penalty parameters
#'
#' Performs a greedy search over a given grid of penalty parameter combinations
#' (lambda0, lambda2), and fits SMI model(s) with best (lowest validation set
#' MSE) penalty parameter combination(s). If a grouping variable is used, penalty
#' parameters are tuned separately for each individual model.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tsibble`.
#' @param val.data Validation data set. (The data set on which the penalty
#'   parameter selection will be performed.) Must be a data set of class
#'   `tsibble`. (Once the penalty parameter selection is completed, the best
#'   model will be re-fitted for the combined data set `data` + `val.data`.)
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
#' @param refit Whether to refit the model combining training and validation
#'   sets after parameter tuning. If `FALSE`, the final model will be estimated
#'   only on the training set.
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
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `val.data` to be filled with forecasts.
#'
#' @examples
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
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
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Validation set
#' sim_val <- sim_data[1001:1200, ]
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#' # Penalty parameter values to search
#' # L0 penalty
#' lambda0 = seq(1, 12, by = 1)
#' # L2 penalty
#' lambda2 = seq(0, 12, by = 1)
#' # Full grid
#' grid1 <- expand.grid(lambda0, lambda2)
#' # Starting point options
#' starting <- grid1[c(1, 6, 12, 73, 78, 84, 145, 150, 156), ]
#' # L0 penalty
#' lambda0_start = as.numeric(unique(unlist(starting[1])))
#' # L2 penalty
#' lambda2_start = as.numeric(unique(unlist(starting[2])))
#' # Model fitting
#' smi_greedy <- greedy_smimodel(data = sim_train,
#'                               val.data = sim_val,
#'                               yvar = "y1",
#'                               index.vars = index.vars,
#'                               initialise = "additive",
#'                               lambda0_seq = lambda0,
#'                               lambda2_seq = lambda2,
#'                               lambda_step = 1,
#'                               lambda0_start_seq = lambda0_start,
#'                               lambda2_start_seq = lambda2_start)
#' smi_greedy$fit[[1]]
#' @export
greedy_smimodel <- function(data, val.data, yvar, neighbour = 0, 
                            family = gaussian(), index.vars, 
                            initialise = c("ppr", "additive", "linear", 
                                           "multiple", "userInput"),
                            num_ind = 5, num_models = 5, seed = 123, 
                            index.ind = NULL, index.coefs = NULL, 
                            s.vars = NULL, linear.vars = NULL, 
                            lambda0_seq, lambda2_seq, lambda_step,
                            lambda0_start_seq, lambda2_start_seq, 
                            refit = TRUE, M = 10, max.iter = 50, 
                            tol = 0.001, tolCoefs = 0.001,
                            TimeLimit = Inf, MIPGap = 1e-4, NonConvex = -1, 
                            verbose = FALSE, parallel = FALSE, workers = NULL,
                            recursive = FALSE, recursive_colRange = NULL){
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
                                      lambda0_seq = lambda0_seq, 
                                      lambda2_seq = lambda2_seq, 
                                      lambda_step = lambda_step,
                                      lambda0_start_seq = lambda0_start_seq, 
                                      lambda2_start_seq = lambda2_start_seq, 
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