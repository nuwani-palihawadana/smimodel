#' SMI Model Estimation through a Greedy Search for Penalty Parameters
#'
#' Performs a greedy search over a given grid of penalty parameter combinations
#' (lambda0, lambda2), and fits SMI model(s) with best (lowest in-sample MSE)
#' penalty parameter combinations. If a grouping variable is used, penalty
#' parameters are tuned separately for each individual model.
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
greedy_smimodel <- function(data, yvar, neighbour = 0, 
                            family = gaussian(), index.vars, 
                            initialise = c("ppr", "additive", "linear", 
                                           "multiple", "userInput"),
                            num_ind = 5, num_models = 5, seed = 123, 
                            index.ind = NULL, index.coefs = NULL, 
                            s.vars = NULL, linear.vars = NULL, 
                            lambda0_seq, lambda2_seq, lambda_step,
                            lambda0_start_seq, lambda2_start_seq, 
                            M = 10, max.iter = 50, 
                            tol = 0.001, tolCoefs = 0.001,
                            TimeLimit = Inf, MIPGap = 1e-4, NonConvex = -1, 
                            verbose = FALSE, parallel = FALSE, workers = NULL){
  stopifnot(tsibble::is_tsibble(data))
  initialise <- match.arg(initialise)
  data1 <- data
  data_index <- index(data1)
  data_key <- key(data1)
  if (length(key(data1)) == 0) {
    data1 <- data1 %>%
      dplyr::mutate(dummy_key = rep(1, NROW(data1))) %>%
      tsibble::as_tsibble(index = data_index, key = dummy_key)
    data_key <- key(data1)
  }
  key11 <- key(data1)[[1]]
  key_unique <- unique(as.character(sort(dplyr::pull((data1[, {{ key11 }}])[, 1]))))
  key_num <- seq_along(key_unique)
  ref <- data.frame(key_unique, key_num)
  data1 <- data1 %>%
    dplyr::mutate(
      num_key = as.numeric(factor(as.character({{ key11 }}), levels = key_unique))
    )
  smimodels_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) 
    smimodels_list[[i]] <- greedy.fit(data = df_cat, yvar = yvar, 
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
                                      M = M, max.iter = max.iter, 
                                      tol = tol, tolCoefs = tolCoefs,
                                      TimeLimit = TimeLimit, MIPGap = MIPGap, 
                                      NonConvex = NonConvex, verbose = verbose, 
                                      parallel = parallel, workers = workers)
  }
  data_list <- list(key_unique, smimodels_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ vctrs::vec_as_names(..., repair = "universal", quiet = TRUE)
  )
  models <- models %>%
    dplyr::rename(key = ...1) %>%
    dplyr::rename(fit = ...2)
  class(models) <- c("smimodel", "tbl_df", "tbl", "data.frame")
  return(models)
}