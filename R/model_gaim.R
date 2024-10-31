#' Groupwise Additive Index Models (GAIM)
#'
#' A wrapper for \code{cgaim::cgaim()} enabling multiple GAIM models based on a
#' grouping variable. Currently does not support Constrained GAIM (CGAIM)s.
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
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices should be estimated.
#' @param index.ind An \code{integer} vector that assigns group index for each
#'   predictor in \code{index.vars}.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted individually (rather than considering as
#'   part of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model.
#' @param ... Other arguments not currently used. (Note that the arguments in
#'   \code{cgaim::cgaim()} related to constrained GAIMs are currently not
#'   supported. Furthermore, the argument \code{subset} is also not supported
#'   due to a bug in \code{cgaim::cgaim()}.)
#' @return An object of class \code{gaimFit}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.} Each row of the column \code{fit} is an
#'   object of class \code{cgaim}. For details refer \code{cgaim::cgaim()}.
#'
#' @details Group-wise Additive Index Model (GAIM) can be written in the form
#' \deqn{y_{i} = \sum_{j = 1}^{p} g_{j}(\boldsymbol{\alpha}_{j}^{T}\boldsymbol{x}_{ij}) +
#' \varepsilon_{i}, \quad i = 1, \dots, n,} where \eqn{y_{i}} is the univariate
#'   response, \eqn{\boldsymbol{x}_{ij} \in \mathbb{R}^{l{j}}}, \eqn{j = 1,
#'   \dots, p} are pre-specified non-overlapping subsets of
#'   \eqn{\boldsymbol{x}_{i}}, and \eqn{\boldsymbol{\alpha}_j} are the
#'   corresponding index coefficients, \eqn{g_{j}} is an unknown (possibly
#'   nonlinear) component function, and \eqn{\varepsilon_{i}} is the random
#'   error, which is independent of \eqn{\boldsymbol{x}_{i}}.
#'
#' @seealso \code{\link{model_smimodel}}, \code{\link{model_backward}},
#'   \code{\link{model_ppr}}, \code{\link{model_gam}}, \code{\link{model_lm}}
#'
#' @examples
#' library(dplyr)
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
#' # Predictors taken as index variables
#' index.vars <- colnames(sim_data)[3:7]
#'
#' # Assign group indices for each predictor
#' index.ind = c(rep(1, 3), rep(2, 2))
#'
#' # Predictors taken as non-linear variables not entering indices
#' s.vars = "x_lag_005"
#'
#' # Model fitting
#' gaimModel <- model_gaim(data = sim_data,
#'                         yvar = "y1",
#'                         index.vars = index.vars,
#'                         index.ind = index.ind,
#'                         s.vars = s.vars)
#' # Fitted model
#' gaimModel$fit[[1]]
#'
#' @export
model_gaim <- function(data, yvar, neighbour = 0, index.vars, index.ind, 
                       s.vars = NULL, linear.vars = NULL, ...){
  stopifnot(tsibble::is_tsibble(data))
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
  # Constructing the formula
  ind_pos <- split(seq_along(index.ind), index.ind)
  var_list <- index.vars[ind_pos[[1]]]
  pre.formula <- lapply(var_list, function(var) paste0(var)) |>
    paste(collapse = ",") |> 
    paste0(")")
  pre.formula <- paste0(yvar, " ~ g(", pre.formula)
  if(length(ind_pos) > 1){
    for(j in 2:length(ind_pos)){
      var_list <- index.vars[ind_pos[[j]]]
      index.formula <- lapply(var_list, function(var) paste0(var)) |>
        paste(collapse = ",") |> 
        paste0(")")
      pre.formula <- paste0(pre.formula, " + g(", index.formula)
    }
  }
  if (!is.null(s.vars)){
    svars.formula <- lapply(s.vars, function(var) paste0("s(", var, ")")) |>
      paste(collapse = "+") 
    pre.formula <- paste(pre.formula, "+", svars.formula)
  }
  if (!is.null(linear.vars)){
    linear.formula <- lapply(linear.vars, function(var) paste0(var)) |>
      paste(collapse = "+") 
    pre.formula <- paste(pre.formula, "+", linear.formula)
  }
  # Model fitting
  gaim_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) 
    gaim_list[[i]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                   data = df_cat, ... = ...)
    modelFrame <- model.frame(formula = as.formula(pre.formula), 
                              data = df_cat)
    add <- df_cat |>
      drop_na() |>
      select({{ data_index }}, {{ key11 }})
    gaim_list[[i]]$model <- bind_cols(add, modelFrame)
    gaim_list[[i]]$model <- as_tsibble(gaim_list[[i]]$model,
                                       index = data_index,
                                       key = all_of(key11))
  }
  # Structuring the output
  data_list <- list(key_unique, gaim_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("gaimFit", "tbl_df", "tbl", "data.frame")
  return(models)
}