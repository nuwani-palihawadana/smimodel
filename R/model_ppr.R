#' Projection Pursuit Regression (PPR) models
#'
#' A wrapper for \code{stats::ppr()} enabling multiple PPR models based on a
#' grouping variable.
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
#' @param num_ind An \code{integer} that specifies the number of indices to be
#'   used in the model(s). (Corresponds to \code{nterms} in
#'   \code{stats::ppr()}.)
#' @param ... Other arguments not currently used. (For more information on other
#'   arguments that can be passed, refer \code{stats::ppr()}.)
#' @return An object of class \code{pprFit}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.} Each row of the column \code{fit} is an
#'   object of class \code{c("ppr.form", "ppr")}. For details refer
#'   \code{stats::ppr()}.
#'
#' @details A Projection Pursuit Regression (PPR) model (Friedman & Stuetzle
#'   (1981)) is given by
#' \deqn{y_{i} = \sum_{j=1}^{p} {g_{j}(\boldsymbol{\alpha}_{j}^{T}\boldsymbol{x}_{i})} +
#' \varepsilon_{i}, \quad i = 1, \dots, n,} where \eqn{y_{i}} is the response,
#'   \eqn{\boldsymbol{x}_{i}} is the \eqn{q}-dimensional predictor vector,
#'   \eqn{\boldsymbol{\alpha}_{j} = ( \alpha_{j1}, \dots, \alpha_{jp} )^{T}},
#'   \eqn{j = 1, \dots, p} are \eqn{q}-dimensional projection vectors (or
#'   vectors of "index coefficients"), \eqn{g_{j}}'s are unknown nonlinear
#'   functions, and \eqn{\varepsilon_{i}} is the random error.
#'
#' @references Friedman, J. H. & Stuetzle, W. (1981). Projection pursuit
#'   regression. *Journal of the American Statistical Association*, 76, 817–823.
#'   \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.1981.10477729}{doi:10.2307/2287576}.
#'
#' @seealso \code{\link{model_smimodel}}, \code{\link{model_backward}},
#'   \code{\link{model_gaim}}, \code{\link{model_gam}}, \code{\link{model_lm}}
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
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' pprModel <- model_ppr(data = sim_data,
#'                       yvar = "y1",
#'                       index.vars = index.vars)
#'
#' # Fitted model
#' pprModel$fit[[1]]
#'
#' @export
model_ppr <- function(data, yvar, neighbour = 0, index.vars, num_ind = 5, ...){
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
  pre.formula <- lapply(index.vars, function(var) paste0(var)) |>
    paste(collapse = "+") 
  pre.formula <- paste(yvar, "~", pre.formula)
  # Model fitting
  ppr_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    df_cat <- df_cat |>
      drop_na()
    ppr_list[[i]] <- stats::ppr(formula = as.formula(pre.formula),
                                data = df_cat, nterms = num_ind, ... = ...)
    dot_args <- list(...)
    if ("subset" %in% names(dot_args)){
      modelFrame <- model.frame(formula = as.formula(pre.formula), 
                                data = df_cat[dot_args$subset, ])
      add <- df_cat[dot_args$subset, ] |>
        drop_na() |>
        select({{ data_index }}, {{ key11 }})
    }else{
      modelFrame <- model.frame(formula = as.formula(pre.formula), 
                                data = df_cat)
      add <- df_cat |>
        drop_na() |>
        select({{ data_index }}, {{ key11 }})
    }
    ppr_list[[i]]$model <- bind_cols(add, modelFrame)
    ppr_list[[i]]$model <- as_tsibble(ppr_list[[i]]$model,
                                      index = data_index,
                                      key = all_of(key11))
  }
  # Structuring the output
  data_list <- list(key_unique, ppr_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("pprFit", "tbl_df", "tbl", "data.frame")
  return(models)
}