#' Linear Regression models
#'
#' A wrapper for \code{\link{lm}} enabling multiple linear models based on a
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
#' @param linear.vars A character vector of names of the predictor variables.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{lmFit}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.} Each row of the column \code{fit} is
#'   an object of class \code{lm}. For details refer \code{stats::lm}.
#'
#' @seealso \code{\link{model_smimodel}}, \code{\link{model_backward}},
#'   \code{\link{model_gaim}}, \code{\link{model_ppr}}, \code{\link{model_gam}}
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
#' # Predictor variables
#' linear.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' lmModel <- model_lm(data = sim_data,
#'                     yvar = "y1",
#'                     linear.vars = linear.vars)
#' # Fitted model
#' lmModel$fit[[1]]
#'
#' @export
model_lm <- function(data, yvar, neighbour = 0, linear.vars, ...){
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
  pre.formula <- lapply(linear.vars, function(var) paste0(var)) |>
    paste(collapse = "+") 
  pre.formula <- paste(yvar, "~", pre.formula)
  # Model fitting
  lm_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    df_cat <- df_cat |>
      drop_na()
    lm_list[[i]] <- stats::lm(formula = as.formula(pre.formula),
                              data = df_cat, ... = ...)
    add <- df_cat |>
      drop_na() |>
      select({{ data_index }}, {{ key11 }})
    lm_list[[i]]$model <- bind_cols(add, lm_list[[i]]$model)
    lm_list[[i]]$model <- as_tsibble(lm_list[[i]]$model,
                                     index = data_index,
                                     key = all_of(key11))
  }
  # Structuring the output
  data_list <- list(key_unique, lm_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("lmFit", "tbl_df", "tbl", "data.frame")
  return(models)
}