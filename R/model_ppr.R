#' Projection Pursuit Regression (PPR) models
#'
#' A wrapper for `stats::ppr()` enabling multiple PPR models based on a grouping
#' variable.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class `tsibble`.(Make sure there are no additional
#'   date/time/date-time/yearmonth/POSIXct/POSIXt variables except for the
#'   `index` of the `tsibble`). If multiple models are fitted, the grouping
#'   variable should be the key of the `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param num_ind An integer that specifies the number of indices to be used in
#'   the model(s). (Corresponds to `nterms` in `stats::ppr()`.)
#' @param ... Other arguments not currently used. (For more information on other
#'   arguments that can be passed, refer `stats::ppr()`.)
#' 
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
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
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#' # Model fitting
#' pprModel <- model_ppr(data = sim_data,
#'                       yvar = "y1",
#'                       index.vars = index.vars)
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