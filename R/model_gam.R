#' Generalised Additive Models
#'
#' A wrapper for \code{mgcv::gam()} enabling multiple GAMs based on a grouping
#' variable.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}). If
#'   multiple models are fitted, the grouping variable should be the \code{key}
#'   of the \code{tsibble}. If a \code{key} is not specified, a dummy key with
#'   only one level will be created.
#' @param yvar Name of the response variable as a character string.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an \code{integer}. If \code{neighbour =
#'   x}, \code{x} number of keys before the key of interest and \code{x} number
#'   of keys after the key of interest are grouped together for model fitting.
#'   The default is \code{neighbour = 0} (i.e. no neighbours are considered for
#'   model fitting).
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted (i.e. non-linear predictors).
#' @param s.basedim Dimension of the bases used to represent the smooth terms
#'   corresponding to \code{s.vars}. (For more information refer
#'   \code{mgcv::s()}.)
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model (i.e. linear
#'   predictors).
#' @param ... Other arguments not currently used.
#' @return An object of class \code{gamFit}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.} Each row of the column \code{fit} is an
#'   object of class \code{gam}. For details refer \code{mgcv::gamObject}.
#'
#' @seealso \code{\link{model_smimodel}}, \code{\link{model_backward}},
#'   \code{\link{model_gaim}}, \code{\link{model_ppr}}, \code{\link{model_lm}}
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
#' # Predictors taken as non-linear variables
#' s.vars <- colnames(sim_data)[3:6]
#'
#' # Predictors taken as linear variables
#' linear.vars <- colnames(sim_data)[7:8]
#'
#' # Model fitting
#' gamModel <- model_gam(data = sim_data,
#'                       yvar = "y1",
#'                       s.vars = s.vars,
#'                       linear.vars = linear.vars)
#'
#' # Fitted model
#' gamModel$fit[[1]]
#'
#' @export
model_gam <- function(data, yvar, family = gaussian(), neighbour = 0, s.vars, 
                      s.basedim = NULL, linear.vars = NULL, ...){
  stopifnot(tsibble::is_tsibble(data))
  if (is.null(c(linear.vars, s.vars))) 
    stop("No predictor variables specified; s.vars = NULL, linear.vars = NULL.")
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
  if(!is.null(s.basedim)){
    pre.formula <- lapply(s.vars, function(var) paste0("s(", var, ', bs="cr",k=', s.basedim, ")")) |>
      paste(collapse = "+") 
  }else{
    pre.formula <- lapply(s.vars, function(var) paste0("s(", var, ', bs="cr")')) |>
      paste(collapse = "+") 
  }
  pre.formula <- paste(yvar, "~", pre.formula)
  if (!is.null(linear.vars)){
    linear.formula <- lapply(linear.vars, function(var) paste0(var)) |>
      paste(collapse = "+")
    pre.formula <- paste(pre.formula, "+", linear.formula)
  }

  # Model fitting
  gam_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    df_cat <- df_cat |>
      drop_na()
    gam_list[[i]] <- mgcv::gam(as.formula(pre.formula), family = family, 
                                method = "REML", data = df_cat, ... = ...)
    add <- df_cat |>
      #drop_na() |>
      select({{ data_index }}, {{ key11 }})
    gam_list[[i]]$model <- bind_cols(add, gam_list[[i]]$model)
    gam_list[[i]]$model <- as_tsibble(gam_list[[i]]$model,
                                     index = data_index,
                                     key = all_of(key11))
  }
  # Structuring the output
  data_list <- list(key_unique, gam_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("gamFit", "tbl_df", "tbl", "data.frame")
  return(models)
}