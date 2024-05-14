#' Groupwise Additive Index Models (GAIM)
#'
#' A wrapper for `cgaim::cgaim()` enabling multiple GAIM models based on a
#' grouping variable. Currently does not support Constrained GAIM (CGAIM)s.
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
#' @param index.ind An integer vector that assigns group index for each
#'   predictor in `index.vars`.
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted individually (rather than considering as a
#'   part of an index considered in `index.vars`).
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#' @param ... Other arguments not currently used. (Note that the arguments in
#'   `cgaim::cgaim()` related to constrained GAIMs are currently not supported.
#'   Furthermore, the argument `subset` is also not supported due to a bug in
#'   `cgaim::cgaim()`.)
#'
#' @importFrom cgaim cgaim
#' 
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) %>%
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) %>%
#'   unpack(x_lag, names_sep = "_") %>%
#'   mutate(
#'     # Response variable
#'     y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) %>%
#'   drop_na() %>%
#'   select(inddd, y1, starts_with("x_lag")) %>%
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#' # Index variables
#' index.vars <- colnames(sim_data)[3:7]
#' # Assign group indices for each predictor
#' index.ind = c(rep(1, 3), rep(2, 2))
#' # Smooth variables not entering indices
#' s.vars = "x_lag_005"
#' # Model fitting
#' gaimModel <- model_gaim(data = sim_data,
#'                         yvar = "y1",
#'                         index.vars = index.vars,
#'                         index.ind = index.ind,
#'                         s.vars = s.vars)
#' gaimModel$fit[[1]]
#' @export
model_gaim <- function(data, yvar, neighbour = 0, index.vars, index.ind, 
                       s.vars = NULL, linear.vars = NULL, ...){
  stopifnot(tsibble::is_tsibble(data))
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
  # Constructing the formula
  ind_pos <- split(seq_along(index.ind), index.ind)
  var_list <- index.vars[ind_pos[[1]]]
  pre.formula <- lapply(var_list, function(var) paste0(var)) %>%
    paste(collapse = ",") %>% 
    paste0(")") %>%
    paste0(yvar, " ~ g(", .)
  if(length(ind_pos) > 1){
    for(j in 2:length(ind_pos)){
      var_list <- index.vars[ind_pos[[j]]]
      pre.formula <- lapply(var_list, function(var) paste0(var)) %>%
        paste(collapse = ",") %>% 
        paste0(")") %>%
        paste0(pre.formula, " + g(", .)
    }
  }
  if (!is.null(s.vars)){
    pre.formula <- lapply(s.vars, function(var) paste0("s(", var, ")")) %>%
      paste(collapse = "+") %>% 
      paste(pre.formula, "+", .)
  }
  if (!is.null(linear.vars)){
    pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
      paste(collapse = "+") %>% 
      paste(pre.formula, "+", .)
  }
  # Model fitting
  gaim_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) 
    gaim_list[[i]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                   data = df_cat, ... = ...)
    modelFrame <- model.frame(formula = as.formula(pre.formula), 
                              data = df_cat)
    add <- df_cat %>%
      drop_na() %>%
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
    .name_repair = ~ vctrs::vec_as_names(..., repair = "universal", quiet = TRUE)
  )
  models <- models %>%
    dplyr::rename(key = ...1) %>%
    dplyr::rename(fit = ...2)
  class(models) <- c("gaimFit", "tbl_df", "tbl", "data.frame")
  return(models)
}
#' @export
cgaim::g
#' @export
cgaim::s