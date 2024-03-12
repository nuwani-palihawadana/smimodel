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
#' @importFrom stats ppr model.frame
#'
#' @export
model_ppr <- function(data, yvar, neighbour = 0, index.vars, num_ind = 5, ...){
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
  pre.formula <- lapply(index.vars, function(var) paste0(var)) %>%
    paste(collapse = "+") %>% 
    paste(yvar, "~", .)
  # Model fitting
  ppr_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) 
    ppr_list[[i]] <- stats::ppr(formula = as.formula(pre.formula),
                                data = df_cat, nterms = num_ind, ... = ...)
    modelFrame <- model.frame(formula = as.formula(pre.formula), data = df_cat)
    add <- df_cat %>%
      drop_na() %>%
      select({{ data_index }}, {{ key11 }})
    ppr_list[[i]]$model <- bind_cols(add, modelFrame)
    ppr_list[[i]]$model <- as_tsibble(ppr_list[[i]]$model,
                                      index = data_index,
                                      key = all_of(key11))
  }
  # Structuring the output
  data_list <- list(key_unique, ppr_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ vctrs::vec_as_names(..., repair = "universal", quiet = TRUE)
  )
  models <- models %>%
    dplyr::rename(key = ...1) %>%
    dplyr::rename(fit = ...2)
  class(models) <- c("pprFit", "tbl_df", "tbl", "data.frame")
  return(models)
}