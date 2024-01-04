#' Augment function for class `groupSmimodel`
#'
#' Generates residuals and fitted values of a fitted `groupSmimodel` object.
#'
#' @param x A `groupSmimodel` object.
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class`tsibble`.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param ... Other arguments not currently used.
#'
#'
#' @method augment groupSmimodel
#'
#' @export
augment.groupSmimodel <- function(x, data, neighbour = 0, ...) {
  if (!tsibble::is_tsibble(data)) stop("data is not a tsibble.")
  yvar <- x$fit[[1]]$var_y
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
  model.resid <- vector(mode = "list", length = NROW(x))
  model.fitted <- vector(mode = "list", length = NROW(x))
  time_variable <- vector(mode = "list", length = NROW(x))
  old_group <- vector(mode = "list", length = NROW(x))
  new_group <- vector(mode = "list", length = NROW(x))
  df <- vector(mode = "list", length = NROW(x))
  for (i in seq(NROW(x))) {
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) %>%
      tibble::as_tibble() %>%
      dplyr::arrange({{data_index}})
    model.resid[[i]] <- x$fit[[i]]$gam$residuals
    model.fitted[[i]] <- x$fit[[i]]$gam$fitted.values
    time_variable[[i]] <- df_cat[, data_index][[1]]
    old_group[[i]] <- df_cat[, key11][[1]]
    new_group[[i]] <- rep(x[[1]][i], nrow(df_cat))
    df[[i]] <- tibble::as_tibble(data.frame(Index = time_variable[[i]], 
                                            New = new_group[[i]],
                                            Old = old_group[[i]],
                                            .resid = model.resid[[i]],
                                            .fitted = model.fitted[[i]])) %>%
      dplyr::filter(New == Old)
  }
  mod_res <- dplyr::bind_rows(df)
  return(mod_res)
}
utils::globalVariables(c("New", "Old"))
#' @export
generics::augment
