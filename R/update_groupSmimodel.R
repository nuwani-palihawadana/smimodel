#' Updating a `groupSmimodel`
#'
#' Optimises and updates a given `groupSmimodel`.
#'
#' @param object A `groupSmimodel` object.
#' @param data Training data set on which models will be trained. Must be a
#'   data set of class `tsibble`.(Make sure there are no additional
#'   date/time/date-time/yearmonth/POSIXct/POSIXt variables except for the
#'   `index` of the `tsibble`). If multiple models are fitted, the grouping
#'   variable should be the key of the `tsibble`.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param tolCoefs Tolerance for coefficients. 
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param verbose The option to print detailed solver output.
#' @param ... Other arguments not currently used.
#'
#' @export
update_groupSmimodel <- function(object, data, neighbour = 0, lambda0 = 1, lambda2 = 1, 
                                 M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                                 TimeLimit = Inf, verbose = FALSE, ...){
  if (!tsibble::is_tsibble(data)) stop("data is not a tsibble.")
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
  smimodels_list <- vector(mode = "list", length = NROW(object))
  loss_list <- vector(mode = "list", length = NROW(object))
  yvar <- object$fit[[1]]$var_y
  for (i in 1:NROW(object)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) %>%
      tibble::as_tibble() %>%
      dplyr::arrange({{data_index}})
    update_temp <- update_smimodel(object = object$fit[[i]], 
                                   data = df_cat, 
                                   lambda0 = lambda0, lambda2 = lambda2, 
                                   M = M, max.iter = max.iter, 
                                   tol = tol, tolCoefs = tolCoefs,
                                   TimeLimit = TimeLimit, verbose = verbose)
    smimodels_list[[i]] <- update_temp
  }
  data_list <- list(object$key, smimodels_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ vctrs::vec_as_names(..., repair = "universal", quiet = TRUE)
  )
  models <- models %>%
    dplyr::rename(key = ...1) %>%
    dplyr::rename(fit = ...2)
  class(models) <- c("groupSmimodel", "tbl_df", "tbl", "data.frame")
  return(models)
}