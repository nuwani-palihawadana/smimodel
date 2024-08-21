#' Generalised Additive Models
#'
#' A wrapper for `mgcv::gam()` enabling multiple GAMs based on a grouping
#' variable.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class `tsibble`.(Make sure there are no additional
#'   date/time/date-time/yearmonth/POSIXct/POSIXt variables except for the
#'   `index` of the `tsibble`). If multiple models are fitted, the grouping
#'   variable should be the key of the `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted (i.e. non-linear predictors). 
#' @param s.basedim Dimension of the bases used to represent the smooth terms
#'   corresponding to `s.vars`. (For more information refer `mgcv::s()`.)
#'   (default: NULL)
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model. (default:
#'   NULL)
#' @param ... Other arguments not currently used.
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
  pre.formula <- paste0(yvar, " ~ ")
  if (!is.null(s.basedim)) {
    pre.formula <- paste0(
      pre.formula, "+s(", paste0(s.vars), ',bs="cr",k=',
      paste0(s.basedim), ")"
    )
  } else {
    pre.formula <- paste0(pre.formula, "+s(", paste0(s.vars), ',bs="cr")')
  }
  if (!is.null(linear.vars)){
    pre.formula <- paste0(pre.formula, "+", paste0(linear.vars))
  }
  # Model fitting
  gam_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) 
     gam_list[[i]] <- mgcv::gam(as.formula(pre.formula), family = family, 
                                method = "REML", data = df_cat, ... = ...)
    add <- df_cat |>
      drop_na() |>
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