#' Sparse Multiple Index (SMI) models based on a grouping variable
#'
#' Fits nonparametric multiple index models, with simultaneous variable
#' selection for each group based on a grouping variable of interest.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class `tsibble`.(Make sure there are no additional
#'   date/time/date-time/yearmonth/POSIXct/POSIXt variables except for the
#'   `index` of the `tsibble`). If multiple models are fitted, the grouping
#'   variable should be the key of the `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is "additive", where the initial model
#'   will be a nonparametric additive model. The other options are "linear" -
#'   linear regression model (i.e. a special case single-index model, where the
#'   initial values of the index coefficients are obtained through a linear
#'   regression), "multiple" - multiple models are fitted starting with
#'   different initial models (single-index (linear), 2-index, 3-index and
#'   5-index models, where the predictor assignment to indices and initial index
#'   coefficients are generated randomly), and the final optimal model with
#'   lowest loss is returned, and "userInput" - user specifies the initial model
#'   structure (i.e. the number of indices and the placement of index variables
#'   among indices) and the initial index coefficients through `index.ind` and
#'   `index.coefs` arguments respectively.
#' @param index.ind If `initialise = "userInput"`: an integer vector that
#'   assigns group index for each predictor in `index.vars`.
#' @param index.coefs If `initialise = "userInput"`: a numeric vector of index
#'   coefficients.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
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
#'
#' @importFrom dplyr arrange filter mutate rename
#' @importFrom tsibble as_tsibble is_tsibble index key
#' @importFrom vctrs vec_as_names
#'
#' @export
groupSmimodel <- function(data, yvar, index.vars, 
                          initialise = c("additive", "linear", "multiple", "userInput"),
                          index.ind = NULL, index.coefs = NULL, 
                          neighbour = 0, linear.vars = NULL, 
                          lambda0 = 1, lambda2 = 1, 
                          M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                          TimeLimit = Inf, verbose = FALSE){
  stopifnot(tsibble::is_tsibble(data))
  initialise <- match.arg(initialise)
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
  smimodels_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) %>%
      tibble::as_tibble() %>%
      dplyr::arrange({{data_index}})
    smimodels_list[[i]] <- smimodel(data = df_cat, yvar = yvar, 
                                    index.vars = index.vars, 
                                    initialise = initialise, 
                                    index.ind = index.ind, 
                                    index.coefs = index.coefs,
                                    linear.vars = linear.vars,
                                    lambda0 = lambda0, lambda2 = lambda2, 
                                    M = M, max.iter = max.iter, 
                                    tol = tol, tolCoefs = tolCoefs,
                                    TimeLimit = TimeLimit, verbose = verbose)
  }
  data_list <- list(key_unique, smimodels_list)
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