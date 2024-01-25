#' Converting a `smimodel` object to a `gam` object
#'
#' Converts a given object of class `smimodel` to an object of class `gam`.
#'
#' @param x A `smimodel` object.
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#'
#' @importFrom dplyr bind_cols
#' @importFrom mgcv gam
#' @importFrom magrittr %>%
#' @importFrom stats as.formula
#' @importFrom tibble as_tibble is_tibble
#' @importFrom tidyr drop_na
#'
#' @export
make_gam <- function(x, data){
  if (!tibble::is_tibble(data)) stop("data is not a tibble.")
  data <- data %>%
    drop_na()
  X_index <- as.matrix(data[ , x$vars_index])
  list_index <- x$alpha[ , 2:NCOL(x$alpha)]
  alpha <- vector(mode = "list", length = length(list_index))
  for(i in 1:length(list_index)){
    alpha[[i]] <- list_index[ , i]
  }
  alpha <- unlist(alpha)
  if(all(alpha == 0)){
    # Constructing the formula and model fitting
    if (!is.null(x$vars_linear)){
      pre.formula <- lapply(x$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(x$var_y, "~", .)
    }else{
      pre.formula <- paste(x$var_y, "~", 1)
    }
    fun1_final <- mgcv::gam(as.formula(pre.formula), data = data, method = "REML")
  }else{
    # Calculating indices
    ind <- vector(mode = "list", length = length(list_index))
    for(i in 1:length(ind)){
      ind[[i]] <- as.numeric(X_index %*% as.matrix(list_index[ , i], ncol = 1))
    }
    names(ind) <- names(list_index)
    dat <- tibble::as_tibble(ind)
    # Fitting a `gam`
    # Constructing the formula
    yvar <- x$var_y
    pre.formula <- lapply(names(list_index), function(var) paste0("s(", var, ',bs="cr")')) %>%
      paste(collapse = "+") %>% 
      paste(yvar, "~", .)
    if (!is.null(x$vars_linear)){
      pre.formula <- lapply(x$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    # Model fitting
    dat_new <- dplyr::bind_cols(data, dat)
    fun1 <- mgcv::gam(as.formula(pre.formula), data = dat_new, method = "REML")
  }
  return(fun1)
}
utils::globalVariables(c("."))