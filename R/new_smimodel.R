#' Constructor function for the class `smimodel`
#'
#' Constructs an object of class `smimodel` using the information passed to
#' arguments.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param yvar Name of the response variable as a character string.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param index.ind An integer vector that assigns group index for each
#'   predictor in `index.vars`. (The default is `NULL`. If `NULL`, the model
#'   will be initialised with an Additive Model.)
#' @param index.coefs A numeric vector of index coefficients (default: NULL).
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#'
#' @export
new_smimodel <- function(data, yvar, index.vars, index.ind = NULL, 
                         index.coefs = NULL, linear.vars = NULL){
  stopifnot(tibble::is_tibble(data))
  data <- data %>%
    drop_na()
  X_index <- as.matrix(data[ , index.vars])
  if(!is.null(index.ind) & !is.null(index.coefs)){
    # Index positions
    ind_pos <- split(seq_along(index.ind), index.ind)
    # Index coefficients
    alpha <- unlist(tapply(index.coefs, index.ind, normalise_alpha))
  }else{
    # Initialise the model with an Additive Model
    index.ind <- seq_along(index.vars)
    # Index positions
    ind_pos <- split(seq_along(index.ind), index.ind)
    # Index coefficients
    alpha <- index.coefs <- rep(1, length(index.vars))
  }
  # Calculating indices
  ind <- vector(length = length(ind_pos), mode = "list")
  for(i in 1:length(ind)){
    ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                             as.matrix(alpha[ind_pos[[i]]], ncol = 1))
  }
  dat_names <- names(ind) <- paste0("index", 1:length(ind))
  dat <- tibble::as_tibble(ind)
  # Fitting a `gam`
  # Constructing the formula
  pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ', bs="cr")')) %>%
    paste(collapse = "+") %>% 
    paste(yvar, "~", .)
  if (!is.null(linear.vars)){
    pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
      paste(collapse = "+") %>% 
      paste(pre.formula, "+", .)
  }
  # Model fitting
  dat_new <- dplyr::bind_cols(data, dat)
  fun1 <- mgcv::gam(as.formula(pre.formula), data = dat_new, method = "REML")
  smimodel <- make_smimodel(x = fun1, yvar = yvar, index.vars = index.vars, 
                            index.ind = index.ind, index.data = dat_new,
                            index.names = dat_names, alpha = alpha,
                            linear.vars = linear.vars)
  return(smimodel)
}