#' Constructor function for the class `smimodelFit`
#'
#' Constructs an object of class `smimodelFit` using the information passed to
#' arguments.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is "additive", where the initial model
#'   will be a nonparametric additive model. The other options are "linear" -
#'   linear regression model (i.e. a special case single-index model, where the
#'   initial values of the index coefficients are obtained through a linear
#'   regression), and "userInput" - user specifies the initial model structure
#'   (i.e. the number of indices and the placement of index variables among
#'   indices) and the initial index coefficients through `index.ind` and
#'   `index.coefs` arguments respectively.
#' @param index.ind If `initialise = "userInput"`: an integer vector that
#'   assigns group index for each predictor in `index.vars`.
#' @param index.coefs If `initialise = "userInput"`: a numeric vector of index
#'   coefficients.
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted individually (rather than considering as a
#'   part of an index considered in `index.vars`).
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.

new_smimodelFit <- function(data, yvar, neighbour = 0, 
                            family = gaussian(), index.vars, 
                            initialise = c("additive", "linear", "userInput"), 
                            index.ind = NULL, index.coefs = NULL, 
                            s.vars = NULL, linear.vars = NULL){
  stopifnot(tsibble::is_tsibble(data))
  data_index <- index(data)
  data_key <- key(data)[[1]]
  initialise <- match.arg(initialise)
  data <- data %>%
    drop_na()
  Y_data <- as.matrix(data[ , yvar])
  X_index <- as.matrix(data[ , index.vars])
  if(initialise == "linear"){
    # Initialise the model with a linear model
    index.ind <- rep(1, length(index.vars))
    ind_pos <- split(seq_along(index.ind), index.ind)
    pre.formula <- lapply(index.vars, function(var) paste0(var)) %>%
      paste(collapse = "+") %>% 
      paste(yvar, "~", .)
    if (!is.null(s.vars)){
      pre.formula <- lapply(s.vars, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    if (!is.null(linear.vars)){
      pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    pre.formula <- paste(pre.formula, "-", 1)
    fun1 <- mgcv::gam(as.formula(pre.formula), data = data, family = family,
                      method = "REML")
    add <- data %>%
      drop_na() %>%
      select({{ data_index }}, {{ data_key }})
    fun1$model <- bind_cols(add, fun1$model)
    fun1$model <- as_tsibble(fun1$model,
                             index = data_index,
                             key = all_of(data_key))
    # Index coefficients
    alpha <- index.coefs <- fun1$coefficients[1:length(index.vars)]
    alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    # Calculating indices
    ind <- vector(length = length(ind_pos), mode = "list")
    for(i in 1:length(ind)){
      if(length(ind_pos[[i]]) == 1){
        temp <- i
      }else{
        temp <- unlist(lapply(1:length(ind_pos[[i]]), function(x) paste0(i, x))) 
      }
      ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                               as.matrix(alpha[match(temp, names(alpha))], ncol = 1))
    }
    dat_names <- names(ind) <- paste0("index", 1:length(ind))
    dat <- tibble::as_tibble(ind)
    dat_new <- dplyr::bind_cols(data, dat)
  }else{
    if(initialise == "additive"){
      # Initialise the model with a nonparametric additive model
      index.ind <- seq_along(index.vars)
      # Index positions
      ind_pos <- split(seq_along(index.ind), index.ind)
      # Index coefficients
      alpha <- index.coefs <- rep(1, length(index.vars))
      alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    }else if(initialise == "userInput"){
      if (is.null(index.ind) | is.null(index.coefs)) stop("index.ind and/or index.coefs are/is not provided.")
      # Number of (index) predictors
      num_pred <- length(index.vars)
      # Index positions
      ind_pos <- split(seq_along(index.ind), index.ind)
      # Number of indices
      num_ind <- length(ind_pos)
      # Index coefficients
      alpha <- unlist(tapply(index.coefs, index.ind, normalise_alpha))
      # Constructing a new index coefficient vector to have all predictors in each index
      newIndex <- allpred_index(num_pred = num_pred,
                                num_ind = num_ind,
                                ind_pos = ind_pos,
                                alpha = alpha)
      alpha <- newIndex$alpha_init_new
      index.ind <- newIndex$index
      ind_pos <- newIndex$index_positions
      # Adjusting X (matrix of predictors) to fit number of indices
      X_index <- do.call(cbind, replicate(num_ind, X_index, simplify = FALSE))
      # Checking for all zero indices
      ind_rm_id <- numeric()
      ind_rm_pos <- numeric()
      for(i in 1:num_ind){
        if(all(alpha[ind_pos[[i]]] == 0)){
          ind_rm_id <- c(ind_rm_id, i)
          ind_rm_pos <- c(ind_rm_pos, ind_pos[[i]])
          warning(paste0('Initial model', ': All coefficients of index', i,
                         ' are zero. Removing index', i, 
                         '. However, the variables in the removed index are considered in subsequent model searches.')) 
        }
      }
      if(length(ind_rm_id) != 0){
        X_index <- as.matrix(X_index[ , -ind_rm_pos])
        alpha <- alpha[-ind_rm_pos]
        num_ind  <- num_ind - length(ind_rm_id)
        index.ind <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          index.ind[[i]] <- rep(i, num_pred)
        }
        index.ind <- unlist(index.ind)
        ind_pos <- split(1:length(index.ind), index.ind)
      }
      names(alpha) <- NULL
      alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    }
    # Calculating indices
    ind <- vector(length = length(ind_pos), mode = "list")
    for(i in 1:length(ind)){
      if(length(ind_pos[[i]]) == 1){
        temp <- i
      }else{
        temp <- unlist(lapply(1:length(ind_pos[[i]]), function(x) paste0(i, x))) 
      }
      ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                               as.matrix(alpha[match(temp, names(alpha))], ncol = 1))
    }
    dat_names <- names(ind) <- paste0("index", 1:length(ind))
    dat <- tibble::as_tibble(ind)
    # Constructing the formula
    pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ', bs="cr")')) %>%
      paste(collapse = "+") %>% 
      paste(yvar, "~", .)
    if (!is.null(s.vars)){
      pre.formula <- lapply(s.vars, function(var) paste0("s(", var, ', bs="cr")')) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    if (!is.null(linear.vars)){
      pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    # Model fitting
    dat_new <- dplyr::bind_cols(data, dat)
    dat_new <- dat_new %>%
      tsibble::as_tsibble(index = {{data_index}}, key = {{data_key}})
    fun1 <- mgcv::gam(as.formula(pre.formula), data = dat_new, family = family,
                      method = "REML")
    add <- dat_new %>%
      drop_na() %>%
      select({{ data_index }}, {{ data_key }})
    fun1$model <- bind_cols(add, fun1$model)
    fun1$model <- as_tsibble(fun1$model,
                             index = data_index,
                             key = all_of(data_key))
  }
  smimodel <- make_smimodelFit(x = fun1, yvar = yvar, 
                               neighbour = neighbour,
                               index.vars = index.vars, 
                               index.ind = index.ind, index.data = dat_new,
                               index.names = dat_names, alpha = alpha,
                               s.vars = s.vars, linear.vars = linear.vars)
  return(smimodel)
}