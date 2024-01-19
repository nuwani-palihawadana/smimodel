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
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#'
#' @export
new_smimodel <- function(data, yvar, index.vars, 
                         initialise = c("additive", "linear", "userInput"), 
                         index.ind = NULL, index.coefs = NULL, 
                         linear.vars = NULL){
  stopifnot(tibble::is_tibble(data))
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
    if (!is.null(linear.vars)){
      pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    pre.formula <- paste(pre.formula, "-", 1)
    fun1 <- mgcv::gam(as.formula(pre.formula), data = data, method = "REML")
    # Index coefficients
    alpha <- index.coefs <- fun1$coefficients[1:length(index.vars)]
    alpha <- unlist(tapply(alpha, index.ind, normalise_alpha))
    # Calculating indices
    ind <- vector(length = length(ind_pos), mode = "list")
    for(i in 1:length(ind)){
      ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                               as.matrix(alpha[startsWith(names(alpha), paste0(i))], ncol = 1))
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
      # Index positions
      ind_pos <- split(seq_along(index.ind), index.ind)
      # Index coefficients
      alpha <- unlist(tapply(index.coefs, index.ind, normalise_alpha))
      # Checking for all zero indices
      ind_rm_id <- numeric()
      ind_rm_pos <- numeric()
      num_ind <- length(ind_pos)
      for(j in 1:num_ind){
        if(all(alpha[ind_pos[[j]]] == 0)){
          ind_rm_id <- c(ind_rm_id, j)
          ind_rm_pos <- c(ind_rm_pos, ind_pos[[j]])
          warning(paste0('Initial model: ', ': All coefficients of index', j,
                         ' are zero. Removing index', j, 
                         ' . The variables in the removed index will be considered in subsequent model searches.')) 
        }
      }
      if(length(ind_rm_id) != 0){
        X_index <- as.matrix(X_index[ , -ind_rm_pos])
        alpha <- alpha[-ind_rm_pos]
        num_ind  <- num_ind - length(ind_rm_id)
        index.ind <- index.ind[-ind_rm_pos]
        ind_pos <- split(1:length(index.ind), index.ind)
      }
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
    if (!is.null(linear.vars)){
      pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    # Model fitting
    dat_new <- dplyr::bind_cols(data, dat)
    fun1 <- mgcv::gam(as.formula(pre.formula), data = dat_new, method = "REML")
  }
  smimodel <- make_smimodel(x = fun1, yvar = yvar, index.vars = index.vars, 
                            index.ind = index.ind, index.data = dat_new,
                            index.names = dat_names, alpha = alpha,
                            linear.vars = linear.vars)
  return(smimodel)
}