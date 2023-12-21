#' Updating a `smimodel` 
#' 
#' Optimises and updates a given `smimodel`. 
#' 
#' @param object A `smimodel` object.
#' @param data Training data set on which models will be trained. Should be a `tibble`. 
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param verbose The option to print detailed solver output.
#' @param ... Other arguments not currently used.
#'
#' @export
update_smimodel <- function(object, data, lambda0 = 1, lambda2 = 1, 
                            M = 10, max.iter = 50, tol = 0.001, 
                            TimeLimit = Inf, verbose = FALSE, ...){
  if (!tibble::is_tibble(data)) stop("data is not a tibble.")
  data <- data %>%
    drop_na()
  gam1 <- make_gam(x = object, data = data)
  # Preparing inputs to `inside_update()`
  list_index <- object[1:(length(object)-3)]
  num_ind <- length(list_index)
  alpha <- vector(mode = "list", length = num_ind)
  dgz <- vector(mode = "list", length = num_ind)
  dgz_names <- character(length = num_ind)
  for(i in 1:num_ind){
    alpha[[i]] <- list_index[[i]]$coefficients
    dgz[[i]] <- list_index[[i]]$derivatives
    dgz_names[i] <- paste0("d", i)
  }
  alpha <- unlist(alpha)
  names(dgz) <- dgz_names
  dgz <- as.matrix(tibble::as_tibble(dgz))
  # Optimising the model
  best_alpha1 <- inner_update(x = gam1, data = data, yvar = object$var_y,
                              index.vars = object$vars_index, 
                              linear.vars = object$vars_linear, 
                              num_ind = num_ind, dgz = dgz, 
                              alpha_old = alpha, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, max.iter = max.iter, 
                              tol = tol, TimeLimit = TimeLimit, verbose = verbose)
  # Checking models with higher number of indices
  alpha_current <- best_alpha1$best_alpha
  MSE_current <- best_alpha1$min_loss
  index_current <- best_alpha1$index.ind
  ind_pos_current <- best_alpha1$ind_pos
  X_new_current <- best_alpha1$X_new
  j <- num_ind + 1
  num_pred <- length(object$vars_index)
  Y_data <- as.matrix(data[ , object$var_y], ncol = 1, nrow = NROW(data))
  X_index <- as.matrix(data[ , object$vars_index])
  while(j <= num_pred){
    split_info <- split_index(num_pred = num_pred, num_ind = j)
    ind_pos <- split_info$index_positions
    index <- split_info$index
    # Initialising index coefficients
    alpha_start <- init_alpha(Y = Y_data, X = X_index, index.ind = index, 
                              init.type = "penalisedReg", 
                              lambda0 = lambda0, lambda2 = lambda2, M = M)$alpha_init
    # Checking for all zero indices
    rm_ind <- numeric(1)
    for(i in 1:j){
      ind_alpha = alpha_start[ind_pos[[i]]]
      ind_zero = ifelse(ind_alpha == 0, TRUE, FALSE)
      if(all(ind_zero)){
        rm_ind <- i
        warning(paste0('Initialising ', paste0(j), '-index model: All coefficients of index', paste0(i), 
                       ' are zero. Reverting to the previous best model.')) 
        break
      }
    }
    if(rm_ind != 0){
      break
    }
    # Calculating indices
    ind <- vector(length = length(ind_pos), mode = "list")
    ind_names <- character(length = length(ind_pos))
    for(i in 1:length(ind)){
      ind[[i]] <- as.numeric(X_index[, ind_pos[[i]]] %*% 
                               as.matrix(alpha_start[ind_pos[[i]]], ncol = 1))
      ind_names[i] <- paste0("index", i)
    }
    names(ind) <- ind_names
    dat <- tibble::as_tibble(ind)
    dat_names <- colnames(dat)
    # Nonlinear function update 
    # Constructing the formula
    pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) %>%
      paste(collapse = "+") %>% 
      paste(object$var_y, "~", .)
    if (!is.null(object$vars_linear)) {
      pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    # Model fitting
    dat <- dplyr::bind_cols(data, dat)
    fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, method = "REML")
    # Derivatives of the fitted smooths
    dgz <- vector(length = length(dat_names), mode = "list")
    dgz_names <- character(length = length(dat_names))
    for (i in seq_along(dat_names)) {
      temp <- gratia::derivatives(fun1, type = "central", data = dat, 
                          term = paste0("s(", paste0(dat_names[i]), ")"))
      dgz[[i]] <- temp$derivative
      dgz_names[i] <- paste0("d", i)
    }
    names(dgz) <- dgz_names
    dgz <- as.matrix(as_tibble(dgz))
    # Constructing new initial values to have all predictors in each index
    new_init <- allpred_index(num_pred = num_pred, num_ind = j, ind_pos = ind_pos, 
                             alpha = alpha_start)
    index <- new_init$index
    ind_pos <- new_init$index_positions
    alpha <- new_init$alpha_init_new
    # Optimising the model
    best_alpha2 <- inner_update(x = fun1, data = data, yvar = object$var_y, 
                                index.vars = object$vars_index, 
                                linear.vars = object$vars_linear, 
                                num_ind = j, dgz = dgz, 
                                alpha_old = alpha, lambda0 = lambda0, 
                                lambda2 = lambda2, M = M, max.iter = max.iter, 
                                tol = tol, TimeLimit = TimeLimit, verbose = verbose)
    if(best_alpha2$min_loss >= MSE_current){
      print("MSE of the model is higher/equal to the previous model; reverting to the previous best model.")
      break
    }else if(length(best_alpha2$ind_pos) < j){
      alpha_current <- best_alpha2$best_alpha
      MSE_current <- best_alpha2$min_loss
      index_current <- best_alpha2$index.ind
      ind_pos_current <- best_alpha2$ind_pos
      X_new_current <- best_alpha2$X_new
      break
    }else{
      alpha_current <- best_alpha2$best_alpha
      MSE_current <- best_alpha2$min_loss
      index_current <- best_alpha2$index.ind
      ind_pos_current <- best_alpha2$ind_pos
      X_new_current <- best_alpha2$X_new
    }
    j = j + 1
  }
  # Calculating indices
  ind <- vector(length = length(ind_pos_current), mode = "list")
  ind_names <- character(length = length(ind_pos_current))
  for(i in 1:length(ind)){
    ind[[i]] <- as.numeric(X_new_current[, ind_pos_current[[i]]] %*% 
                             as.matrix(alpha_current[ind_pos_current[[i]]], ncol = 1))
    ind_names[i] <- paste0("index", i)
  }
  names(ind) <- ind_names
  dat <- tibble::as_tibble(ind)
  dat_names <- colnames(dat)
  ## Nonlinear function update 
  # Constructing the formula
  pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) %>%
    paste(collapse = "+") %>% 
    paste(object$var_y, "~", .)
  if (!is.null(object$vars_linear)) {
    pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
      paste(collapse = "+") %>% 
      paste(pre.formula, "+", .)
  }
  # Model fitting
  dat <- dplyr::bind_cols(data, dat)
  fun1_final <- mgcv::gam(as.formula(pre.formula), data = dat, method = "REML")
  print("Final model fitted!")
  final_smimodel <- make_smimodel(x = fun1_final, yvar = object$var_y, 
                                  index.vars = object$vars_index, 
                                  index.ind = index_current, 
                                  index.data = dat, index.names = dat_names,
                                  alpha = alpha_current, 
                                  linear.vars = object$vars_linear)
  return(final_smimodel)
}
