#' Updating a `smimodel` - Alternative function (trial for changes in algorithm) 
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
#' @param tolCoefs Tolerance for coefficients. 
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param verbose The option to print detailed solver output.
#' @param ... Other arguments not currently used.
#' 
#' @importFrom dplyr bind_rows
#'
#' @export
alt_update_smimodel <- function(object, data, lambda0 = 1, lambda2 = 1, 
                            M = 10, max.iter = 50, 
                            tol = 0.001, tolCoefs = 0.001,
                            TimeLimit = Inf, verbose = FALSE, ...){
  if (!tibble::is_tibble(data)) stop("data is not a tibble.")
  data <- data %>% drop_na()
  # Preparing inputs to `inner_update()`
  list_index <- object[1:(length(object)-4)]
  num_ind <- length(list_index)
  alpha <- vector(mode = "list", length = num_ind)
  dgz <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    alpha[[i]] <- list_index[[i]]$coefficients
    dgz[[i]] <- list_index[[i]]$derivatives
  }
  alpha <- unlist(alpha)
  names(dgz) <- paste0("d", 1:num_ind)
  dgz <- as.matrix(tibble::as_tibble(dgz))
  # Optimising the model
  best_alpha1 <- alt_inner_update(x = object$gam, data = data, yvar = object$var_y,
                              index.vars = object$vars_index, 
                              linear.vars = object$vars_linear, 
                              num_ind = num_ind, dgz = dgz, 
                              alpha_old = alpha, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, max.iter = max.iter, 
                              tol = tol, TimeLimit = TimeLimit, verbose = verbose)
  # Checking models with higher number of indices
  alpha_current <- best_alpha1$best_alpha
  loss_current <- best_alpha1$min_loss
  index_current <- best_alpha1$index.ind
  ind_pos_current <- best_alpha1$ind_pos
  X_new_current <- best_alpha1$X_new
  if(all(alpha_current == 0)){
    # Constructing the formula and model fitting
    if (!is.null(object$vars_linear)){
      pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(object$var_y, "~", .)
    }else{
      pre.formula <- paste(object$var_y, "~", 1)
    }
    fun1_final <- mgcv::gam(as.formula(pre.formula), data = data, method = "REML")
    index.names <- paste0("index", 1:length(ind_pos_current))
    print("Final model fitted!")
    final_smimodel <- make_smimodel(x = fun1_final, yvar = object$var_y, 
                                    index.vars = object$vars_index, 
                                    index.ind = index_current, 
                                    index.data = NULL, index.names = index.names,
                                    alpha = alpha_current, 
                                    linear.vars = object$vars_linear)
  }else{
    j <- length(ind_pos_current)+1
    num_pred <- length(object$vars_index)
    while(j <= num_pred){
      # Checking whether a new index can be added without variable repetition;
      # if not, terminating the loop.
      index_list <- split(alpha_current, index_current)
      index_mat <- do.call(rbind, index_list)
      drop_pred_ind <- which(colSums(index_mat) == 0)
      if(length(drop_pred_ind) == 0){
        break
      }else{
        # Initialising the new index to be added to the current model
        drop_pred_name <- object$vars_index[drop_pred_ind]
        X_init <- as.matrix(data[, drop_pred_name])
        index.ind <- rep(1, length(drop_pred_name))
        # Calculating indices
        ind <- vector(length = length(ind_pos_current), mode = "list")
        for(i in 1:length(ind)){
          ind[[i]] <- as.numeric(X_new_current[, ind_pos_current[[i]]] %*% 
                                   as.matrix(alpha_current[ind_pos_current[[i]]], ncol = 1))
        }
        dat_names <- names(ind) <- paste0("index", 1:length(ind))
        dat <- as_tibble(ind)
        # Nonlinear function update 
        # Constructing the formula
        pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) %>%
          paste(collapse = "+") %>% 
          paste(object$var_y, "~", .)
        if (!is.null(object$vars_linear)){
          pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
            paste(collapse = "+") %>% 
            paste(pre.formula, "+", .)
        }
        # Model fitting
        dat <- dplyr::bind_cols(data, dat)
        fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, method = "REML")
        Yhat <- as.matrix(fun1$fitted.values, ncol = 1, nrow = length(fun1$fitted.values))
        # Residuals (R) matrix
        R <- as.matrix(data[ , object$var_y] - Yhat)
        coefs_init <- init_alpha(Y = R, X = X_init, index.ind = index.ind,
                                 init.type = "reg")$alpha_init
        new_ind <- numeric(length = num_pred)
        new_ind[drop_pred_ind] <- coefs_init
        alpha <- c(alpha_current, new_ind)
        X_index <- data[ , object$vars_index]
        # Adjusting X (matrix of predictors) to fit number of indices
        X_new <- as.matrix(do.call(cbind, replicate(j, X_index, simplify = FALSE)))
        index <- vector(mode = "list", length = j)
        for(i in 1:j){
          index[[i]] <- rep(i, num_pred)
        }
        index <- unlist(index)
        ind_pos <- split(1:length(index), index)
        # Calculating indices
        ind <- vector(length = length(ind_pos), mode = "list")
        for(i in 1:length(ind)){
          ind[[i]] <- as.numeric(X_new[, ind_pos[[i]]] %*% 
                                   as.matrix(alpha[ind_pos[[i]]], ncol = 1))
        }
        dat_names <- names(ind) <- paste0("index", 1:length(ind))
        dat <- as_tibble(ind)
        # Nonlinear function update 
        # Constructing the formula
        pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) %>%
          paste(collapse = "+") %>% 
          paste(object$var_y, "~", .)
        if (!is.null(object$vars_linear)){
          pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
            paste(collapse = "+") %>% 
            paste(pre.formula, "+", .)
        }
        # Model fitting
        dat <- dplyr::bind_cols(data, dat)
        gam2 <- mgcv::gam(as.formula(pre.formula), data = dat, method = "REML")
        # Derivatives of the fitted smooths
        dgz <- vector(length = length(dat_names), mode = "list")
        for (i in seq_along(dat_names)) {
          temp <- gratia::derivatives(gam2, type = "central", data = dat, 
                                      term = paste0("s(", paste0(dat_names[i]), ")"))
          dgz[[i]] <- temp$derivative
        }
        names(dgz) <- paste0("d", seq_along(dat_names))
        dgz <- as.matrix(as_tibble(dgz))
        # Optimising the new model
        best_alpha2 <- alt_inner_update(x = gam2, data = data, yvar = object$var_y, 
                                        index.vars = object$vars_index, 
                                        linear.vars = object$vars_linear, 
                                        num_ind = j, dgz = dgz, 
                                        alpha_old = alpha, lambda0 = lambda0, 
                                        lambda2 = lambda2, M = M, max.iter = max.iter,
                                        tol = tol, TimeLimit = TimeLimit, verbose = verbose)
        # Termination/continuation checks
        if(best_alpha2$min_loss >= loss_current){
          break
        }else{
          alpha_current <- best_alpha2$best_alpha
          loss_current <- best_alpha2$min_loss
          index_current <- best_alpha2$index.ind
          ind_pos_current <- best_alpha2$ind_pos
          X_new_current <- best_alpha2$X_new
        }
        j <- length(ind_pos_current)+1
      }
    }
    # Calculating indices
    ind <- vector(length = length(ind_pos_current), mode = "list")
    for(i in 1:length(ind)){
      ind[[i]] <- as.numeric(X_new_current[, ind_pos_current[[i]]] %*% 
                               as.matrix(alpha_current[ind_pos_current[[i]]], ncol = 1))
    }
    dat_names <- names(ind) <- paste0("index", 1:length(ind))
    dat <- tibble::as_tibble(ind)
    ## Nonlinear function update 
    # Constructing the formula
    pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) %>%
      paste(collapse = "+") %>% 
      paste(object$var_y, "~", .)
    if (!is.null(object$vars_linear)){
      pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }
    dat <- dplyr::bind_cols(data, dat)
    # Model fitting
    fun1_final <- mgcv::gam(as.formula(pre.formula), data = dat, method = "REML")
    print("Final model fitted!")
    final_smimodel <- make_smimodel(x = fun1_final, yvar = object$var_y, 
                                    index.vars = object$vars_index, 
                                    index.ind = index_current, 
                                    index.data = dat, index.names = dat_names,
                                    alpha = alpha_current, 
                                    linear.vars = object$vars_linear)
  }
  return(final_smimodel)
}
