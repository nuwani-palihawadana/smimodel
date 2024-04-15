#' Updating a `smimodelFit`
#'
#' Optimises and updates a given `smimodelFit`.
#'
#' @param object A `smimodelFit` object.
#' @param data Training data set on which models will be trained. Should be a
#'   `tsibble`.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @param ... Other arguments not currently used.
#'
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom tibble as_tibble is_tibble
#' @importFrom tidyr drop_na

update_smimodelFit <- function(object, data, lambda0 = 1, lambda2 = 1, 
                            M = 10, max.iter = 50, 
                            tol = 0.001, tolCoefs = 0.001,
                            TimeLimit = Inf, MIPGap = 1e-4, 
                            NonConvex = -1, verbose = FALSE, ...){
  if (!tsibble::is_tsibble(data)) stop("data is not a tsibble.")
  data_index <- index(data)
  data_key <- key(data)[[1]]
  data <- data %>% drop_na()
  # Preparing inputs to `inner_update()`
  list_index <- object$alpha
  num_ind <- NCOL(list_index)
  alpha <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    alpha[[i]] <- list_index[ , i]
  }
  alpha <- unlist(alpha)
  names(alpha) <- NULL
  dgz <- as.matrix(object$derivatives)
  # Optimising the model
  best_alpha1 <- inner_update(x = object$gam, data = data, yvar = object$var_y,
                              family = object$gam$family$family,
                              index.vars = object$vars_index, 
                              s.vars = object$vars_s,
                              linear.vars = object$vars_linear, 
                              num_ind = num_ind, dgz = dgz, 
                              alpha_old = alpha, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, max.iter = max.iter, 
                              tol = tol, TimeLimit = TimeLimit, 
                              MIPGap = MIPGap, NonConvex = NonConvex, 
                              verbose = verbose)
  if(all(best_alpha1$best_alpha == 0)){
    # Constructing the formula and model fitting
    if (!is.null(object$vars_s) & !is.null(object$vars_linear)){
      pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) %>%
        paste(collapse = "+") %>% 
        paste(object$var_y, "~", .)
      pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(pre.formula, "+", .)
    }else if(!is.null(object$vars_s)){
      pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) %>%
        paste(collapse = "+") %>% 
        paste(object$var_y, "~", .)
    }else if(!is.null(object$vars_linear)){
      pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
        paste(collapse = "+") %>% 
        paste(object$var_y, "~", .)
    }else{
      pre.formula <- paste(object$var_y, "~", 1)
    }
    fun1_final <- mgcv::gam(as.formula(pre.formula), data = data, 
                            family = object$gam$family$family, method = "REML")
    add <- data %>%
      drop_na() %>%
      select({{ data_index }}, {{ data_key }})
    fun1_final$model <- bind_cols(add, fun1_final$model)
    fun1_final$model <- as_tsibble(fun1_final$model,
                                   index = data_index,
                                   key = all_of(data_key))
    index.names <- paste0("index", 1:length(best_alpha1$ind_pos))
    print("Final model fitted!")
    final_smimodel <- make_smimodelFit(x = fun1_final, yvar = object$var_y, 
                                       neighbour = object$neighbour,
                                       index.vars = object$vars_index, 
                                       index.ind = best_alpha1$index.ind, 
                                       index.data = NULL, index.names = index.names,
                                       alpha = best_alpha1$best_alpha, 
                                       s.vars = object$vars_s,
                                       linear.vars = object$vars_linear)
  }else{
    # Checking models with higher number of indices
    alpha_current <- best_alpha1$best_alpha
    loss_current <- best_alpha1$min_loss
    index_current <- best_alpha1$index.ind
    ind_pos_current <- best_alpha1$ind_pos
    X_new_current <- best_alpha1$X_new
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
        if (!is.null(object$vars_s)){
          pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ',bs="cr")')) %>%
            paste(collapse = "+") %>% 
            paste(pre.formula, "+", .)
        }
        if (!is.null(object$vars_linear)){
          pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
            paste(collapse = "+") %>% 
            paste(pre.formula, "+", .)
        }
        # Model fitting
        dat <- dplyr::bind_cols(data, dat)
        fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, 
                          family = object$gam$family$family, method = "REML")
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
        if (!is.null(object$vars_s)){
          pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ',bs="cr")')) %>%
            paste(collapse = "+") %>% 
            paste(pre.formula, "+", .)
        }
        if (!is.null(object$vars_linear)){
          pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
            paste(collapse = "+") %>% 
            paste(pre.formula, "+", .)
        }
        # Model fitting
        dat <- dplyr::bind_cols(data, dat)
        gam2 <- mgcv::gam(as.formula(pre.formula), data = dat, 
                          family = object$gam$family$family, method = "REML")
        # Derivatives of the fitted smooths
        dgz <- vector(length = length(dat_names), mode = "list")
        for (i in seq_along(dat_names)) {
          temp <- gratia::derivatives(gam2, type = "central", data = dat, 
                                      select = paste0("s(", paste0(dat_names[i]), ")"))
          dgz[[i]] <- temp$.derivative
        }
        names(dgz) <- paste0("d", seq_along(dat_names))
        dgz <- as.matrix(as_tibble(dgz))
        # Optimising the new model
        best_alpha2 <- inner_update(x = gam2, data = data, yvar = object$var_y,
                                    family = object$gam$family$family,
                                    index.vars = object$vars_index, 
                                    s.vars = object$vars_s,
                                    linear.vars = object$vars_linear, 
                                    num_ind = j, dgz = dgz, 
                                    alpha_old = alpha, lambda0 = lambda0, 
                                    lambda2 = lambda2, M = M, max.iter = max.iter,
                                    tol = tol, TimeLimit = TimeLimit, 
                                    MIPGap = MIPGap, NonConvex = NonConvex,
                                    verbose = verbose)
        if(all(best_alpha2$best_alpha == 0)){
          # Constructing the formula and model fitting
          if (!is.null(object$vars_s) & !is.null(object$vars_linear)){
            pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) %>%
              paste(collapse = "+") %>% 
              paste(object$var_y, "~", .)
            pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
              paste(collapse = "+") %>% 
              paste(pre.formula, "+", .)
          }else if(!is.null(object$vars_s)){
            pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ', bs="cr")')) %>%
              paste(collapse = "+") %>% 
              paste(object$var_y, "~", .)
          }else if(!is.null(object$vars_linear)){
            pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
              paste(collapse = "+") %>% 
              paste(object$var_y, "~", .)
          }else{
            pre.formula <- paste(object$var_y, "~", 1)
          }
          fun_null <- mgcv::gam(as.formula(pre.formula), data = data, 
                                family = object$gam$family$family, 
                                method = "REML")
          add <- data %>%
            drop_na() %>%
            select({{ data_index }}, {{ data_key }})
          fun_null$model <- bind_cols(add, fun_null$model)
          fun_null$model <- as_tsibble(fun_null$model,
                                       index = data_index,
                                       key = all_of(data_key))
          Y <- as.matrix(data[ , object$var_y])
          Yhat <- as.matrix(fun_null$fitted.values, 
                            ncol = 1, nrow = length(fun_null$fitted.values))
          loss_null <- loss(Y = Y, Yhat = Yhat, alpha = best_alpha2$best_alpha,
                            lambda0 = lambda0, lambda2 = lambda2)
          if(loss_null >= loss_current){
            break
          }else{
            index.names <- paste0("index", 1:length(best_alpha2$ind_pos))
            alpha_current <- best_alpha2$best_alpha
            print("Final model fitted!")
            final_smimodel <- make_smimodelFit(x = fun_null, 
                                               yvar = object$var_y, 
                                               neighbour = object$neighbour,
                                               index.vars = object$vars_index, 
                                               index.ind = best_alpha2$index.ind, 
                                               index.data = NULL, 
                                               index.names = index.names,
                                               alpha = best_alpha2$best_alpha, 
                                               s.vars = object$vars_s,
                                               linear.vars = object$vars_linear)
            break
          }
        }else{
          # Termination/continuation checks
          if((length(ind_pos_current) == length(best_alpha2$ind_pos))){
            comparison <- abs(alpha_current - best_alpha2$best_alpha) <= tolCoefs
            if(all(comparison)){
              if(best_alpha2$min_loss >= loss_current){
                break
              }else{
                alpha_current <- best_alpha2$best_alpha
                loss_current <- best_alpha2$min_loss
                index_current <- best_alpha2$index.ind
                ind_pos_current <- best_alpha2$ind_pos
                X_new_current <- best_alpha2$X_new
                break
              }
            }
            if(best_alpha2$min_loss >= loss_current){
              break
            }else{
              alpha_current <- best_alpha2$best_alpha
              loss_current <- best_alpha2$min_loss
              index_current <- best_alpha2$index.ind
              ind_pos_current <- best_alpha2$ind_pos
              X_new_current <- best_alpha2$X_new
            }
          }else{
            if(best_alpha2$min_loss >= loss_current){
              break
            }else{
              alpha_current <- best_alpha2$best_alpha
              loss_current <- best_alpha2$min_loss
              index_current <- best_alpha2$index.ind
              ind_pos_current <- best_alpha2$ind_pos
              X_new_current <- best_alpha2$X_new
            }
          }
          j <- length(ind_pos_current)+1
        }
      }
    }
    if(any(alpha_current != 0)){
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
      if (!is.null(object$vars_s)){
        pre.formula <- lapply(object$vars_s, function(var) paste0("s(", var, ',bs="cr")')) %>%
          paste(collapse = "+") %>% 
          paste(pre.formula, "+", .)
      }
      if (!is.null(object$vars_linear)){
        pre.formula <- lapply(object$vars_linear, function(var) paste0(var)) %>%
          paste(collapse = "+") %>% 
          paste(pre.formula, "+", .)
      }
      dat <- dplyr::bind_cols(data, dat)
      # Model fitting
      fun1_final <- mgcv::gam(as.formula(pre.formula), data = dat, 
                              family = object$gam$family$family, method = "REML")
      add <- data %>%
        drop_na() %>%
        select({{ data_index }}, {{ data_key }})
      fun1_final$model <- bind_cols(add, fun1_final$model)
      fun1_final$model <- as_tsibble(fun1_final$model,
                                     index = data_index,
                                     key = all_of(data_key))
      print("Final model fitted!")
      final_smimodel <- make_smimodelFit(x = fun1_final, 
                                         yvar = object$var_y,
                                         neighbour = object$neighbour,
                                         index.vars = object$vars_index, 
                                         index.ind = index_current, 
                                         index.data = dat, 
                                         index.names = dat_names,
                                         alpha = alpha_current, 
                                         s.vars = object$vars_s,
                                         linear.vars = object$vars_linear)
    }
  }
  return(final_smimodel)
}
