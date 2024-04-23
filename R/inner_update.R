#' Updating index coefficients and non-linear functions iteratively
#'
#' Iteratively updates index coefficients and non-linear functions using mixed
#' integer programming. (Used within `update_smimodelFit()`; users are not expected
#' to directly call this function.)
#'
#' @param x Fitted `gam`.
#' @param data Training data set on which models will be trained. Should be a
#'   `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted individually (rather than considering as a
#'   part of an index considered in `index.vars`).
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model.
#' @param num_ind Number of indices.
#' @param dgz The `tibble` of derivatives of the estimated smooths.
#' @param alpha_old Current vector of index coefficients.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for loss.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.

inner_update <- function(x, data, yvar, family = gaussian(), index.vars, 
                         s.vars, linear.vars, num_ind, dgz, alpha_old, 
                         lambda0 = 1, lambda2 = 1, M = 10, max.iter = 50, 
                         tol = 0.001, TimeLimit = Inf,
                         MIPGap = 1e-4, NonConvex = -1, verbose = FALSE){
  data <- data %>%
    drop_na()
  data.Y <- as.matrix(data[ , yvar])
  X_index <- as.matrix(data[ , index.vars])
  num_pred <- NCOL(X_index)
  index <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    index[[i]] <- rep(i, num_pred)
  }
  index <- unlist(index)
  ind_pos <- split(1:length(index), index)
  Yhat <- as.matrix(x$fitted.values, ncol = 1)
  # Residuals (R) matrix
  R <- as.matrix(data.Y - Yhat)
  # Loss
  l2 <- loss(Y = as.matrix(data.Y), Yhat = Yhat, 
             alpha = alpha_old, 
             lambda0 = lambda0, lambda2 = lambda2)
  # Setting up a counter for increases in loss (If the loss increases for 
  # 3 consecutive iterations, the algorithm will terminate.)
  increase_count <- 0
  # Adjusting X (matrix of predictors) to fit number of indices
  X_new <- do.call(cbind, replicate(num_ind, X_index, simplify = FALSE))
  # Initial minimum loss and best index coefficient estimates
  best_l2 <- l2
  best_alpha <- alpha_old
  best_index <- index
  best_ind_pos <- ind_pos
  best_X_new <- X_new
  # Iteratively update index coefficients
  maxIt <- 1
  while(maxIt <= max.iter){
    alpha_new <- update_alpha(Y = R, X = X_new, num_pred = num_pred, 
                              num_ind = num_ind, index.ind = index, dgz = dgz, 
                              alpha_old = alpha_old, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, TimeLimit = TimeLimit,
                              MIPGap = MIPGap, NonConvex = NonConvex, 
                              verbose = verbose)
    if(all(alpha_new == 0)){
      best_l2 <- NULL
      best_alpha <- alpha_new
      best_index <- index
      best_ind_pos <- ind_pos
      best_X_new <- NULL
      print("Null indices are produced!")
      break
    }else{
      # Checking for all zero indices
      ind_rm_id <- numeric()
      ind_rm_pos <- numeric()
      for(i in 1:num_ind){
        if(all(alpha_new[ind_pos[[i]]] == 0)){
          ind_rm_id <- c(ind_rm_id, i)
          ind_rm_pos <- c(ind_rm_pos, ind_pos[[i]])
          message(paste0('Iteration ', maxIt, ': All coefficients of index', i,
                         ' are zero. Removing index', i, ' from subsequent iterations.')) 
        }
      }
      if(length(ind_rm_id) != 0){
        X_new <- as.matrix(X_new[ , -ind_rm_pos])
        alpha_new <- alpha_new[-ind_rm_pos]
        num_ind  <- num_ind - length(ind_rm_id)
        index <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          index[[i]] <- rep(i, num_pred)
        }
        index <- unlist(index)
        ind_pos <- split(1:length(index), index)
      }
      # Calculating indices
      ind <- vector(length = length(ind_pos), mode = "list")
      for(i in 1:length(ind)){
        ind[[i]] <- as.numeric(X_new[, ind_pos[[i]]] %*% 
                                 as.matrix(alpha_new[ind_pos[[i]]], ncol = 1))
      }
      dat_names <- names(ind) <- paste0("index", 1:length(ind))
      dat <- as_tibble(ind)
      # Nonlinear function update 
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
      dat <- dplyr::bind_cols(data, dat)
      fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, family = family, 
                        method = "REML")
      Yhat <- as.matrix(fun1$fitted.values, ncol = 1)
      # Derivatives of the fitted smooths
      dgz <- vector(length = length(dat_names), mode = "list")
      for (i in seq_along(dat_names)) {
        temp <- gratia::derivatives(fun1, type = "central", data = dat, 
                            select = paste0("s(", dat_names[i], ")"))
        dgz[[i]] <- temp$.derivative
      }
      names(dgz) <- paste0("d", seq_along(dat_names))
      dgz <- as.matrix(as_tibble(dgz))
      # Residuals (R) matrix
      R <- as.matrix(data.Y - Yhat)
      # Loss
      l2_new <- loss(Y = as.matrix(data.Y), Yhat = Yhat, 
                             alpha = alpha_new, 
                             lambda0 = lambda0, lambda2 = lambda2)
      eps <- (l2 - l2_new)/l2
      # Update loss and estimates
      alpha_old <- alpha_new
      l2 <- l2_new
      # Update minimum loss and best estimates
      if(l2_new < best_l2){
        best_l2 <- l2_new
        best_alpha <- alpha_new
        best_index <- index
        best_ind_pos <- ind_pos
        best_X_new <- X_new
      }
      # Check tolerance conditions
      if(eps <= 0){
        increase_count <- increase_count + 1
      }else{
        increase_count <- 0
      }
      if ((eps > 0) & (eps < tol)) { 
        print("Tolerance for loss reached!")
        break
      }else if(increase_count >= 3){
        print("Loss increased for 3 consecutive iterations!")
        break
      }
      maxIt <- maxIt + 1
      if (maxIt > max.iter) { 
        print("Maximum iterations reached!")
      }
    }
  }
  output <- list("best_alpha" = best_alpha, "min_loss" = best_l2, 
                 "index.ind" = best_index, "ind_pos" = best_ind_pos, 
                 "X_new" = best_X_new)
  return(output)
}