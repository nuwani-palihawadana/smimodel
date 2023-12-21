#' Updating index coefficients and non-linear functions iteratively
#'
#' Iteratively updates index coefficients and non-linear functions using mixed
#' integer programming. (Used within `update_smimodel()`; users are not expected
#' to directly call this function.)
#'
#' @param x Fitted `gam`
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param yvar Name of the response variable as a character string.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
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
#' @param verbose The option to print detailed solver output.
#'
#' @export
inner_update <- function(x, data, yvar, index.vars, linear.vars, 
                         num_ind, dgz, alpha_old, lambda0 = 1, lambda2 = 1, 
                         M = 10, max.iter = 50, tol = 0.001, TimeLimit = Inf,
                         verbose = FALSE){
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
  Yhat <- as.matrix(x$fitted.values, ncol = 1, nrow = length(x$fitted.values))
  # Residuals (R) matrix
  R <- as.matrix(data.Y - Yhat)
  # Loss
  l2 <- LossFunction(Y = as.matrix(data.Y), Yhat = Yhat, 
                     alpha = alpha_old, 
                     lambda0 = lambda0, lambda2 = lambda2)
  # Initial minimum loss and best index coefficient estimates
  min_l2 <- l2
  best_alpha <- alpha_old
  # Setting up a counter for increases in loss (If the loss increases for 
  # 3 consecutive iterations, the algorithm will terminate.)
  increase_count <- 0
  # Setting up a counter for same results repeating (If the same loss values 
  # appears for 3 times, the algorithm will terminate.)
  similar_count <- 1
  # Adjusting X (matrix of predictors) to fit number of indices
  X_new <- do.call(cbind, replicate(num_ind, X_index, simplify = FALSE))
  # Iteratively update index coefficients
  maxIt <- 1
  while(maxIt <= max.iter){
    alpha_new <- update_alpha(Y = R, X = X_new, num_pred = num_pred, 
                              num_ind = num_ind, index.ind = index, dgz = dgz, 
                              alpha_old = alpha_old, lambda0 = lambda0, 
                              lambda2 = lambda2, M = M, TimeLimit = TimeLimit,
                              verbose = verbose)
    ##### TO DO: You must remove this if condition! ############################
    if(is.null(alpha_new)){
      print("TimeLimit exceeded; intial estimates are chosen!")
      break
    }else{
      # Checking for all zero indices
      ind_rm_id <- numeric()
      ind_rm_pos <- numeric()
      for(i in 1:num_ind){
        ind_alpha = alpha_new[ind_pos[[i]]]
        ind_zero = ifelse(ind_alpha == 0, TRUE, FALSE)
        if(all(ind_zero)){
          ind_rm_id <- c(ind_rm_id, i)
          ind_rm_pos <- c(ind_rm_pos, ind_pos[[i]])
          warning(paste0('Iteration ', paste0(maxIt), ': All coefficients of index', paste0(i), 
                         ' are zero. Removing index', paste0(i), ' from subsequent iterations.')) 
        }
      }
      if(length(ind_rm_id) != 0){
        X_new = as.matrix(X_new[ , -ind_rm_pos])
        alpha_new = alpha_new[-ind_rm_pos]
        num_ind  = num_ind - length(ind_rm_id)
        index <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          index[[i]] <- rep(i, num_pred)
        }
        index <- unlist(index)
        ind_pos <- split(1:length(index), index)
      }
      # Calculating indices
      ind <- vector(length = length(ind_pos), mode = "list")
      ind_names <- character(length = length(ind_pos))
      for(i in 1:length(ind)){
        ind[[i]] <- as.numeric(X_new[, ind_pos[[i]]] %*% 
                                 as.matrix(alpha_new[ind_pos[[i]]], ncol = 1))
        ind_names[i] <- paste0("index", i)
      }
      names(ind) <- ind_names
      dat <- as_tibble(ind)
      dat_names <- colnames(dat)
      # Nonlinear function update 
      # Constructing the formula
      pre.formula <- lapply(dat_names, function(var) paste0("s(", var, ',bs="cr")')) %>%
        paste(collapse = "+") %>% 
        paste(yvar, "~", .)
      if (!is.null(linear.vars)) {
        pre.formula <- lapply(linear.vars, function(var) paste0(var)) %>%
          paste(collapse = "+") %>% 
          paste(pre.formula, "+", .)
      }
      # Model fitting
      dat <- dplyr::bind_cols(data, dat)
      fun1 <- mgcv::gam(as.formula(pre.formula), data = dat, method = "REML")
      Yhat <- as.matrix(fun1$fitted.values, ncol = 1, 
                        nrow = length(fun1$fitted.values))
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
      # Residuals (R) matrix
      R <- as.matrix(data.Y - Yhat)
      # Loss
      l2_new <- LossFunction(Y = as.matrix(data.Y), Yhat = Yhat, 
                             alpha = alpha_new, 
                             lambda0 = lambda0, lambda2 = lambda2)
      eps <- (l2 - l2_new)/l2
      alpha_old <- alpha_new
      l2 <- l2_new
      if(all(eps < 0)){
        increase_count = increase_count + 1
      }else{
        increase_count = 0
      }
      if(l2_new == min_l2){
        similar_count = similar_count + 1
      }else{
        similar_count = 1
      }
      # Update minimum loss and best estimates
      if(l2_new < min_l2){
        min_l2 <- l2_new
        best_alpha <- alpha_new
      }
      if (all(eps >= 0) & all(eps < tol)) { 
        print("Tolerance reached!")
        break
      }else if(increase_count >= 2){
        print("Alternative termination condition 1 reached!")
        break
      }else if(similar_count >= 3){
        print("Alternative termination condition 2 reached!")
        break
      }
      maxIt <- maxIt + 1
      if (maxIt > max.iter) { 
        print("Maximum iterations reached!")
      }
    }
  }
  output <- list("best_alpha" = best_alpha, "min_loss" = min_l2, 
                 "index.ind" = index, "ind_pos" = ind_pos, "X_new" = X_new)
  return(output)
}