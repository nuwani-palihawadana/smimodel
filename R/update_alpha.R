#' Updating index coefficients using MIP
#'
#' Updates index coefficients by solving a mixed integer program.
#'
#' @param Y Column matrix of response.
#' @param X Matrix of predictors (size adjusted to number of indices).
#' @param num_pred Number of predictors.
#' @param num_ind Number of indices.
#' @param index.ind An integer vector that assigns group index for each
#'   predictor.
#' @param dgz The `tibble` of derivatives of the estimated smooths from previous
#'   iteration.
#' @param alpha_old Vector of index coefficients from previous iteration.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param verbose The option to print detailed solver output.
#'
#' @export
update_alpha <- function(Y, X, num_pred, num_ind, index.ind, dgz, alpha_old, 
                         lambda0 = 1, lambda2 = 1, M = 10, TimeLimit = Inf,
                         MIPGap = 1e-4, verbose = FALSE){
  dgz <- as.matrix(dgz[, index.ind])
  V <- X * dgz
  print("Time (in seconds) taken for the calculation of VtV:")
  time <- system.time(VtV <- t(V) %*% V)
  print(time)
  VtR <- t(V) %*% Y
  VtR_t <- t(VtR)
  VtV_alpha_old_t <- t(VtV %*% alpha_old)
  Q <- as.matrix(2 * Matrix::bdiag(VtV + lambda2*diag(num_pred*num_ind), 
                                   diag(0, num_pred*num_ind)))
  L <- t(as.matrix(c((-2*VtV_alpha_old_t - 2*VtR_t),
                     rep(lambda0, num_ind*num_pred)), 
                   ncol = 2*num_ind*num_pred, nrow = 1))
  # Objective function
  print("Time (in seconds) taken for the construction of the quadratic objective function:")
  time <- system.time(obj <- ROI::Q_objective(Q = Q, L = L))
  print(time)
  # Constraints
  C1 <- Matrix::bdiag(replicate(num_ind, diag(num_pred), simplify = FALSE))
  C2 <- Matrix::bdiag(replicate(num_ind, diag(-M, num_pred), simplify = FALSE))
  cons1 <- as.matrix(cbind(C1, C2))
  C3 <- Matrix::bdiag(replicate(num_ind, diag(-1, num_pred), simplify = FALSE))
  C4 <- Matrix::bdiag(replicate(num_ind, diag(-M, num_pred), simplify = FALSE))
  cons2 <- as.matrix(cbind(C3, C4))
  if(num_ind > 1){
    C5 <- do.call(cbind, replicate(num_ind, diag(0, num_pred), simplify = FALSE))
    C6 <- do.call(cbind, replicate(num_ind, diag(num_pred), simplify = FALSE))
    cons3 <- as.matrix(cbind(C5, C6))
    cons <- ROI::L_constraint(
      L = rbind(cons1, cons2, cons3),
      dir = c(rep("<=", (2*num_ind*num_pred)), rep("<=", num_pred)),
      rhs = c(rep(0, (2*num_ind*num_pred)), rep(1, num_pred))
    )
  }else{
    cons <- ROI::L_constraint(
      L = rbind(cons1, cons2),
      dir = rep("<=", (2*num_ind*num_pred)),
      rhs = rep(0, (2*num_ind*num_pred))
    )
  }
  # Optimization problem
  init <- ROI::OP(obj, cons, types = c(rep("C", num_ind*num_pred), rep("B", num_ind*num_pred)),
                  bounds = ROI::V_bound(li = 1:(num_ind*num_pred), lb = rep(-M, num_ind*num_pred),
                                        ui = 1:(num_ind*num_pred), ub = rep(M, num_ind*num_pred),
                                        nobj = 2*num_ind*num_pred), # lower default bound is 0
                  maximum = FALSE)
  # Solving
  print("Time (in seconds) taken for solving the MIP:")
  time <- system.time(sol <- ROI::ROI_solve(init, solver = "gurobi", TimeLimit = TimeLimit, 
                        MIPGap = MIPGap, verbose = verbose))
  print(time)
  # Binary variables check
  coefs <- sol$solution[1:(num_ind*num_pred)]
  indicators <- sol$solution[((num_ind*num_pred)+1):(num_ind*num_pred*2)]
  zero_pos <- which(indicators < 1e-5)
  coefs[zero_pos] <- 0
  # Index coefficients
  alpha <- unlist(tapply(coefs, index.ind, normalise_alpha))
  return(alpha)
}