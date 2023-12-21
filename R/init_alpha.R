#' Initialising index coefficients
#'
#' Initialises index coefficient vector through linear regression or penalised
#' linear regression.
#'
#' @param Y Column matrix of response.
#' @param X Matrix of predictors entering indices.
#' @param index.ind An integer vector that assigns group index for each
#'   predictor.
#' @param init.type Type of initialisation for index coefficients.
#'   ("penalisedReg" - Penalised linear regression; "reg" - Linear regression)
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' 
#' @importFrom Matrix bdiag
#' @importFrom ROI Q_objective L_constraint OP V_bound ROI_solve 
#'
#' @export
init_alpha <- function(Y, X, index.ind, init.type = "penalisedReg", 
                       lambda0 = 1, lambda2 = 1, M = 10){
  p = NCOL(X)
  XtX = t(X) %*% X
  XtY = t(X) %*% Y
  if(init.type == "reg"){
    Q <- 2 * XtX
    L <- t(as.matrix(-2 * t(XtY), nrow = 1, ncol = p))
    # Objective function
    obj = ROI::Q_objective(Q = Q, L = L)
    # Optimization problem
    init <- ROI::OP(obj, types = rep("C", p),
               bounds = V_bound(li = 1:p, lb = rep(-Inf, p)),
               maximum = FALSE)
    # Solving
    sol <- ROI::ROI_solve(init, solver = "gurobi")
    # Initial alpha values
    alpha_raw <- sol$solution[1:p]
    alpha_0 <- unlist(tapply(sol$solution[1:p], index.ind, normalise_alpha))
  }else if(init.type == "penalisedReg"){
    Q <- as.matrix(2 * Matrix::bdiag(XtX + (lambda2*diag(p)), diag(0, p)))
    L <- t(as.matrix(c(-2 * t(XtY), rep(lambda0, p)), nrow = 1, ncol = 2*p))
    # Objective function
    obj = ROI::Q_objective(Q = Q, L = L)
    # Constraints
    cons1 <- cbind(diag(p), - M * diag(p))
    cons2 <- cbind(- diag(p), - M * diag(p))
    cons = ROI::L_constraint(
      L = rbind(cons1, cons2),
      dir = rep("<=", (2*p)),
      rhs = rep(0, 2*p)
    )
    # Optimization problem
    init <- ROI::OP(obj, cons, types = c(rep("C", p), rep("B", p)), 
               bounds = ROI::V_bound(li = 1:p, lb = rep(- M, p),
                                ui = 1:p, ub = rep(M, p),
                                nobj = 2*p), # lower default bound is 0
               maximum = FALSE)
    # Solving
    sol <- ROI::ROI_solve(init, solver = "gurobi")
    # Binary variables check
    coefs <- sol$solution[1:p]
    indicators <- sol$solution[(p+1):(p*2)]
    zero_pos <- which(indicators < 1e-5)
    coefs[zero_pos] <- 0
    # Initial alpha values
    alpha_raw <- coefs
    alpha_0 <- unlist(tapply(coefs, index.ind, normalise_alpha))
  }
  output <- list("alpha_init" = alpha_0, "alpha_nonNormalised" = alpha_raw)
  return(output)
}