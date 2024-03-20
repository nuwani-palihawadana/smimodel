#' Calculating the loss in MIP
#'
#' Calculates the value of the objective function (loss function) of the
#' corresponding mixed integer program.
#'
#' @param Y Column matrix of response.
#' @param Yhat Predicted value of the response.
#' @param alpha Vector of index coefficients.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.

loss <- function(Y, Yhat, alpha, lambda0, lambda2){
  sqdError <- sum((Y - Yhat)^2)
  L0Penalty <- lambda0*sum(alpha != 0)
  nonzeroAlpha <- alpha[alpha != 0]
  L2Penalty <- lambda2*sum(nonzeroAlpha^2)
  return((sqdError + L0Penalty + L2Penalty))
}
