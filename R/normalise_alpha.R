#' Scaling index coefficient vectors to have unit norm
#'
#' Scales a coefficient vector of a particular index to have unit norm.
#'
#' @param alpha A vector of index coefficients.
#'
#' @export
normalise_alpha <- function (alpha) {
  anorm <- norm(matrix(alpha, ncol = 1))
  if (!(is.na(anorm) | all(alpha == 0))) 
    alpha <- alpha/anorm
  return(alpha)
}