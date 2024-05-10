#' @importFrom stats ts
#' 
#' @method forecast gaimFit
#'
#' @export
forecast.gaimFit <- function(object, newdata, 
                             recursive = FALSE, recursive_colRange = NULL, ...){
  method <- "Groupwise Additive Index Model"
  pred <- predict(object = object, newdata = newdata, recursive = recursive,
                  recursive_colRange = recursive_colRange)$.predict
  fitted <- augment(object)$.fitted
  residuals <- augment(object)$.resid
  return(structure(list(method = method, 
                        model = object, 
                        mean = ts(pred), 
                        fitted = ts(fitted), 
                        residuals = ts(residuals)), 
                   class = "forecast"))
}
#' @export
generics::forecast