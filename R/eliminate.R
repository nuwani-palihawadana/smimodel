#' Function to eliminate a specified variable and fit a nonparametric additive
#' model with remaining variables
#'
#' This is an internal function of the package `smimodel`, and designed to be
#' called from `model_backward()`.
#'
#' @param ind An integer corresponding to the position of the predictor variable
#'   to be eliminated when fitting the model. (i.e. the function will combine
#'   `s.vars` and `linear.vars` in a single vector and eliminate the element
#'   corresponding to `ind`.)
#' @param train The data set on which the model(s) will be trained. Must be a
#'   data set of class `tsibble`.
#' @param val Validation data set. (The data set on which the model selection
#'   will be performed.) Must be a data set of class `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see `glm` and `family`).
#' @param log.transformed Whether the response is log-transformed or not.
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted (i.e. non-linear predictors). (default:
#'   NULL)
#' @param s.basedim Dimension of the bases used to represent the smooth terms
#'   corresponding to `s.vars`. (For more information refer `mgcv::s()`.)
#'   (default: NULL)
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model (i.e. linear predictors).
#'   (default: NULL)

eliminate <- function(ind, train, val, yvar, family = gaussian(), 
                      log.transformed = FALSE, s.vars = NULL, s.basedim = NULL, 
                      linear.vars = NULL){
  allVars = c(s.vars, linear.vars)
  pre.formula <- paste0(yvar, " ~ ")
  temp.var1 <- allVars[-ind]
  for(j in seq_along(temp.var1)){
    if(temp.var1[j] %in% s.vars){
      if (!is.null(s.basedim)) {
        pre.formula <- paste0(
          pre.formula, "+s(", paste0(temp.var1[j]), ',bs="cr",k=',
          paste0(s.basedim), ")"
        )
      } else {
        pre.formula <- paste0(pre.formula, "+s(", paste0(temp.var1[j]), ',bs="cr")')
      }
    }else if(temp.var1[j] %in% linear.vars){
      pre.formula <- paste0(pre.formula, "+", paste0(temp.var1[j]))
    }
  }
  my.formula <- as.formula(pre.formula)
  # Model fitting
  model1 <- mgcv::bam(my.formula, family = family, method = "REML", data = train)
  # Validation set MAPE
  # Predictions
  pred <- predict(object = model1, newdata = val, type = "response")
  if(log.transformed == TRUE){
    # Back-transformation and bias-adjustment
    res_output <- model1$residuals
    pred_final <- as.vector((exp(pred)*(1+(var(res_output)/2))))
    # MAPE
    mape1 = MAPE(.resid = (exp(as.numeric(val[,{{yvar}}])) - pred_final), 
                 .actual = exp(as.numeric(val[,{{yvar}}])))
  }else{
    # MAPE
    mape1 = MAPE(.resid = (as.numeric(as.matrix(val[,{{yvar}}], ncol = 1)) - pred), 
                 .actual = as.numeric(as.matrix(val[,{{yvar}}])))
  }
  return(mape1)
}