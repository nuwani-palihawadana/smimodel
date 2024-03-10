#' Nonparametric Additive Model with Backward Elimination
#'
#' Fits a nonparametric additive model, with simultaneous variable selection
#' through a backward elimination procedure as proposed by Fan and Hyndman
#' (2012).
#'
#' Fits a nonparametric multiple index model to the data, with simultaneous
#' variable selection (hence "sparse").
#'
#' @param data The data set on which the model(s) will be trained. Must be a
#'   data set of class `tsibble`.
#' @param val.data Validation data set. (The data set on which the model
#'   selection will be performed.) Must be a data set of class `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param family A description of the error distribution and link function to be
#'   used in the model (see `glm` and `family`).
#' @param log.transformed Whether the response is log-transformed or not.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param s.vars A character vector of names of the predictor variables for
#'   which splines should be fitted (i.e. non-linear predictors). (default:
#'   NULL)
#' @param s.basedim Dimension of the bases used to represent the smooth terms
#'   corresponding to `s.vars`. (For more information refer `mgcv::s()`.)
#'   (default: NULL)
#' @param linear.vars A character vector of names of the predictor variables
#'   that should be included linearly into the model (i.e. linear predictors).
#'   (default: NULL)
#' @param tol Tolerance for the ratio of relative change in MAPE, used in model
#'   selection. (default: 0.001)
#' @param parallel Whether to use parallel computing in model selection or not.
#'   (default: FALSE)
#' @param workers If `parallel = TRUE`, number of workers to use. (default:
#'   NULL)
#'
#' @importFrom dplyr pull
#' @importFrom fabletools MAPE
#' @importFrom future multisession
#' @importFrom stats predict var
#' @importFrom tidyselect all_of
#'
#' @references Fan, S. & Hyndman, R.J. (2012). Short-Term Load Forecasting Based
#'   on a Semi-Parametric Additive Model. *IEEE Transactions on Power Systems*,
#'   27(1), 134-141. \url{http://doi.org/10.1109/TPWRS.2011.2162082}
#'
#' @export
model_backward <- function(data, val.data, yvar, family = gaussian(),
                     log.transformed = FALSE, neighbour = 0, 
                     s.vars = NULL, s.basedim = NULL, 
                     linear.vars = NULL,  tol = 0.001, 
                     parallel = FALSE, workers = NULL){
  if (!is_tsibble(data)) stop("data is not a tsibble.")
  if (!is_tsibble(val.data)) stop("val.data is not a tsibble.")
  if (is.null(c(linear.vars, s.vars))) 
    stop("No predictor variables specified; s.vars = NULL, linear.vars = NULL.")
  if(parallel){
    future::plan(multisession, workers = workers)
    map_f <- furrr::future_map
  } else {
    map_f <- purrr::map
  }
  # Data Preparation
  data1 <- data
  data2 <- val.data
  index_train <- index(data)
  if (length(key(data1)) == 0) {
    data1 <- data1 %>%
      mutate(dummy_key = rep(1, NROW(data1))) %>%
      as_tsibble(index = index_train, key = dummy_key)
    data2 <- data2 %>%
      mutate(dummy_key = rep(1, NROW(data2))) %>%
      as_tsibble(index = index_train, key = dummy_key)
  }
  key_train <- key(data1)[[1]]
  key_unique <- unique(as.character(sort(pull((data1[, {{ key_train }}])[, 1]))))
  key_num <- seq_along(key_unique)
  ref <- data.frame(key_unique, key_num)
  data1 <- data1 %>%
    mutate(
      num_key = as.numeric(factor(as.character({{ key_train }}), 
                                  levels = key_unique))
    )
  data2 <- data2 %>%
    mutate(
      num_key = as.numeric(factor(as.character({{ key_train }}), 
                                  levels = key_unique))
    )
  models_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)) {
    # Separating data for each element of the key
    df_cat <- data1 %>%
      filter((abs(num_key - ref$key_num[i]) <= neighbour) |
               (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
               (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    df_cat_val <- data2 %>%
      filter((abs(num_key - ref$key_num[i]) <= neighbour) |
               (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
               (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour))
    # Full Model
    allVars = c(s.vars, linear.vars)
    pre.formula <- paste0(yvar, " ~ ")
    for(d in seq_along(allVars)){
      if(allVars[d] %in% s.vars){
        if (!is.null(s.basedim)) {
          pre.formula <- paste0(
            pre.formula, "+s(", paste0(allVars[d]), ',bs="cr",k=',
            paste0(s.basedim), ")"
          )
        } else {
          pre.formula <- paste0(pre.formula, "+s(", paste0(allVars[d]), ',bs="cr")')
        }
      }else if(allVars[d] %in% linear.vars){
        pre.formula <- paste0(pre.formula, "+", paste0(allVars[d]))
      }
    }
    my.formula <- as.formula(pre.formula)
    # Model fitting
    model1 <- mgcv::gam(my.formula, family = family, method = "REML", 
                        data = df_cat)
    # Validation set MAPE
    # Predictions
    pred <- predict(object = model1, newdata = df_cat_val, type = "response")
    if(log.transformed == TRUE){
      # Back-transformation and bias-adjustment
      res_output <- model1$residuals
      pred_final <- as.vector((exp(pred)*(1+(var(res_output)/2))))
      # MAPE
      mape1 = MAPE(.resid = (exp(as.numeric(df_cat_val[,{{yvar}}])) - pred_final), 
                   .actual = exp(as.numeric(df_cat_val[,{{yvar}}])))
    }else{
      # MAPE
      mape1 = MAPE(.resid = (as.numeric(as.matrix(df_cat_val[,{{yvar}}], ncol = 1)) - pred), 
                   .actual = as.numeric(as.matrix(df_cat_val[,{{yvar}}])))
    }
    mape_old <- mape1
    mapeMinRatio <- 1
    Temp_s.vars <- s.vars
    Temp_linear.vars <- linear.vars
    # Testing reduced models
    while(mapeMinRatio > tol){
      allVars = c(Temp_s.vars, Temp_linear.vars)
      MAPE_list <- seq_along(allVars) %>%
        map_f(~ eliminate(ind = ., train = df_cat, val = df_cat_val, 
                          yvar = yvar, log.transformed = log.transformed, 
                          s.vars = Temp_s.vars, s.basedim = s.basedim, 
                          linear.vars = Temp_linear.vars))
      # Selecting best model
      best_model_pos <- which.min(unlist(MAPE_list))
      # MAPE check
      mape_new <- MAPE_list[[best_model_pos]]
      mapeMinRatio <- (mape_old - mape_new)/mape_old
      mape_old <- mape_new
      # Update Temp_s.vars and Temp_linear.vars
      if(mapeMinRatio >= 0){
        if(allVars[best_model_pos] %in% Temp_s.vars){
          Temp_s.vars <- Temp_s.vars[-best_model_pos]
        }else if(allVars[best_model_pos] %in% Temp_linear.vars){
          Temp_linear.vars <- Temp_linear.vars[-best_model_pos]
        }
      }
    }
    allVars = c(Temp_s.vars, Temp_linear.vars)
    pre.formula <- paste0(yvar, " ~ ")
    for(p in seq_along(allVars)){
      if(allVars[p] %in% Temp_s.vars){
        if (!is.null(s.basedim)) {
          pre.formula <- paste0(
            pre.formula, "+s(", paste0(allVars[p]), ',bs="cr",k=',
            paste0(s.basedim), ")"
          )
        } else {
          pre.formula <- paste0(pre.formula, "+s(", paste0(allVars[p]), ',bs="cr")')
        }
      }else if(allVars[p] %in% Temp_linear.vars){
        pre.formula <- paste0(pre.formula, "+", paste0(allVars[p]))
      }
    }
    my.formula <- as.formula(pre.formula)
    models_list[[i]] <- mgcv::gam(my.formula, family = family, method = "REML",
                                  data = df_cat)
    add <- df_cat %>%
      drop_na() %>%
      select({{ index_train }}, {{ key_train }})
    models_list[[i]]$model <- bind_cols(add, models_list[[i]]$model)
    models_list[[i]]$model <- as_tsibble(models_list[[i]]$model,
                                         index = index_train,
                                         key = all_of(key_train))
    print(paste0("Model ", i, " fitted!"))
  }
  data_list <- list(key_unique, models_list)
  models <- as_tibble(x = data_list, 
                      .rows = length(data_list[[1]]),
                      .name_repair = ~ vctrs::vec_as_names(..., 
                                                           repair = "universal", 
                                                           quiet = TRUE))
  models <- models %>%
    rename(key = ...1) %>%
    rename(fit = ...2)
  class(models) <- c("backward", "tbl_df", "tbl", "data.frame")
  return(models)
}
utils::globalVariables(c("...1", "...2"))