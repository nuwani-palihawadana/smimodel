#' Conformal bootstrap prediction intervals through time series cross-validation
#' forecasting
#'
#' Compute prediction intervals by applying the conformal bootstrap method to
#' subsets of time series data using a rolling forecast origin.
#'
#' @param object Fitted model object of class \code{smimodel}, \code{backward},
#'   \code{gaimFit} or \code{pprFit}.
#' @param data Data set. Must be a data set of class \code{tsibble}.(Make sure
#'   there are no additional date or time related variables except for the
#'   \code{index} of the \code{tsibble}). If multiple models are fitted, the
#'   grouping variable should be the \code{key} of the \code{tsibble}. If a
#'   \code{key} is not specified, a dummy key with only one level will be
#'   created.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an \code{integer}. If \code{neighbour =
#'   x}, \code{x} number of keys before the key of interest and \code{x} number
#'   of keys after the key of interest are grouped together for model fitting.
#'   The default is \code{neighbour = 0} (i.e. no neighbours are considered for
#'   model fitting).
#' @param predictor.vars A character vector of names of the predictor variables.
#' @param h Forecast horizon.
#' @param ncal Length of a calibration window.
#' @param num.futures Number of possible future sample paths to be generated in
#'   bootstrap.
#' @param level Confidence level for prediction intervals.
#' @param forward If \code{TRUE}, the final forecast origin for forecasting is
#'   \eqn{y_T}. Otherwise, the final forecast origin is \eqn{y_{T-1}}.
#' @param initial Initial period of the time series where no cross-validation
#'   forecasting is performed.
#' @param window Length of the rolling window. If \code{NULL}, a rolling window
#'   will not be used.
#' @param roll.length Number of observations by which each rolling/expanding
#'   window should be rolled forward.
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable is treated non-linearly in the
#'   estimated model, will be truncated to be in the in-sample range before
#'   obtaining predictions. If any variables are listed here will be excluded
#'   from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colNames If \code{recursive = TRUE}, a character vector
#'   giving the names of the columns in test data to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{data}, with no break in the lagged variable sequence even if
#'   some of the intermediate lags are not used as predictors.
#' @param na.rm logical; if \code{TRUE} (default), any \code{NA} and
#'   \code{NaN}'s are removed from the sample before the quantiles are computed.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{cb_cvforecast}, which is a list that
#'   contains following elements: \item{x}{The original time series.}
#'   \item{method}{A character string "cb_cvforecast".}
#'   \item{fit_times}{The number of times the model is fitted in
#'   cross-validation.} \item{mean}{Point forecasts as a multivariate time
#'   series, where the \eqn{h^{th}} column holds the point forecasts for
#'   forecast horizon \eqn{h}. The time index corresponds to the period for
#'   which the forecast is produced.}
#'   \item{error}{Forecast errors given by \eqn{e_{t+h|t} = y_{t+h} -
#'   \hat{y}_{t+h|t}}.} \item{level}{The confidence values associated with the
#'   prediction intervals.}
#'   \item{lower}{A list containing lower bounds for prediction intervals for
#'   each level. Each element within the list will be a multivariate time series
#'    with the same dimensional characteristics as \code{mean}.}
#'    \item{upper}{A list containing upper bounds for prediction intervals for
#'    each level. Each element within the list will be a multivariate time
#'    series with the same dimensional characteristics as \code{mean}.}
#'
#' @seealso \code{\link{bb_cvforecast}}
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1105
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 +
#'     (0.35*x_lag_002 + 0.7*x_lag_005)^2 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1100, ]
#'
#' # Model fitting
#' smimodel_ppr <- model_smimodel(data = sim_train,
#'                                yvar = "y",
#'                                index.vars = index.vars,
#'                                initialise = "ppr")
#'
#' # Conformal bootstrap prediction intervals (3-steps-ahead interval forecasts)
#' set.seed(12345)
#' smimodel_ppr_cb <- cb_cvforecast(object = smimodel_ppr,
#'                                  data = sim_data,
#'                                  yvar = "y",
#'                                  predictor.vars = index.vars,
#'                                  h = 3,
#'                                  ncal = 30,
#'                                  num.futures = 100,
#'                                  window = 1000)
#' }
#'
#' @export
cb_cvforecast <- function(object, data, yvar, neighbour = 0, predictor.vars,
                          h = 1, ncal = 100, num.futures = 1000,
                          level = c(80, 95), forward = TRUE,
                          initial = 1, window = NULL, roll.length = 1,
                          exclude.trunc = NULL,
                          recursive = FALSE, recursive_colNames = NULL, 
                          na.rm = TRUE, ...) {
  # Check input data
  if (!is_tsibble(data)) stop("data is not a tsibble.")

  index_data <- index(data)
  if (length(key(data)) == 0) {
    data <- data |>
      dplyr::mutate(dummy_key = rep(1, NROW(data))) |>
      tsibble::as_tsibble(index = index_data, key = dummy_key)
  }
  key_data1 <- key(data)[[1]]

  data1 <- data |>
    as_tibble() |>
    arrange({{index_data}})
  y <- as.ts(data1[ , yvar][[1]])
  n <- length(y) - forward * h

  # Check confidence levels
  if (min(level) > 0 && max(level) < 1) {
    level <- 100 * level
  } else if (min(level) < 0 || max(level) > 99.99) {
    stop("confidence limit out of range")
  }
  level <- sort(level)

  # Check other inputs
  if (h <= 0)
    stop("forecast horizon out of bounds")
  if (initial < 1 | initial > n)
    stop("initial period out of bounds")
  if (initial == n && !forward)
    stop("initial period out of bounds")
  if (!is.null(window)) {
    if (window < 1 | window > n)
      stop("window out of bounds")
    if (window == n && !forward)
      stop("window out of bounds")
  }

  # cvforecast
  N <- ifelse(forward, n + h, n + h - 1L)
  nlast <- ifelse(forward, n, n - 1L)
  nfirst <- ifelse(is.null(window), initial, max(window, initial))
  indx <- seq(nfirst, nlast, by = roll.length)
  fit_times <- length(indx)

  pf <- err <- `colnames<-` (ts(matrix(NA_real_, nrow = N, ncol = h),
                                start = start(y), frequency = frequency(y)),
                             paste0("h=", 1:h))
  lower <- upper <- `names<-` (rep(list(pf), length(level)),
                               paste0(level, "%"))

  modelFit = vector(mode = "list", length = indx[length(indx)])
  for (i in indx) {
    print(paste("This is", i))
    train_start <- ifelse(is.null(window), 1L, i - window + 1L)
    suppressWarnings(train <- data1[train_start:i, ] |>
                       as_tibble() |>
                       arrange({{index_data}}) |>
                       select(all_of(index_data), all_of(key_data1), all_of(yvar), all_of(predictor.vars)) |>
                       as_tsibble(index = index_data, key = key_data1))
    suppressWarnings(test <- data1[(i+1):(i+h), ] |>
                       as_tibble() |>
                       arrange({{index_data}}) |>
                       select(all_of(index_data), all_of(key_data1), all_of(predictor.vars)) |>
                       as_tsibble(index = index_data, key = key_data1))
    
    if(any(is.na(test))){
      print(paste0("Skipping a window due to missing values in test set!"))
      next
    }
    
    if(recursive == TRUE){
      recursive_colRange <- suppressWarnings(which(colnames(test) %in% recursive_colNames))
    }else{
      recursive_colRange <- NULL
    }

    # Model fitting in each expanding/rolling window
    key_unique <- unique(as.character(sort(dplyr::pull((train[, {{ key_data1 }}])[, 1]))))
    key_num <- seq_along(key_unique)
    ref <- data.frame(key_unique, key_num)
    train <- train |>
      dplyr::mutate(
        num_key = as.numeric(factor(as.character({{ key_data1 }}), levels = key_unique))
      )

    if ("smimodel" %in% class(object)){
      indexStr <- vector(mode = "list", length = NROW(ref))
      model_list <- vector(mode = "list", length = NROW(ref))
    }else{
      model_list <- vector(mode = "list", length = NROW(ref))
    }

    for (b in seq_along(ref$key_num)){
      # Data filtered by the relevant key
      df_cat <- train |>
        dplyr::filter((abs(num_key - ref$key_num[b]) <= neighbour) |
                        (abs(num_key - ref$key_num[b] + NROW(ref)) <= neighbour) |
                        (abs(num_key - ref$key_num[b] - NROW(ref)) <= neighbour))
      df_cat <- df_cat |>
        drop_na()
      if ("smimodel" %in% class(object)){
        # Index structure for the relevant key
        for(d in 1:ncol(object$fit[[b]]$best$alpha)){
          nonzero <- which(object$fit[[b]]$best$alpha[ , d] != 0)
          indexStr[[b]]$index.vars <- c(indexStr[[b]]$index.vars, names(nonzero))
          indexStr[[b]]$index.ind <- c(indexStr[[b]]$index.ind, rep(d, length(nonzero)))
        }
        # Formula
        ind_pos <- split(seq_along(indexStr[[b]]$index.ind), indexStr[[b]]$index.ind)
        temp <- vector(mode = "list", length = length(ind_pos))
        for(a in 1:length(ind_pos)){
          var_list <- indexStr[[b]]$index.vars[ind_pos[[a]]]
          temp[[a]] <- ifelse(length(var_list) == 1, 1, 2)
        }
        temp <- unlist(temp)
        if(all(temp == 1)){
          pre.formula <- lapply(indexStr[[b]]$index.vars, function(var) paste0("s(", var, ")")) |>
            paste(collapse = "+")
          pre.formula <- paste0(yvar, " ~", pre.formula)
        }else{
          var_list <- indexStr[[b]]$index.vars[ind_pos[[1]]]
          if(length(var_list) > 1){
            pre.formula <- lapply(var_list, function(var) paste0(var)) |>
              paste(collapse = ",") |>
              paste0(")")
            pre.formula <- paste0(yvar, " ~ g(", pre.formula)
          }else{
            pre.formula <- lapply(var_list, function(var) paste0(var)) |>
              paste0(")")
            pre.formula <- paste0(yvar, " ~ s(", pre.formula)
          }
          if(length(ind_pos) > 1){
            for(g in 2:length(ind_pos)){
              var_list <- indexStr[[b]]$index.vars[ind_pos[[g]]]
              if(length(var_list) > 1){
                add.formula <- lapply(var_list, function(var) paste0(var)) |>
                  paste(collapse = ",") |>
                  paste0(")")
                pre.formula <- paste0(pre.formula, " + g(", add.formula)
              }else{
                add.formula <- lapply(var_list, function(var) paste0(var)) |>
                  paste0(")")
                pre.formula <- paste0(pre.formula, " +s(", add.formula)
              }
            }
          }
        }
        if (!is.null(object$fit[[1]]$best$vars_s)){
          s.formula <- lapply(object$fit[[1]]$best$vars_s, function(var) paste0("s(", var, ")")) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", s.formula)
        }
        if (!is.null(object$fit[[1]]$best$vars_linear)){
          linear.formula <- lapply(object$fit[[1]]$best$vars_linear, function(var) paste0(var)) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", linear.formula)
        }
        if(all(temp == 1)){
          # Model fitting
          model_list[[b]] <- mgcv::gam(formula = as.formula(pre.formula),
                                       method = "REML",
                                       data = df_cat)
          add <- df_cat |>
            #drop_na() |>
            select({{ index_data }}, {{ key_data1 }})
          model_list[[b]]$model <- bind_cols(add, model_list[[b]]$model)
          model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                              index = index_data,
                                              key = all_of(key_data1))
        }else{
          # Model fitting
          model_list[[b]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                          data = df_cat)
          modelFrame <- model.frame(formula = as.formula(pre.formula),
                                    data = df_cat)
          add <- df_cat |>
            #drop_na() |>
            select({{ index_data }}, {{ key_data1 }})
          model_list[[b]]$model <- bind_cols(add, modelFrame)
          model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                              index = index_data,
                                              key = all_of(key_data1))
        }
      }else if("backward" %in% class(object)){
        # Model fitting
        model_list[[b]] <- mgcv::gam(formula = object$fit[[b]]$formula,
                                     family = object$fit[[b]]$family$family,
                                     method = "REML",
                                     data = df_cat)
        add <- df_cat |>
          #drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, model_list[[b]]$model)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }else if("gaimFit" %in% class(object)){
        pre.formula <- object$fit[[b]]$terms
        attributes(pre.formula) <- NULL
        model_list[[b]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                        data = df_cat)
        modelFrame <- model.frame(formula = as.formula(pre.formula),
                                  data = df_cat)
        add <- df_cat |>
          #drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, modelFrame)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }else if("pprFit" %in% class(object)){
        pre.formula <- object$fit[[b]]$terms
        attributes(pre.formula) <- NULL
        model_list[[b]] <- stats::ppr(formula = as.formula(pre.formula),
                                      data = df_cat, nterms = object$fit[[b]]$mu)
        add <- df_cat |>
          #drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, model_list[[b]]$model)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }
    }
    # Structuring the output
    data_list <- list(key_unique, model_list)
    modelFit[[i]] <- tibble::as_tibble(
      x = data_list, .rows = length(data_list[[1]]),
      .name_repair = ~ make.names(names = c("key", "fit"))
    )
    if ("smimodel" %in% class(object)){
      class(modelFit[[i]]) <- c("gaimFit", "tbl_df", "tbl", "data.frame")
    }else if("backward" %in% class(object)){
      class(modelFit[[i]]) <- c("backward", "tbl_df", "tbl", "data.frame")
    }else if("gaimFit" %in% class(object)){
      class(modelFit[[i]]) <- c("gaimFit", "tbl_df", "tbl", "data.frame")
    }else if("pprFit" %in% class(object)){
      class(modelFit[[i]]) <- c("pprFit", "tbl_df", "tbl", "data.frame")
    }

    # Obtain predictions on test set
    preds <- predict(object = modelFit[[i]],
                     newdata = test,
                     exclude.trunc = exclude.trunc,
                     recursive = recursive,
                     recursive_colRange = recursive_colRange)
    # Convert to a tibble
    preds <- preds |> 
      tibble::as_tibble() |>
      dplyr::arrange({{index_data}})
    # Store predictions and non-conformity scores
    pf[i, ] <- as.numeric(preds$.predict)
    err[i, ] <- y[i + 1:h] - as.numeric(preds$.predict)
  }

  ### Bootstrapping non-conformity scores with rolling calibration sets
  errors_temp <- na.omit(err)
  errors <- ts(errors_temp[1:NROW(errors_temp) - 1, ], start = start(errors_temp), frequency = frequency(errors_temp))
  if (ncal > NROW(errors))
    stop("`ncal` is larger than the number of rows in the matrix of non-conformity scores.")

  indx_cal <- seq(ncal, NROW(errors), by = 1)
  for(j in indx_cal){
    print(paste("This is", j))
    # Calibration set
    errors_subset <- subset(errors,
                            start = j - ncal + 1L,
                            end = j)
    ## Generate the matrix of bootstrapped series
    # Bootstrapped row indices
    sample_row_indx <- sample(seq(NROW(errors_subset)), num.futures, replace = TRUE)
    # Matrix of bootstrapped series (h*num.futures)
    bootstraps <- t(errors_subset[sample_row_indx, ])
    colnames(bootstraps) <- seq(1, num.futures)
    ## Possible futures through bootstrapping non-conformity scores
    if(recursive == FALSE){
      # Predictions for the rolling test set (of length h)
      preds <- pf[end(errors_subset)[1] + 1, ]
      npreds <- length(preds)
      # Generate matrix of possible futures
      possibleFutures <- vector(mode = "list", length = npreds)
      for(m in 1:npreds){
        possibleFutures[[m]] <- preds[m] + bootstraps[m, ]
      }
      possibleFutures_mat <- as.matrix(bind_rows(possibleFutures))
    }else if(recursive == TRUE){
      # Rolling test set (of length h)
      newdata <- data1[(end(errors_subset)[1] + 2):(end(errors_subset)[1] + 1 + h), ] |>
        as_tsibble(index = index_data, key = key_data1)
      # recursive_colRange
      recursive_colRange <- suppressWarnings(which(colnames(newdata) %in% recursive_colNames))
      # Forecasting model estimated on the most recent training window
      forecastModel <- modelFit[[end(errors_subset)[1] + 1]]
      # Objects of class "smimodel" do not occur here. Hence, using
      # possibleFutures_benchmark()
      futures <- possibleFutures_benchmark(object = forecastModel,
                                           newdata = newdata,
                                           bootstraps = bootstraps,
                                           exclude.trunc = exclude.trunc,
                                           recursive_colRange = recursive_colRange)
      names(futures$future_cols) <- 1:length(futures$future_cols)
      possibleFutures_part2 <- as.matrix(bind_cols(futures$future_cols))
      possibleFutures_part1 <- matrix(futures$firstFuture, nrow = 1,
                                      ncol = length(futures$firstFuture))
      possibleFutures_mat <- rbind(possibleFutures_part1, possibleFutures_part2)
    }

    # Prediction interval bounds
    lower_q <- (1 - (level/100))/2
    upper_q <- lower_q + (level/100)
    for(p in 1:NROW(possibleFutures_mat)){
      intervalsHilo <- quantile(possibleFutures_mat[p, ], probs = c(lower_q, upper_q), na.rm = na.rm)
      nint <- length(level)
      lowerBound <- matrix(NA, ncol = nint, nrow = 1)
      upperBound <- lowerBound
      for (k in 1:nint) {
        lowerBound[1, k] <- intervalsHilo[1:nint][[k]]
        upperBound[1, k] <- intervalsHilo[(nint + 1):length(intervalsHilo)][[k]]
      }
      colnames(lowerBound) <- colnames(upperBound) <- paste(level, "%", sep = "")
      for (l in level) {
        levelname <- paste0(l, "%")
        lower[[levelname]][end(errors_subset)[1] + 1, p] <- lowerBound[ , levelname]
        upper[[levelname]][end(errors_subset)[1] + 1, p] <- upperBound[ , levelname]
      }
    }
  }

  # Prepare return
  out <- list(x = y)
  out$method <- paste("cb_cvforecast")
  out$fit_times <- fit_times
  out$mean <- leadlagMat(pf, 1:h) |> window(start = time(pf)[nfirst + 1L])
  out$error <- leadlagMat(err, 1:h) |> window(start = time(err)[nfirst + 1L], end = time(err)[n])
  out$level <- level
  out$lower <- lapply(lower,
                      function(low) leadlagMat(low, 1:h) |>
                        window(start = time(low)[nfirst + 1L]))
  out$upper <- lapply(upper,
                      function(up) leadlagMat(up, 1:h) |>
                        window(start = time(up)[nfirst + 1L]))

  return(structure(out, class = "cb_cvforecast"))
}


#' Single season block bootstrap prediction intervals through time series
#' cross-validation forecasting
#'
#' Compute prediction intervals by applying the single season block bootstrap
#' method to subsets of time series data using a rolling forecast origin.
#'
#' @param object Fitted model object of class \code{smimodel}, \code{backward},
#'   \code{gaimFit} or \code{pprFit}.
#' @param data Data set. Must be a data set of class \code{tsibble}.(Make sure
#'   there are no additional date or time related variables except for the
#'   \code{index} of the \code{tsibble}). If multiple models are fitted, the
#'   grouping variable should be the \code{key} of the \code{tsibble}. If a
#'   \code{key} is not specified, a dummy key with only one level will be
#'   created.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an \code{integer}. If \code{neighbour =
#'   x}, \code{x} number of keys before the key of interest and \code{x} number
#'   of keys after the key of interest are grouped together for model fitting.
#'   The default is \code{neighbour = 0} (i.e. no neighbours are considered for
#'   model fitting).
#' @param predictor.vars A character vector of names of the predictor variables.
#' @param h Forecast horizon.
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = \code{NULLseason.period * m})
#' @param num.futures Number of possible future sample paths to be generated.
#' @param level Confidence level for prediction intervals.
#' @param forward If \code{TRUE}, the final forecast origin for forecasting is
#'   \eqn{y_T}. Otherwise, the final forecast origin is \eqn{y_{T-1}}.
#' @param initial Initial period of the time series where no cross-validation
#'   forecasting is performed.
#' @param window Length of the rolling window. If \code{NULL}, a rolling window
#'   will not be used.
#' @param roll.length Number of observations by which each rolling/expanding
#'   window should be rolled forward.
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable is treated non-linearly in the
#'   estimated model, will be truncated to be in the in-sample range before
#'   obtaining predictions. If any variables are listed here will be excluded
#'   from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colNames If \code{recursive = TRUE}, a character vector
#'   giving the names of the columns in test data to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{data}, with no break in the lagged variable sequence even if
#'   some of the intermediate lags are not used as predictors.
#' @param na.rm logical; if \code{TRUE} (default), any \code{NA} and
#'   \code{NaN}'s are removed from the sample before the quantiles are computed.
#' @param ... Other arguments not currently used.
#' @return An object of class \code{bb_cvforecast}, which is a list that
#'   contains following elements: \item{x}{The original time series.}
#'   \item{method}{A character string "bb_cvforecast".} \item{fit_times}{The
#'   number of times the model is fitted in cross-validation.}
#' \item{mean}{Point forecasts as a multivariate time series, where the
#' \eqn{h^{th}} column holds the point forecasts for forecast horizon \eqn{h}.
#' The time index corresponds to the period for which the forecast is produced.}
#' \item{res}{The matrix of in-sample residuals produced in cross-validation.
#' The number of rows corresponds to \code{fit_times}, and the row names
#' corresponds the time index of the forecast origin of the corresponding
#' cross-validation iteration.} \item{level}{The confidence values
#'  associated with the prediction intervals.} \item{lower}{A list containing
#'  lower bounds for prediction intervals for each level. Each element within
#'  the list will be a multivariate time series with the same dimensional
#'  characteristics as \code{mean}.} \item{upper}{A list containing upper bounds
#'   for prediction intervals for each level. Each element within the list will
#'   be a multivariate time series with the same dimensional characteristics as
#'   \code{mean}.}
#'
#' @seealso \code{\link{cb_cvforecast}}
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1105
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 +
#'     (0.35*x_lag_002 + 0.7*x_lag_005)^2 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1100, ]
#'
#' # Model fitting
#' smimodel_ppr <- model_smimodel(data = sim_train,
#'                                yvar = "y",
#'                                index.vars = index.vars,
#'                                initialise = "ppr")
#'
#' # Block bootstrap prediction intervals (3-steps-ahead interval forecasts)
#' set.seed(12345)
#' smimodel_ppr_bb <- bb_cvforecast(object = smimodel_ppr,
#'                                  data = sim_data,
#'                                  yvar = "y",
#'                                  predictor.vars = index.vars,
#'                                  h = 3,
#'                                  num.futures = 50,
#'                                  window = 1000)
#' }
#'
#' @export
bb_cvforecast <- function(object, data,
                          yvar, neighbour = 0, predictor.vars,
                          h = 1, season.period = 1, m = 1,
                          num.futures = 1000, level = c(80, 95), forward = TRUE,
                          initial = 1, window = NULL, roll.length = 1,
                          exclude.trunc = NULL,
                          recursive = FALSE, recursive_colNames = NULL, 
                          na.rm = TRUE, ...) {
  # Check input data
  if (!is_tsibble(data)) stop("data is not a tsibble.")

  index_data <- index(data)
  if (length(key(data)) == 0) {
    data <- data |>
      dplyr::mutate(dummy_key = rep(1, NROW(data))) |>
      tsibble::as_tsibble(index = index_data, key = dummy_key)
  }
  key_data1 <- key(data)[[1]]

  data1 <- data |>
    as_tibble() |>
    arrange({{index_data}})
  y <- as.ts(data1[ , yvar][[1]])
  n <- length(y) - forward * h

  # Check confidence levels
  if (min(level) > 0 && max(level) < 1) {
    level <- 100 * level
  } else if (min(level) < 0 || max(level) > 99.99) {
    stop("confidence limit out of range")
  }
  level <- sort(level)

  # Check other inputs
  if (h <= 0)
    stop("forecast horizon out of bounds")
  if (initial < 1 | initial > n)
    stop("initial period out of bounds")
  if (initial == n && !forward)
    stop("initial period out of bounds")
  if (!is.null(window)) {
    if (window < 1 | window > n)
      stop("window out of bounds")
    if (window == n && !forward)
      stop("window out of bounds")
  }

  # cvforecast
  N <- ifelse(forward, n + h, n + h - 1L)
  nlast <- ifelse(forward, n, n - 1L)
  nfirst <- ifelse(is.null(window), initial, max(window, initial))
  indx <- seq(nfirst, nlast, by = roll.length)
  fit_times <- length(indx)

  pf <- `colnames<-` (ts(matrix(NA_real_, nrow = N, ncol = h),
                         start = start(y), frequency = frequency(y)),
                      paste0("h=", 1:h))
  if(is.null(window)){
    res <- ts(matrix(NA_real_, nrow = N, ncol = max(indx)),
              start = start(y), frequency = frequency(y))
  }else{
    res <- ts(matrix(NA_real_, nrow = N, ncol = window),
              start = start(y), frequency = frequency(y))
  }

  lower <- upper <- `names<-` (rep(list(pf), length(level)),
                               paste0(level, "%"))

  modelFit = vector(mode = "list", length = length(indx))
  pFutures = vector(mode = "list", length = length(indx))
  for (i in indx) {
    print(paste("This is", i))
    train_start <- ifelse(is.null(window), 1L, i - window + 1L)
    suppressWarnings(train <- data1[train_start:i, ] |>
                       as_tibble() |>
                       arrange({{index_data}}) |>
                       select(all_of(index_data), all_of(key_data1), all_of(yvar), all_of(predictor.vars)) |>
                       as_tsibble(index = index_data, key = key_data1))
    suppressWarnings(test <- data1[(i+1):(i+h), ] |>
                       as_tibble() |>
                       arrange({{index_data}}) |>
                       select(all_of(index_data), all_of(key_data1), all_of(predictor.vars)) |>
                       as_tsibble(index = index_data, key = key_data1))
    
    if(any(is.na(test))){
      print(paste0("Skipping a window due to missing values in test set!"))
      next
    }
    
    if(recursive == TRUE){
      recursive_colRange <- suppressWarnings(which(colnames(test) %in% recursive_colNames))
    }else{
      recursive_colRange <- NULL
    }

    # Model fitting in each expanding/rolling window
    key_unique <- unique(as.character(sort(dplyr::pull((train[, {{ key_data1 }}])[, 1]))))
    key_num <- seq_along(key_unique)
    ref <- data.frame(key_unique, key_num)
    train <- train |>
      dplyr::mutate(
        num_key = as.numeric(factor(as.character({{ key_data1 }}), levels = key_unique))
      )

    if ("smimodel" %in% class(object)){
      indexStr <- vector(mode = "list", length = NROW(ref))
      model_list <- vector(mode = "list", length = NROW(ref))
    }else{
      model_list <- vector(mode = "list", length = NROW(ref))
    }

    for (b in seq_along(ref$key_num)){
      # Data filtered by the relevant key
      df_cat <- train |>
        dplyr::filter((abs(num_key - ref$key_num[b]) <= neighbour) |
                        (abs(num_key - ref$key_num[b] + NROW(ref)) <= neighbour) |
                        (abs(num_key - ref$key_num[b] - NROW(ref)) <= neighbour))
      df_cat <- df_cat |>
        drop_na()
      if ("smimodel" %in% class(object)){
        # Index structure for the relevant key
        for(d in 1:ncol(object$fit[[b]]$best$alpha)){
          nonzero <- which(object$fit[[b]]$best$alpha[ , d] != 0)
          indexStr[[b]]$index.vars <- c(indexStr[[b]]$index.vars, names(nonzero))
          indexStr[[b]]$index.ind <- c(indexStr[[b]]$index.ind, rep(d, length(nonzero)))
        }
        # Formula
        ind_pos <- split(seq_along(indexStr[[b]]$index.ind), indexStr[[b]]$index.ind)
        temp <- vector(mode = "list", length = length(ind_pos))
        for(a in 1:length(ind_pos)){
          var_list <- indexStr[[b]]$index.vars[ind_pos[[a]]]
          temp[[a]] <- ifelse(length(var_list) == 1, 1, 2)
        }
        temp <- unlist(temp)
        if(all(temp == 1)){
          pre.formula <- lapply(indexStr[[b]]$index.vars, function(var) paste0("s(", var, ")")) |>
            paste(collapse = "+")
          pre.formula <- paste0(yvar, " ~", pre.formula)
        }else{
          var_list <- indexStr[[b]]$index.vars[ind_pos[[1]]]
          if(length(var_list) > 1){
            pre.formula <- lapply(var_list, function(var) paste0(var)) |>
              paste(collapse = ",") |>
              paste0(")")
            pre.formula <- paste0(yvar, " ~ g(", pre.formula)
          }else{
            pre.formula <- lapply(var_list, function(var) paste0(var)) |>
              paste0(")")
            pre.formula <- paste0(yvar, " ~ s(", pre.formula)
          }
          if(length(ind_pos) > 1){
            for(g in 2:length(ind_pos)){
              var_list <- indexStr[[b]]$index.vars[ind_pos[[g]]]
              if(length(var_list) > 1){
                add.formula <- lapply(var_list, function(var) paste0(var)) |>
                  paste(collapse = ",") |>
                  paste0(")")
                pre.formula <- paste0(pre.formula, " + g(", add.formula)
              }else{
                add.formula <- lapply(var_list, function(var) paste0(var)) |>
                  paste0(")")
                pre.formula <- paste0(pre.formula, " +s(", add.formula)
              }
            }
          }
        }
        if (!is.null(object$fit[[1]]$best$vars_s)){
          s.formula <- lapply(object$fit[[1]]$best$vars_s, function(var) paste0("s(", var, ")")) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", s.formula)
        }
        if (!is.null(object$fit[[1]]$best$vars_linear)){
          linear.formula <- lapply(object$fit[[1]]$best$vars_linear, function(var) paste0(var)) |>
            paste(collapse = "+")
          pre.formula <- paste(pre.formula, "+", linear.formula)
        }
        if(all(temp == 1)){
          # Model fitting
          model_list[[b]] <- mgcv::gam(formula = as.formula(pre.formula),
                                       method = "REML",
                                       data = df_cat)
          add <- df_cat |>
            #drop_na() |>
            select({{ index_data }}, {{ key_data1 }})
          model_list[[b]]$model <- bind_cols(add, model_list[[b]]$model)
          model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                              index = index_data,
                                              key = all_of(key_data1))
        }else{
          # Model fitting
          model_list[[b]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                          data = df_cat)
          modelFrame <- model.frame(formula = as.formula(pre.formula),
                                    data = df_cat)
          add <- df_cat |>
            #drop_na() |>
            select({{ index_data }}, {{ key_data1 }})
          model_list[[b]]$model <- bind_cols(add, modelFrame)
          model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                              index = index_data,
                                              key = all_of(key_data1))
        }
      }else if("backward" %in% class(object)){
        # Model fitting
        model_list[[b]] <- mgcv::gam(formula = object$fit[[b]]$formula,
                                     family = object$fit[[b]]$family$family,
                                     method = "REML",
                                     data = df_cat)
        add <- df_cat |>
          #drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, model_list[[b]]$model)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }else if("gaimFit" %in% class(object)){
        pre.formula <- object$fit[[b]]$terms
        attributes(pre.formula) <- NULL
        model_list[[b]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                        data = df_cat)
        modelFrame <- model.frame(formula = as.formula(pre.formula),
                                  data = df_cat)
        add <- df_cat |>
          #drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, modelFrame)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }else if("pprFit" %in% class(object)){
        pre.formula <- object$fit[[b]]$terms
        attributes(pre.formula) <- NULL
        model_list[[b]] <- stats::ppr(formula = as.formula(pre.formula),
                                      data = df_cat, nterms = object$fit[[b]]$mu)
        add <- df_cat |>
          #drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, model_list[[b]]$model)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }
    }
    # Structuring the output
    data_list <- list(key_unique, model_list)
    modelFit[[i]] <- tibble::as_tibble(
      x = data_list, .rows = length(data_list[[1]]),
      .name_repair = ~ make.names(names = c("key", "fit"))
    )
    if ("smimodel" %in% class(object)){
      class(modelFit[[i]]) <- c("gaimFit", "tbl_df", "tbl", "data.frame")
    }else if("backward" %in% class(object)){
      class(modelFit[[i]]) <- c("backward", "tbl_df", "tbl", "data.frame")
    }else if("gaimFit" %in% class(object)){
      class(modelFit[[i]]) <- c("gaimFit", "tbl_df", "tbl", "data.frame")
    }else if("pprFit" %in% class(object)){
      class(modelFit[[i]]) <- c("pprFit", "tbl_df", "tbl", "data.frame")
    }

    # Obtain predictions on test set
    preds <- predict(object = modelFit[[i]],
                     newdata = test,
                     exclude.trunc = exclude.trunc,
                     recursive = recursive,
                     recursive_colRange = recursive_colRange)
    # Store predictions
    pf[i, ] <- as.numeric(preds$.predict)
    # Obtain in-sample residuals
    resids <- augment(x = modelFit[[i]])
    # Store residuls
    res[i, 1:length(resids$.resid)] <- as.numeric(resids$.resid)
    # Possible futures through block bootstrapping on in-sample residuals
    pFutures[[i]] <- blockBootstrap(object = modelFit[[i]], newdata = test,
                                    resids = na.omit(res[i, ]), preds = pf[i, ],
                                    season.period = season.period, m = m,
                                    num.futures = num.futures,
                                    exclude.trunc = exclude.trunc,
                                    recursive = recursive,
                                    recursive_colRange = recursive_colRange)
    # Prediction interval bounds
    lower_q <- (1 - (level/100))/2
    upper_q <- lower_q + (level/100)
    for(j in 1:NROW(pFutures[[i]])){
      intervalsHilo <- quantile(pFutures[[i]][j, ], probs = c(lower_q, upper_q), na.rm = na.rm)
      nint <- length(level)
      lowerBound <- matrix(NA, ncol = nint, nrow = 1)
      upperBound <- lowerBound
      for (k in 1:nint) {
        lowerBound[1, k] <- intervalsHilo[1:nint][[k]]
        upperBound[1, k] <- intervalsHilo[(nint + 1):length(intervalsHilo)][[k]]
      }
      colnames(lowerBound) <- colnames(upperBound) <- paste(level, "%", sep = "")
      for (l in level) {
        levelname <- paste0(l, "%")
        lower[[levelname]][i, j] <- lowerBound[ , levelname]
        upper[[levelname]][i, j] <- upperBound[ , levelname]
      }
    }
  }

  # Prepare return
  out <- list(x = y)
  out$method <- paste("bb_cvforecast")
  out$fit_times <- fit_times
  out$mean <- leadlagMat(pf, 1:h) |> window(start = time(pf)[nfirst + 1L])
  out$res <- res[rowSums(is.na(res)) != ncol(res), ]
  row.names(out$res) <- seq(nfirst, nlast, by = 1)
  out$level <- level
  out$lower <- lapply(lower,
                      function(low) leadlagMat(low, 1:h) |>
                        window(start = time(low)[nfirst + 1L]))
  out$upper <- lapply(upper,
                      function(up) leadlagMat(up, 1:h) |>
                        window(start = time(up)[nfirst + 1L]))

  return(structure(out, class = "bb_cvforecast"))
}
utils::globalVariables(c("indexLag", "indexDiff", "row_idx", "grp"))


#' Prepare a data set for recursive forecasting
#'
#' Prepare a test data for recursive forecasting by appropriately removing
#' exisiting (actual) values from a specified range of columns (lagged response
#' columns) of the data set. Handles seasonal data with gaps.
#'
#' @param newdata Data set to be prepared. Should be a \code{tsibble}.
#' @param recursive_colRange The range of column numbers (lagged response
#'   columns) in \code{newdata} from which existing values should be removed.
#'   Make sure such columns are positioned together in increasing lag order
#'   (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag used) in
#'   \code{newdata}, with no break in the lagged variable sequence even if some
#'   of the intermediate lags are not used as predictors.
#' @return A \code{tibble}.
prep_newdata <- function(newdata, recursive_colRange){
    # Index
    index_data <- index(newdata)
    # Convert to a tibble
    newdata <- newdata |>
      tibble::as_tibble() |>
      dplyr::arrange({{index_data}})
  # Constructing a column to store time difference between two observations
  newdata <- newdata |>
    mutate(indexLag = lag({{index_data}}),
           indexDiff = {{index_data}} - indexLag)
  # Identify times series with seasonal gaps
  # Regular time difference
  timediff <- as.numeric(names(table(newdata$indexDiff))[which.max(table(newdata$indexDiff))])
  # Identify gap locations
  idx <- which(newdata$indexDiff > timediff, arr.ind = TRUE)
  if(length(idx) == 0){
    # Adjusting the test set data to remove future response lags
    newdata <- remove_lags(data = newdata, recursive_colRange = recursive_colRange)
  }else{
    # Split test set at gap locations and save splits as a list
    newdata_list_temp <- newdata |>
      mutate(row_idx = row_number(),
             grp = cumsum(row_idx %in% idx)) |>
      group_split(grp)
    # Remove future response lags from each list element appropriately
    newdata_list <- newdata_list_temp |>
      purrr::map(~ remove_lags(data = ., recursive_colRange = recursive_colRange))
    newdata <- bind_rows(newdata_list)
  }
  return(newdata)
}


#' Remove actual values from a data set for recursive forecasting
#'
#' Appropriately removes exisiting (actual) values from the specified column
#' range (lagged response columns) of a given data set (typicall a test set for
#' which recursive forecasting is required).
#'
#' @param data Data set (a \code{tibble}) from which the actual lagged values
#'   should be removed.
#' @param recursive_colRange The range of column numbers in \code{data} from
#'   which lagged values should be removed.
#' @return A \code{tibble}.
remove_lags <- function(data, recursive_colRange){
  if(NROW(data) > 1){
    if(NROW(data) <= length(recursive_colRange)){
      for(a in recursive_colRange[1:(NROW(data) - 1)]){
        data[(a - (recursive_colRange[1] - 2)):NROW(data), a] <- NA
      }
    }else{
      for(a in recursive_colRange){
        data[(a - (recursive_colRange[1] - 2)):NROW(data), a] <- NA
      }
    }
  }else{
    data
  }
  return(data)
}


#' Futures through single season block bootstrapping
#'
#' Gerenates possible future sample paths by applying the single season block
#' bootstrap method.
#'
#' @param object Fitted model object.
#' @param newdata Test data set. Must be a data set of class \code{tsibble}.
#' @param resids In-sample residuals from the fitted model.
#' @param preds Predictions for the test set (i.e. data for the forecast
#'   horizon).
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = \code{season.period * m})
#' @param num.futures Number of possible future sample paths to be generated.
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If \code{recursive = TRUE}, The range of column
#'   numbers in test data to be filled with forecasts.
#' @return A matrix of simulated future sample paths.
blockBootstrap <- function(object, newdata, resids, preds, season.period = 1,
                           m = 1, num.futures = 1000, exclude.trunc = NULL,
                           recursive = FALSE, recursive_colRange = NULL){

  # Generate the matrix of bootstrapped series
  bootstraps <- residBootstrap(x = resids, season.period = season.period, m = m,
                                          num.bootstrap = num.futures)
  # Generate possible futures
  npreds <- length(preds)
  if(recursive == FALSE){
    possibleFutures <- vector(mode = "list", length = npreds)
    for(i in 1:npreds){
      possibleFutures[[i]] <- preds[i] + bootstraps[i, ]
    }
    possibleFutures_mat <- as.matrix(bind_rows(possibleFutures))
  }else if(recursive == TRUE){
    if ("smimodel" %in% class(object)){
      futures <- possibleFutures_smimodel(object = object,
                                          newdata = newdata,
                                          bootstraps = bootstraps,
                                          exclude.trunc = exclude.trunc,
                                          recursive_colRange = recursive_colRange)
    }else{
      futures <- possibleFutures_benchmark(object = object,
                                           newdata = newdata,
                                           bootstraps = bootstraps,
                                           exclude.trunc = exclude.trunc,
                                           recursive_colRange = recursive_colRange)
    }
    names(futures$future_cols) <- 1:length(futures$future_cols)
    possibleFutures_part2 <- as.matrix(bind_cols(futures$future_cols))
    possibleFutures_part1 <- matrix(futures$firstFuture, nrow = 1,
                                    ncol = length(futures$firstFuture))
    possibleFutures_mat <- rbind(possibleFutures_part1, possibleFutures_part2)
  }
  return(possibleFutures_mat)
}


#' Generate multiple single season block bootstrap series
#'
#' Generates multiple replications of single season block bootstrap series.
#'
#' @param x A series of residuals from which bootstrap series to be generated.
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = \code{season.period * m})
#' @param num.bootstrap Number of bootstrap series to be generated.
#' @return A matrix of bootstrapped series.
residBootstrap <- function(x, season.period = 1, m = 1, num.bootstrap = 1000){
  bootstraps <- vector(mode = "list", length = num.bootstrap)
  for(i in 1:num.bootstrap){
    bootstraps[[i]] <- seasonBootstrap(x = x,
                                       season.period = season.period, m = m)
  }
  names(bootstraps) <- seq(1, num.bootstrap)
  bootstraps_mat <- as.matrix(bind_cols(bootstraps))
  return(bootstraps_mat)
}


#' Single season block bootstrap
#'
#' Generates a single replication of single season block bootstrap series.
#'
#' @param x A series of residuals from which bootstrap series to be generated.
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = \code{season.period * m})
#' @return A \code{numeric} vector.
seasonBootstrap <- function(x, season.period = 1, m = 1){
  n <- length(x)
  # Block size
  block.size <- season.period*m
  # Number of (complete) blocks
  nblocks <- trunc(n/season.period/m)
  # Construct a simulated residual series through random re-sampling of blocks
  newx <- as.vector(replicate(nblocks, randomBlock(series = x, block.size = block.size)))
  return(na.omit(newx))
}



#' Randomly sampling a block
#'
#' Samples a block of specified size from a given series starting form a random
#' point in the series.
#'
#' @param series A series from which a block should be sampled.
#' @param block.size Size of the block to be sampled.
#' @return A \code{numeric} vector.
randomBlock <- function(series, block.size){
  start_ind <- sample(seq(length(series) - block.size + 1), 1)
  block <- series[start_ind:(start_ind + block.size -1)]
  return(block)
}


#' Possible future sample paths (multi-step) from \code{smimodel} residuals
#'
#' Generates possible future sample paths (multi-step) using residuals of a
#' fitted \code{smimodel} through recursive forecasting.
#'
#' @param object A \code{smimodel} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param bootstraps Generated matrix of bootstrapped residual series.
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive_colRange The range of column numbers in \code{newdata} to be
#'   filled with forecasts.
#' @return A list containing the following components: \item{firstFuture}{A
#'   \code{numeric} vector of 1-step-ahead simulated futures.}
#'   \item{future_cols}{A list of multi-steps-ahead simulated futures, where
#'   each list element corresponds to each 1-step-ahead simulated future in
#'   \code{firstFuture}.} 
possibleFutures_smimodel <- function(object, newdata, bootstraps,
                                     exclude.trunc = NULL, recursive_colRange){
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key11 <- key(newdata)[[1]]

  # Names of the columns to be filled with forecasts
  recursive_colNames <- colnames(newdata)[recursive_colRange]

  # # Predict function to be used
  # predict_fn <- mgcv::predict.gam

  # Prepare newdata for recursive forecasting
  newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)

  # Recursive possible futures
  # First, one-step-ahead possible futures
  data_temp = newdata[1, ]
  key22 = data_temp[ , {{ key11 }}][[1]]
  key22_pos = which(object$key == key22)
  list_index <- object$fit[[key22_pos]]$best$alpha
  numInd <- NCOL(list_index)
  alpha <- vector(mode = "list", length = numInd)
  for(a in 1:numInd){
    alpha[[a]] <- list_index[ , a]
  }
  alpha <- unlist(alpha)
  names(alpha) <- NULL
  if(all(alpha == 0)){
    data_list1 <- data_temp
  }else{
    X_test <- as.matrix(newdata[1, object$fit[[key22_pos]]$best$vars_index])
    # Calculating indices
    ind <- vector(mode = "list", length = numInd)
    for(b in 1:numInd){
      ind[[b]] <- as.numeric(X_test %*% as.matrix(list_index[ , b], ncol = 1))
    }
    names(ind) <- colnames(list_index)
    dat <- tibble::as_tibble(ind)
    data_list1 <- dplyr::bind_cols(data_temp, dat)
  }
  #pred1 <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list1, type = "response")
  pred1 <- predict(object$fit[[key22_pos]]$best, data_list1, 
                   exclude.trunc = exclude.trunc)$.predict
  possibleFutures1 <- as.numeric(pred1) + bootstraps[1, ]

  # Should fill the missing response lags in newdata using each of the
  # possible future values separately
  newdata_list <- vector(mode = "list", length = length(possibleFutures1))
  future_cols <- vector(mode = "list", length = length(possibleFutures1))
  # Separate the columns in recursive_colRange
  # This is done for extracting only the numerical columns (where values should
  # be replaced with possible futures) so that it can converted to a matrix
  # without any unwanted coercion imposed by R (as matrices can contain only one
  # data type.)
  fill_data_temp <- newdata[ , recursive_colRange]
  # Separate the columns not in recursive_colRange
  remain_data_temp <- newdata[ , -recursive_colRange]
  # Column numbers in fill_data_temp
  colRange <- 1:NCOL(fill_data_temp)
  # Cell positions in fill_data_temp to be filled with 1-step ahead possible futures
  x_seq <- seq((1+1), (1+((max(colRange) - min(colRange)) + 1)))
  x_seq <- x_seq[x_seq <= NROW(fill_data_temp)]
  x_seq <- na.omit(x_seq[1:length(colRange)])
  y_seq <- colRange[1:length(x_seq)]
  ref <- data.frame(x_seq, y_seq)
  mat_indx <- apply(ref, 1, function(val)(NROW(fill_data_temp)*(val[2] - 1)) + val[1])
  # Filling
  for(s in 1:length(possibleFutures1)){
    # Convert to a matrix
    fill_data_temp <- as.matrix(fill_data_temp)
    # Identify NA cells to be filled-in
    mat_indx <- mat_indx[which(is.na(fill_data_temp[mat_indx]))]
    # Fill-in
    fill_data_temp[mat_indx] <- possibleFutures1[s]
    # Convert back to a tibble
    fill_data_temp <- as_tibble(fill_data_temp)
    # Combine the columns separated back into a single tibble
    newdata_temp <- bind_cols(remain_data_temp, fill_data_temp)
    # Store in list
    newdata_list[[s]] <- as_tibble(newdata_temp)
    # From 2-steps ahead
    temp_Futures <- vector(mode = "list", length = length(2:NROW(newdata_list[[s]])))
    for(t in 2:(NROW(newdata_list[[s]]) - 1)){
      data_temp = newdata_list[[s]][t, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      list_index <- object$fit[[key22_pos]]$best$alpha
      numInd <- NCOL(list_index)
      alpha <- vector(mode = "list", length = numInd)
      for(d in 1:numInd){
        alpha[[d]] <- list_index[ , d]
      }
      alpha <- unlist(alpha)
      names(alpha) <- NULL
      if(all(alpha == 0)){
        data_list1 <- data_temp
      }else{
        X_test <- as.matrix(newdata_list[[s]][t, object$fit[[key22_pos]]$best$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = numInd)
        for(h in 1:numInd){
          ind[[h]] <- as.numeric(X_test %*% as.matrix(list_index[ , h], ncol = 1))
        }
        names(ind) <- colnames(list_index)
        dat <- tibble::as_tibble(ind)
        data_list1 <- dplyr::bind_cols(data_temp, dat)
      }
      #pred1 <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list1, type = "response")
      pred1 <- predict(object$fit[[key22_pos]]$best, data_list1,
                       exclude.trunc = exclude.trunc)$.predict
      temp_Futures[[t-1]] <- as.numeric(pred1) + bootstraps[t, s]
      recursive_colRange_new <- which(colnames(newdata_list[[s]]) %in% recursive_colNames)
      fill_data_temp <- newdata_list[[s]][ , recursive_colRange_new]
      remain_data_temp <- newdata_list[[s]][ , -recursive_colRange_new]
      colRange <- 1:NCOL(fill_data_temp)
      x_seq <- seq((t+1), (t+((max(colRange) - min(colRange)) + 1)))
      x_seq <- x_seq[x_seq <= NROW(fill_data_temp)]
      x_seq <- na.omit(x_seq[1:length(colRange)])
      y_seq <- colRange[1:length(x_seq)]
      ref <- data.frame(x_seq, y_seq)
      mat_indx <- apply(ref, 1, function(val)(NROW(fill_data_temp)*(val[2] - 1)) + val[1])
      fill_data_temp <- as.matrix(fill_data_temp)
      mat_indx <- mat_indx[which(is.na(fill_data_temp[mat_indx]))]
      fill_data_temp[mat_indx] <- temp_Futures[[t-1]]
      fill_data_temp <- as_tibble(fill_data_temp)
      newdata_temp <- bind_cols(remain_data_temp, fill_data_temp)
      newdata_list[[s]] <- as_tibble(newdata_temp)
    }
    data_temp = newdata_list[[s]][NROW(newdata_list[[s]]), ]
    key22 = data_temp[ , {{ key11 }}][[1]]
    key22_pos = which(object$key == key22)
    list_index <- object$fit[[key22_pos]]$best$alpha
    numInd <- NCOL(list_index)
    alpha <- vector(mode = "list", length = numInd)
    for(w in 1:numInd){
      alpha[[w]] <- list_index[ , w]
    }
    alpha <- unlist(alpha)
    names(alpha) <- NULL
    if(all(alpha == 0)){
      data_list1 <- data_temp
    }else{
      X_test <- as.matrix(newdata_list[[s]][NROW(newdata_list[[s]]),
                                            object$fit[[key22_pos]]$best$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = numInd)
      for(z in 1:numInd){
        ind[[z]] <- as.numeric(X_test %*% as.matrix(list_index[ , z], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list1 <- dplyr::bind_cols(data_temp, dat)
    }
    # pred1 <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list1,
    #                     type = "response")
    pred1 <- predict(object$fit[[key22_pos]]$best, data_list1,
                     exclude.trunc = exclude.trunc)$.predict
    temp_Futures[[NROW(newdata_list[[s]]) - 1]] <- as.numeric(pred1) + bootstraps[NROW(newdata_list[[s]]), s]
    future_cols[[s]] <- unlist(temp_Futures)
  }
  output <- list("firstFuture" = possibleFutures1, "future_cols" = future_cols)
  return(output)
}


#' Possible future sample paths (multi-step) from residuals of a fitted
#' benchmark model
#'
#' Generates possible future sample paths (multi-step) using residuals of a
#' fitted benchmark model through recursive forecasting.
#'
#' @param object A fitted model object of the class \code{backward},
#'   \code{pprFit}, or \code{gaimFit}.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param bootstraps Generated matrix of bootstrapped residual series.
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string.
#' @param recursive_colRange The range of column numbers in \code{newdata} to be
#'   filled with forecasts.
#' @return A list containing the following components: \item{firstFuture}{A
#'   \code{numeric} vector of 1-step-ahead simulated futures.}
#'   \item{future_cols}{A list of multi-steps-ahead simulated futures, where
#'   each list element corresponds to each 1-step-ahead simulated future in
#'   \code{firstFuture}.}  
possibleFutures_benchmark <- function(object, newdata, bootstraps,
                                      exclude.trunc = NULL,
                                      recursive_colRange){
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key11 <- key(newdata)[[1]]

  # Names of the columns to be filled with forecasts
  recursive_colNames <- colnames(newdata)[recursive_colRange]

  # Prepare newdata for recursive forecasting
  newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)

  # Recursive possible futures
  # First, one-step-ahead possible futures
  data_temp = newdata[1, ]
  key22 = data_temp[ , {{ key11 }}][[1]]
  key22_pos = which(object$key == key22)
  object_temp <- object$fit[[key22_pos]]
  if("backward" %in% class(object)){
    ## Avoid extrapolation; truncate non-linear predictors to match in-sample
    gam_cols <- colnames(object_temp$model)
    remove_temp <- as.character(attributes(object_temp$pterms)$variables)
    no_trunc_cols <- unique(c(as.character(index_n), as.character(key11),
                              remove_temp[2:length(remove_temp)], 
                              exclude.trunc))
    trunc_indx <- !(gam_cols %in% no_trunc_cols)
    trunc_cols <- gam_cols[trunc_indx]
    if(length(trunc_cols) != 0){
      data_temp <- truncate_vars(object.data = object_temp$model,
                                 data = data_temp,
                                 cols.trunc = trunc_cols)
    }
    pred1 <- predict(object_temp, data_temp, type = "response")
  }else if("pprFit" %in% class(object)){
    col_retain <- !(colnames(data_temp) %in% c("indexLag", "indexDiff", "row_idx", "grp"))
    if(any(is.na(data_temp[ , col_retain]))){
      pred1 <- NA
    }else{
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      trunc_indx <- !(object_temp$xnames %in% exclude.trunc)
      trunc_cols <- object_temp$xnames[trunc_indx]
      if(length(trunc_cols) != 0){
        data_temp <- truncate_vars(object.data = object_temp$model,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred1 <- predict(object_temp, data_temp, type = "response")
    }
  }else if("gaimFit" %in% class(object)){
    pred1 <- predict(object_temp, data_temp, type = "response")
  }
  #pred1 <- predict(object$fit[[key22_pos]], data_temp, type = "response")
  possibleFutures1 <- as.numeric(pred1) + bootstraps[1, ]
  # Should fill the missing response lags in newdata using each of the
  # possible future values separately
  newdata_list <- vector(mode = "list", length = length(possibleFutures1))
  future_cols <- vector(mode = "list", length = length(possibleFutures1))

  # Separate the columns in recursive_colRange
  # This is done for extracting only the numerical columns (where values should
  # be replaced with possible futures) so that it can converted to a matrix
  # without any unwanted coercion imposed by R (as matrices can contain only one
  # data type.)
  fill_data_temp <- newdata[ , recursive_colRange]
  # Separate the columns not in recursive_colRange
  remain_data_temp <- newdata[ , -recursive_colRange]
  # Column numbers in fill_data_temp
  colRange <- 1:NCOL(fill_data_temp)
  # Cell positions in fill_data_temp to be filled with 1-step ahead possible futures
  x_seq <- seq((1+1), (1+((max(colRange) - min(colRange)) + 1)))
  x_seq <- x_seq[x_seq <= NROW(fill_data_temp)]
  x_seq <- na.omit(x_seq[1:length(colRange)])
  y_seq <- colRange[1:length(x_seq)]
  ref <- data.frame(x_seq, y_seq)
  mat_indx <- apply(ref, 1, function(val)(NROW(fill_data_temp)*(val[2] - 1)) + val[1])
  # Filling
  for(s in 1:length(possibleFutures1)){
    # Convert to a matrix
    fill_data_temp <- as.matrix(fill_data_temp)
    # Identify NA cells to be filled-in
    mat_indx <- mat_indx[which(is.na(fill_data_temp[mat_indx]))]
    # Fill-in
    fill_data_temp[mat_indx] <- possibleFutures1[s]
    # Convert back to a tibble
    fill_data_temp <- as_tibble(fill_data_temp)
    # Combine the columns separated back into a single tibble
    newdata_temp <- bind_cols(remain_data_temp, fill_data_temp)
    # Store in list
    newdata_list[[s]] <- as_tibble(newdata_temp)
    # From 2-steps ahead
    temp_Futures <- vector(mode = "list", length = length(2:NROW(newdata_list[[s]])))
    for(t in 2:(NROW(newdata_list[[s]]) - 1)){
      data_temp = newdata_list[[s]][t, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      object_temp <- object$fit[[key22_pos]]
      if("backward" %in% class(object)){
        ## Avoid extrapolation; truncate non-linear predictors to match in-sample
        gam_cols <- colnames(object_temp$model)
        remove_temp <- as.character(attributes(object_temp$pterms)$variables)
        no_trunc_cols <- unique(c(as.character(index_n), as.character(key11),
                                  remove_temp[2:length(remove_temp)], 
                                  exclude.trunc))
        trunc_indx <- !(gam_cols %in% no_trunc_cols)
        trunc_cols <- gam_cols[trunc_indx]
        if(length(trunc_cols) != 0){
          data_temp <- truncate_vars(object.data = object_temp$model,
                                     data = data_temp,
                                     cols.trunc = trunc_cols)
        }
        pred1 <- predict(object_temp, data_temp, type = "response")
      }else if("pprFit" %in% class(object)){
        col_retain <- !(colnames(data_temp) %in% c("indexLag", "indexDiff", "row_idx", "grp"))
        if(any(is.na(data_temp[ , col_retain]))){
          pred1 <- NA
        }else{
          ## Avoid extrapolation; truncate non-linear predictors to match in-sample
          trunc_indx <- !(object_temp$xnames %in% exclude.trunc)
          trunc_cols <- object_temp$xnames[trunc_indx]
          if(length(trunc_cols) != 0){
            data_temp <- truncate_vars(object.data = object_temp$model,
                                       data = data_temp,
                                       cols.trunc = trunc_cols)
          }
          pred1 <- predict(object_temp, data_temp, type = "response")
        }
      }else if("gaimFit" %in% class(object)){
        pred1 <- predict(object_temp, data_temp, type = "response")
      }
      #pred1 <- predict(object$fit[[key22_pos]], data_temp, type = "response")
      temp_Futures[[t-1]] <- pred1 + bootstraps[t, s]
      recursive_colRange_new <- which(colnames(newdata_list[[s]]) %in% recursive_colNames)
      fill_data_temp <- newdata_list[[s]][ , recursive_colRange_new]
      remain_data_temp <- newdata_list[[s]][ , -recursive_colRange_new]
      colRange <- 1:NCOL(fill_data_temp)
      x_seq <- seq((t+1), (t+((max(colRange) - min(colRange)) + 1)))
      x_seq <- x_seq[x_seq <= NROW(fill_data_temp)]
      x_seq <- na.omit(x_seq[1:length(colRange)])
      y_seq <- colRange[1:length(x_seq)]
      ref <- data.frame(x_seq, y_seq)
      mat_indx <- apply(ref, 1, function(val)(NROW(fill_data_temp)*(val[2] - 1)) + val[1])
      fill_data_temp <- as.matrix(fill_data_temp)
      mat_indx <- mat_indx[which(is.na(fill_data_temp[mat_indx]))]
      fill_data_temp[mat_indx] <- temp_Futures[[t-1]]
      fill_data_temp <- as_tibble(fill_data_temp)
      newdata_temp <- bind_cols(remain_data_temp, fill_data_temp)
      newdata_list[[s]] <- as_tibble(newdata_temp)
    }
    data_temp = newdata_list[[s]][NROW(newdata_list[[s]]), ]
    key22 = data_temp[ , {{ key11 }}][[1]]
    key22_pos = which(object$key == key22)
    object_temp <- object$fit[[key22_pos]]
    if("backward" %in% class(object)){
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      gam_cols <- colnames(object_temp$model)
      remove_temp <- as.character(attributes(object_temp$pterms)$variables)
      no_trunc_cols <- unique(c(as.character(index_n), as.character(key11),
                                remove_temp[2:length(remove_temp)], 
                                exclude.trunc))
      trunc_indx <- !(gam_cols %in% no_trunc_cols)
      trunc_cols <- gam_cols[trunc_indx]
      if(length(trunc_cols) != 0){
        data_temp <- truncate_vars(object.data = object_temp$model,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred1 <- predict(object_temp, data_temp, type = "response")
    }else if("pprFit" %in% class(object)){
      col_retain <- !(colnames(data_temp) %in% c("indexLag", "indexDiff", "row_idx", "grp"))
      if(any(is.na(data_temp[ , col_retain]))){
        pred1 <- NA
      }else{
        ## Avoid extrapolation; truncate non-linear predictors to match in-sample
        trunc_indx <- !(object_temp$xnames %in% exclude.trunc)
        trunc_cols <- object_temp$xnames[trunc_indx]
        if(length(trunc_cols) != 0){
          data_temp <- truncate_vars(object.data = object_temp$model,
                                     data = data_temp,
                                     cols.trunc = trunc_cols)
        }
        pred1 <- predict(object_temp, data_temp, type = "response")
      }
    }else if("gaimFit" %in% class(object)){
      pred1 <- predict(object_temp, data_temp, type = "response")
    }
    #pred1 <- predict(object$fit[[key22_pos]], data_temp, type = "response")
    temp_Futures[[NROW(newdata_list[[s]]) - 1]] <- pred1 + bootstraps[NROW(newdata_list[[s]]), s]
    future_cols[[s]] <- unlist(temp_Futures)
  }
  output <- list("firstFuture" = possibleFutures1, "future_cols" = future_cols)
  return(output)
}


#' Calculate interval forecast coverage
#'
#' This is a wrapper for the function \code{conformalForecast::coverage}.
#' Calculates the mean coverage and the ifinn matrix for prediction intervals on
#' validation set. If \code{window} is not \code{NULL}, a matrix of the rolling
#' means of interval forecast coverage is also returned.
#'
#' @param object An object of class \code{bb_cvforecast}.
#' @param level Target confidence level for prediction intervals.
#' @param window If not \code{NULL}, the rolling mean matrix for coverage is
#'   also returned.
#' @param na.rm A \code{logical} indicating whether \code{NA} values should be
#'   stripped before the rolling mean computation proceeds.
#'
#' @return A list of class \code{coverage} with the following components:
#'   \item{mean}{Mean coverage across the validation set.}
#' \item{ifinn}{A indicator matrix as a multivariate time series, where the
#' \eqn{h}th column holds the coverage for forecast horizon \eqn{h}. The time
#' index corresponds to the period for which the forecast is produced.}
#' \item{rollmean}{If \code{window} is not \code{NULL}, a matrix of the rolling
#' means of interval forecast coverage will be returned.}
#'
#' @export
avgCoverage <- function(object, level = 95, window = NULL, na.rm = FALSE) {
  object$LOWER <- object$lower
  object$UPPER <- object$upper
  output <- coverage(object = object, level = level,
                                        window = window, na.rm = na.rm)
  return(output)
}


#' Calculate interval forecast width
#'
#' This is a wrapper for the function \code{conformalForecast::width}.
#' Calculates the mean width of prediction intervals on the validation set. If
#' \code{window} is not \code{NULL}, a matrix of the rolling means of interval
#' width is also returned. If \code{includemedian} is \code{TRUE}, the
#' information of the median interval width will be returned.
#'
#' @param object An object of class  \code{bb_cvforecast}.
#' @param level Target confidence level for prediction intervals.
#' @param includemedian If \code{TRUE}, the median interval width will also be
#'   returned.
#' @param window If not \code{NULL}, the rolling mean (and rolling median if
#'   applicable) matrix for interval width will also be returned.
#' @param na.rm A logical indicating whether \code{NA} values should be stripped
#'   before the rolling mean and rolling median computation proceeds.
#'
#' @return A list of class \code{width} with the following components:
#' \item{width}{Forecast interval width as a multivariate time series, where the
#'  \eqn{h}th column holds the interval width for the forecast horizon \eqn{h}.
#'  The time index corresponds to the period for which the forecast is
#'  produced.} \item{mean}{Mean interval width across the validation set.}
#' \item{rollmean}{If \code{window} is not \code{NULL}, a matrix of the rolling
#' means of interval width will be returned.} \item{median}{Median interval
#' width across the validation set.} \item{rollmedian}{If \code{window} is not
#' \code{NULL}, a matrix of the rolling medians of interval width will be
#' returned.}
#'
#' @export
avgWidth <- function(object, level = 95, includemedian = FALSE, window = NULL,
                     na.rm = FALSE) {
  object$LOWER <- object$lower
  object$UPPER <- object$upper
  output <- width(object = object, level = level,
                     includemedian = includemedian,
                     window = window, na.rm = na.rm)
  return(output)
}


#' Create lags or leads of a matrix
#'
#' This is a wrapper for the function \code{conformalForecast::lagmatrix}.Find a
#' shifted version of a matrix, adjusting the time base backward (lagged) or
#' forward (leading) by a specified number of observations for each column.
#'
#' @param x A matrix or multivariate time series.
#' @param lag A vector of lags (positive values) or leads (negative values) with
#'   a length equal to the number of columns of \code{x}.
#'
#' @return A matrix with the same class and size as \code{x}.
#'
#' @examples
#' x <- matrix(rnorm(20), nrow = 5, ncol = 4)
#'
#' # Create lags of a matrix
#' leadlagMat(x, c(0, 1, 2, 3))
#'
#' # Create leads of a matrix
#' leadlagMat(x, c(0, -1, -2, -3))
#'
#' @export
leadlagMat <- function(x, lag) {
  lmat <- lagmatrix(x = x, lag = lag)
  return(lmat)
}
