#' Single season block bootstrap prediction intervals through cross-validation
#'
#' Compute prediction intervals by applying the single season block bootstrap
#' method to subsets of time series data using a rolling forecast origin.
#'
#' @param object Fitted model object of class `smimodel`, `backward`, `gaimFit`
#'   or `pprFit`.
#' @param data Data set. Must be a data set of class `tsibble`.(Make
#'   sure there are no additional date/time/date-time/yearmonth/POSIXct/POSIXt
#'   variables except for the `index` of the `tsibble`). If multiple models are
#'   fitted, the grouping variable should be the key of the `tsibble`.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param predictor.vars A character vector of names of the predictor variables.
#' @param h Forecast horizon.
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = `season.period` * `m`)
#' @param num.futures Number of possible future sample paths to be generated.
#' @param level Confidence level for prediction intervals.
#' @param forward If \code{TRUE}, the final forecast origin for forecasting is
#'   \eqn{y_T}. Otherwise, the final forecast origin is \eqn{y_{T-1}}.
#' @param initial Initial period of the time series where no cross-validation
#'   forecasting is performed.
#' @param window Length of the rolling window. If \code{NULL}, a rolling window
#'   will not be used.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colNames If `recursive = TRUE`, a character vector giving
#'   the names of the columns in test data to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @export
crossVal_bb <- function(object, data, #newdata, 
                        yvar, neighbour = 0, predictor.vars, 
                        h = 1, season.period = 1, m = 1, 
                        num.futures = 1000, level = c(80, 95), forward = TRUE, 
                        initial = 1, window = NULL, 
                        recursive = FALSE, recursive_colNames = NULL, ...) {
  # Check input data 
  if (!is_tsibble(data)) stop("data is not a tsibble.")
  #if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_data <- index(data)
  key_data <- key(data)
  if (length(key(data)) == 0) {
    data <- data |>
      dplyr::mutate(dummy_key = rep(1, NROW(data))) |>
      tsibble::as_tsibble(index = index_data, key = dummy_key)
    key_data <- key(data)
    # newdata <- newdata |>
    #   dplyr::mutate(dummy_key = rep(1, NROW(newdata))) |>
    #   tsibble::as_tsibble(index = index_data, key = dummy_key)
  }
  key_data1 <- key(data)[[1]]
  data1 <- data |> 
    as_tibble() |> 
    arrange({{index_data}})
  test.length <- forward * h
  y <- as.ts(data1[ , yvar][[1]][1:(NROW(data1) - test.length)])
  #y <- as.ts(data1[ , yvar][[1]])
  n <- length(y)
  
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
  # forwardData <- bind_rows(data, newdata)
  # forwardData <- forwardData |> 
  #   as_tibble() |> 
  #   arrange({{index_data}})
  # xreg <- forwardData[ , c(as.character(index_data), as.character(key_data1), 
  #                          predictor.vars)]
  xreg <- data1[ , c(as.character(index_data), as.character(key_data1), 
                           predictor.vars)]
  xreg <- ts(xreg,
             start = start(y),
             frequency = frequency(y))
  if (nrow(xreg) < n)
    stop("xreg should be at least of the same size as y")
  if (nrow(xreg) < n + forward * h)
    # Pad xreg with NAs
    xreg <- ts(rbind(xreg, matrix(NA, nrow=n+forward*h-nrow(xreg), ncol=ncol(xreg))),
               start = start(y),
               frequency = frequency(y))
  if (nrow(xreg) > n + forward * h) {
    warning(sprintf("only first %s rows of xreg are being used", n + forward * h))
    xreg <- subset(xreg, start = 1L, end = n + forward * h)
  }
  if (is.null(colnames(xreg))) {
    colnames(xreg) <- if (ncol(xreg) == 1) "xreg" else paste0("xreg", 1:ncol(xreg))
  }
  
  N <- ifelse(forward, n + h, n + h - 1L)
  nlast <- ifelse(forward, n, n - 1L)
  nfirst <- ifelse(is.null(window), initial, max(window, initial))
  indx <- seq(nfirst, nlast, by = 1L)
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
  out <- list(x = y)
  
  modelFit = vector(mode = "list", length = length(indx))
  pFutures = vector(mode = "list", length = length(indx))
  for (i in indx) {
    y_subset <- subset(y,
                       start = ifelse(is.null(window), 1L, i - window + 1L),
                       end = i)
    xreg_subset <- subset(xreg,
                          start = ifelse(is.null(window), 1L, i - window + 1L),
                          end = i)
    xreg_future <- subset(xreg,
                          start = i + 1L,
                          end = i + h)
    suppressWarnings(train <- xreg_subset |>
                       as_tibble() |> 
                       arrange({{index_data}}) |> 
                       mutate(!!yvar := as.numeric(y_subset)) |> 
                       select(all_of(index_data), all_of(key_data1), all_of(yvar), all_of(predictor.vars)) |> 
                       as_tsibble(index = index_data, key = key_data1))
    suppressWarnings(test <- xreg_future |>
                       as_tibble() |> 
                       arrange({{index_data}}) |> 
                       select(all_of(index_data), all_of(key_data1), all_of(predictor.vars)) |> 
                       as_tsibble(index = index_data, key = key_data1))
    if(recursive == TRUE){
      recursive_colRange <- suppressWarnings(which(colnames(test) %in% recursive_colNames))
      # Convert to a tibble
      test <- test |> 
        as_tibble() |> 
        arrange({{index_data}})
      ## Adjusting the test set data to remove future response lags
      for(a in recursive_colRange){
        test[(a - (recursive_colRange[1] - 2)):NROW(test), a] <- NA
      }
      # Convert back to a tsibble
      test <- test |>
        as_tsibble(index = index_data, key = key_data1)
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
      if ("smimodel" %in% class(object)){
        # Index structure for the relevant key
        for(d in 1:ncol(object$fit[[b]]$best$alpha)){
          nonzero <- which(object$fit[[b]]$best$alpha[ , d] != 0)
          indexStr[[b]]$index.vars <- c(indexStr[[b]]$index.vars, names(nonzero))
          indexStr[[b]]$index.ind <- c(indexStr[[b]]$index.ind, rep(d, length(nonzero)))
        }
        # Formula
        ind_pos <- split(seq_along(indexStr[[b]]$index.ind), indexStr[[b]]$index.ind)
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
          for(j in 2:length(ind_pos)){
            var_list <- indexStr[[b]]$index.vars[ind_pos[[j]]]
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
        # Model fitting
        model_list[[b]] <- cgaim::cgaim(formula = as.formula(pre.formula),
                                        data = df_cat)
        modelFrame <- model.frame(formula = as.formula(pre.formula), 
                                  data = df_cat)
        add <- df_cat |>
          drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, modelFrame)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }else if("backward" %in% class(object)){
        # Model fitting
        model_list[[b]] <- mgcv::gam(formula = object$fit[[b]]$formula, 
                                     family = object$fit[[b]]$family$family, 
                                     method = "REML",
                                     data = df_cat)
        add <- df_cat |>
          drop_na() |>
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
          drop_na() |>
          select({{ index_data }}, {{ key_data1 }})
        model_list[[b]]$model <- bind_cols(add, modelFrame)
        model_list[[b]]$model <- as_tsibble(model_list[[b]]$model,
                                            index = index_data,
                                            key = all_of(key_data1))
      }else if("pprFit" %in% class(object)){
        pre.formula <- object$fit[[b]]$terms
        attributes(pre.formula) <- NULL
        model_list[[b]] <- stats::ppr(formula = as.formula(pre.formula),
                                      data = df_cat)
        add <- df_cat |>
          drop_na() |>
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
    
    # Obtain predictions on validation set
    preds <- predict(object = modelFit[[i]], 
                     newdata = test, 
                     recursive = recursive, 
                     recursive_colRange = recursive_colRange)
    # Store predictions
    pf[i, ] <- as.numeric(preds$.predict)
    # Obtain in-sample residuals
    resids <- augment(x = modelFit[[i]])
    # Store residuls
    res[i, 1:length(resids$.resid)] <- as.numeric(resids$.resid)
    # Possible futures through block bootstrapping on in-sample residuals
    pFutures[[i]] <- smimodel::blockBootstrap(object = modelFit[[i]], newdata = test,
                                              resids = na.omit(res[i, ]), preds = pf[i, ], 
                                              season.period = season.period, m = m, 
                                              num.futures = num.futures, 
                                              recursive = recursive, 
                                              recursive_colRange = recursive_colRange)
    # Prediction interval bounds
    lower_q <- (1 - (level/100))/2
    upper_q <- lower_q + (level/100)
    for(j in 1:NROW(pFutures[[i]])){
      intervalsHilo <- quantile(pFutures[[i]][j, ], probs = c(lower_q, upper_q))
      nint <- length(level)
      lowerBound <- matrix(NA, ncol = nint, nrow = 1)
      upperBound <- lowerBound
      for (k in 1:nint) {
        lowerBound[1, k] <- intervalsHilo[[k]]
        upperBound[1, k] <- intervalsHilo[[(k + 2)]]
      }
      colnames(lowerBound) <- colnames(upperBound) <- paste(level, "%", sep = "")
      for (l in level) {
        levelname <- paste0(l, "%")
        lower[[levelname]][i, j] <- lowerBound[ , levelname]
        upper[[levelname]][i, j] <- upperBound[ , levelname]
      }
    }
  }
  
  out$method <- paste("crossVal")
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
  
  return(structure(out, class = "cvbb"))
}


#' Futures through single season block bootstrapping
#'
#' Gerenates possible future sample paths by applying the single season block
#' bootstrap method.
#'
#' @param object Fitted model object.
#' @param newdata Test data set. Must be a data set of class `tsibble`.
#' @param resids In-sample residuals from the fitted model.
#' @param preds Predictions for the test set (i.e. data for the forecast
#'   horizon).
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = `season.period` * `m`)
#' @param num.futures Number of possible future sample paths to be generated.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in test data to be filled with forecasts.
#'
#' @export
blockBootstrap <- function(object, newdata, resids, preds, season.period = 1, 
                           m = 1, num.futures = 1000, 
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
                                                     recursive_colRange = recursive_colRange)
    }else{
      futures <- possibleFutures_benchmark(object = object, 
                                                      newdata = newdata, 
                                                      bootstraps = bootstraps, 
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


#' Prediction intervals using single season block bootstrapping
#'
#' Generates multi-step prediction intervals corresponding to the forecasts
#' obtained on a given test set (i.e. `newdata`), using single season block
#' bootstrapping method.
#'
#' @param object A fitted model object of the class `smimodel`, `backward`,
#'   `pprFit`, or `gaimFit`.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = `season.period` * `m`)
#' @param num.futures Number of possible future sample paths to be generated.
#' @param confidence Confidence level.
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, the range of column numbers
#'   in `newdata` to be filled with forecasts.
#' 
#' @examples
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#' n = 1205
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y1, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Test set
#' sim_test <- sim_data[1001:1200, ]
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#' # Model fitting
#' model1 <- model_smimodel(data = sim_train,
#'                          yvar = "y1",
#'                          index.vars = index.vars,
#'                          initialise = "ppr")
#' # Calculating lower and upper bounds for 95% prediction interval
#' PI_model1 <- pi_bootstrap(object = model1, 
#'                          newdata = sim_test)
#' # Lower and upper bounds
#' PI_model1$bounds
#'                          
#' @references Hyndman, R.J. & Fan, S. (2010). Density Forecasting for Long-Term
#'   Peak Electricity Demand. *IEEE Transactions on Power Systems*, 25(2),
#'   1142–1153. \url{http://dx.doi.org/10.1109/TPWRS.2009.2036017}
#'
#' @export
pi_bootstrap <- function(object, newdata, season.period = 1, 
                         m = 1, num.futures = 1000, confidence = 95,
                         recursive = FALSE, recursive_colRange = NULL){
  # In-sample residuals
  x <- augment(x = object)$.resid
  # Generate the matrix of bootstrapped series
  bootstraps <- residBootstrap(x = x, season.period = season.period, m = m,
                               num.bootstrap = num.futures)
  # Predictions on test set
  preds <- predict(object = object, newdata = newdata, recursive = recursive, 
                   recursive_colRange = recursive_colRange)$.predict
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
                                          recursive_colRange = recursive_colRange)
    }else{
      futures <- possibleFutures_benchmark(object = object, 
                                           newdata = newdata, 
                                           bootstraps = bootstraps, 
                                           recursive_colRange = recursive_colRange)
    }
    names(futures$future_cols) <- 1:length(futures$future_cols)
    possibleFutures_part2 <- as.matrix(bind_cols(futures$future_cols))
    possibleFutures_part1 <- matrix(futures$firstFuture, nrow = 1, 
                                    ncol = length(futures$firstFuture))
    possibleFutures_mat <- rbind(possibleFutures_part1, possibleFutures_part2)
  }
  lower <- vector(mode = "list", length = NROW(possibleFutures_mat))
  upper <- vector(mode = "list", length = NROW(possibleFutures_mat))
  lower_q <- (1 - (confidence/100))/2
  upper_q <- lower_q + (confidence/100)
  for(j in 1:NROW(possibleFutures_mat)){
    intervalsHilo <- quantile(possibleFutures_mat[j, ], probs = c(lower_q, upper_q))
    # intervalsHilo <- distributional::hilo(distributional::dist_sample(list(possibleFutures_mat[j, ])), 
    #                                       confidence)
    lower[[j]] <- intervalsHilo[[1]]
    upper[[j]] <- intervalsHilo[[2]]
  }
  lower <- unlist(lower)
  upper <- unlist(upper)
  bounds <- tibble(
    predictions = preds,
    lower = lower,
    upper = upper
  )
  output <- list("bounds" = bounds, "sample_paths" = possibleFutures_mat)
  return(output)
}


#' Generate multiple single season block bootstrap series
#' 
#' Generates multiple replications of single season block bootstrap series.
#' 
#' @param x A series of residuals from which bootstrap series to be generated.
#' @param season.period Length of the seasonal period.
#' @param m Multiplier. (Block size = `season.period` * `m`)
#' @param num.bootstrap Number of bootstrap series to be generated.
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
#' @param m Multiplier. (Block size = `season.period` * `m`)

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

randomBlock <- function(series, block.size){
  start_ind <- sample(seq(length(series) - block.size + 1), 1)
  block <- series[start_ind:(start_ind + block.size -1)]
  return(block)
}



#' Possible future sample paths (multi-step) from `smimodel` residuals
#'
#' Generates possible future sample paths (multi-step) using residuals of a
#' fitted `smimodel` through recursive forecasting.
#'
#' @param object A `smimodel` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param bootstraps Generated matrix of bootstrapped residual series.
#' @param recursive_colRange The range of column numbers in `newdata` to be
#'   filled with forecasts.

possibleFutures_smimodel <- function(object, newdata, bootstraps, 
                                     recursive_colRange){
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  newdata <- newdata |>
    tibble::as_tibble() |>
    dplyr::arrange({{index_n}})
  predict_fn <- mgcv::predict.gam
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
  pred1 <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list1, type = "response")
  possibleFutures1 <- as.numeric(pred1) + bootstraps[1, ]
  # Should fill the missing response lags in newdata using each of the 
  # possible future values separately
  newdata_list <- vector(mode = "list", length = length(possibleFutures1))
  future_cols <- vector(mode = "list", length = length(possibleFutures1))
  for(s in 1:length(possibleFutures1)){
    newdata_list[[s]] <- newdata
    x_seq = seq((1+1), (1+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
    y_seq = recursive_colRange
    for(l in 1:length(recursive_colRange)){
      if((x_seq[l] <= NROW(newdata_list[[s]])) & is.na(newdata_list[[s]][x_seq[l], y_seq[l]])){
        newdata_list[[s]][x_seq[l], y_seq[l]] = possibleFutures1[s]
      }
    }
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
      pred1 <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list1, type = "response")
      temp_Futures[[t-1]] <- pred1 + bootstraps[t, s]
      x_seq = seq((t+1), (t+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(p in 1:length(recursive_colRange)){
        if((x_seq[p] <= NROW(newdata_list[[s]])) & is.na(newdata_list[[s]][x_seq[p], y_seq[p]])){
          newdata_list[[s]][x_seq[p], y_seq[p]] = temp_Futures[[t-1]]
        }
      }
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
    pred1 <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list1, 
                        type = "response")
    temp_Futures[[NROW(newdata_list[[s]]) - 1]] <- pred1 + bootstraps[NROW(newdata_list[[s]]), s]
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
#' @param object A fitted model object of the class `backward`, `pprFit`, or
#'   `gaimFit`.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param bootstraps Generated matrix of bootstrapped residual series.
#' @param recursive_colRange The range of column numbers in `newdata` to be
#'   filled with forecasts.

possibleFutures_benchmark <- function(object, newdata, bootstraps, 
                                      recursive_colRange){
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  newdata <- newdata |>
    tibble::as_tibble() |>
    dplyr::arrange({{index_n}})
  # Recursive possible futures
  # First, one-step-ahead possible futures
  data_temp = newdata[1, ]
  key22 = data_temp[ , {{ key11 }}][[1]]
  key22_pos = which(object$key == key22)
  pred1 <- predict(object$fit[[key22_pos]], data_temp, type = "response")
  possibleFutures1 <- as.numeric(pred1) + bootstraps[1, ]
  # Should fill the missing response lags in newdata using each of the 
  # possible future values separately
  newdata_list <- vector(mode = "list", length = length(possibleFutures1))
  future_cols <- vector(mode = "list", length = length(possibleFutures1))
  for(s in 1:length(possibleFutures1)){
    newdata_list[[s]] <- newdata
    x_seq = seq((1+1), (1+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
    y_seq = recursive_colRange
    for(l in 1:length(recursive_colRange)){
      if((x_seq[l] <= NROW(newdata_list[[s]])) & is.na(newdata_list[[s]][x_seq[l], y_seq[l]])){
        newdata_list[[s]][x_seq[l], y_seq[l]] = possibleFutures1[s]
      }
    }
    temp_Futures <- vector(mode = "list", length = length(2:NROW(newdata_list[[s]])))
    for(t in 2:(NROW(newdata_list[[s]]) - 1)){
      data_temp = newdata_list[[s]][t, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      pred1 <- predict(object$fit[[key22_pos]], data_temp, type = "response")
      temp_Futures[[t-1]] <- pred1 + bootstraps[t, s]
      x_seq = seq((t+1), (t+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata_list[[s]])) & is.na(newdata_list[[s]][x_seq[l], y_seq[l]])){
          newdata_list[[s]][x_seq[l], y_seq[l]] = temp_Futures[[t-1]]
        }
      }
    }
    data_temp = newdata_list[[s]][NROW(newdata_list[[s]]), ]
    key22 = data_temp[ , {{ key11 }}][[1]]
    key22_pos = which(object$key == key22)
    pred1 <- predict(object$fit[[key22_pos]], data_temp, type = "response")
    temp_Futures[[NROW(newdata_list[[s]]) - 1]] <- pred1 + bootstraps[NROW(newdata_list[[s]]), s]
    future_cols[[s]] <- unlist(temp_Futures)
  }
  output <- list("firstFuture" = possibleFutures1, "future_cols" = future_cols)
  return(output)
}


#' Calculate interval forecast coverage
#'
#' This is a wrapper for the function `conformalForecast::coverage()`.
#' Calculates the mean coverage and the ifinn matrix for prediction intervals on
#' validation set. If \code{window} is not \code{NULL}, a matrix of the rolling
#' means of interval forecast coverage is also returned.
#'
#' @param object An object of class \code{cvbb}.
#' @param level Target confidence level for prediction intervals.
#' @param window If not \code{NULL}, the rolling mean matrix for coverage is
#'   also returned.
#' @param na.rm A logical indicating whether \code{NA} values should be stripped
#'   before the rolling mean computation proceeds.
#'
#' @return A list of class \code{coverage} with the following components:
#'   \item{mean}{Mean coverage across the validation set.}
#' \item{ifinn}{A indicator matrix as a multivariate time series, where the \eqn{h}th column
#' holds the coverage for forecast horizon \eqn{h}. The time index
#' corresponds to the period for which the forecast is produced.}
#' \item{rollmean}{If \code{window} is not NULL, a matrix of the rolling means
#' of interval forecast coverage will be returned.}
#'
#' @export
avgCoverage <- function(object, level = 95, window = NULL, na.rm = FALSE) {
  object$LOWER <- object$lower
  object$UPPER <- object$upper
  output <- conformalForecast::coverage(object = object, level = level,
                                        window = window, na.rm = na.rm)
  return(output)
}


#' Create lags or leads of a matrix
#'
#' This is a wrapper for the function `conformalForecast::lagmatrix()`.Find a
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
  lmat <- conformalForecast::lagmatrix(x = x, lag = lag)
  return(lmat)
}
