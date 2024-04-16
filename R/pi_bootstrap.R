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
#' @importFrom stats quantile
#' @importFrom tibble tibble
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
  return(bounds)
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
#' 
#' @importFrom stats na.omit
#' 
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
    newdata <- newdata %>%
      mutate(dummy_key = rep(1, NROW(newdata))) %>%
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  newdata <- newdata %>%
    tibble::as_tibble() %>%
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
    newdata <- newdata %>%
      mutate(dummy_key = rep(1, NROW(newdata))) %>%
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  newdata <- newdata %>%
    tibble::as_tibble() %>%
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



#' Coverage of a calculated prediction interval
#' 
#' Computes the actual coverage probability of a calculated prediction interval.
#' 
#' @param bounds The `tibble` returned from `pi_bootstrap()`.
#' @param newdata The set of new data on for which the forecasts are required 
#' (i.e. test set; should be a `tsibble`).
#' @param yvar Response variable as a character string.
#' 
#' @export
coverage <- function(bounds, newdata, yvar){
  actual <- newdata[ , {{yvar}}][[1]]
  within_count <- 0
  for(i in 1:length(actual)){
    if((actual[i] > bounds$lower[i]) & (actual[i] < bounds$upper[i])){
      within_count <- within_count + 1
    }
  }
  cover <- within_count/length(actual)
  return(cover)
}