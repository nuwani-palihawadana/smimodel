#' Obtaining forecasts on a test set from a fitted \code{smimodel}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{smimodel} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tsibble} with forecasts on test set.
#'
#' @method predict smimodel
#'
#' @export
predict.smimodel <- function(object, newdata, exclude.trunc = NULL, 
                             recursive = FALSE,
                             recursive_colRange = NULL, ...) {
  if (!tsibble::is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      dplyr::mutate(dummy_key = rep(1, NROW(newdata))) |>
      tsibble::as_tsibble(index = index_n, key = dummy_key)
  }
  key_n <- key(newdata)[[1]]
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    data_list <- vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key_n }}][[1]]
      key22_pos = which(object$key == key22)
      object_temp <- object$fit[[key22_pos]]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      gam_cols <- c(object_temp$best$vars_index, object_temp$best$vars_s)
      trunc_indx <- !(gam_cols %in% exclude.trunc)
      trunc_cols <- gam_cols[trunc_indx]
      # Truncate
      if(length(trunc_cols != 0)){
        data_temp <- truncate_vars(range.object = object_temp$best$vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      list_index <- object_temp$best$alpha
      numInd <- NCOL(list_index)
      alpha <- vector(mode = "list", length = numInd)
      for(b in 1:numInd){
        alpha[[b]] <- list_index[ , b]
      }
      alpha <- unlist(alpha)
      names(alpha) <- NULL
      if(all(alpha == 0)){
        data_list[[m]] <- data_temp
      }else{
        X_test <- as.matrix(data_temp[ , object_temp$best$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = numInd)
        for(i in 1:numInd){
          ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
        }
        names(ind) <- colnames(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[m]] <- dplyr::bind_cols(data_temp, dat)
      }
      
      # Prediction
      pred <- predict(object_temp$best$gam, data_list[[m]], type = "response")
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    key22 = data_temp[ , {{ key_n }}][[1]]
    key22_pos = which(object$key == key22)
    object_temp <- object$fit[[key22_pos]]
    ## Avoid extrapolation; truncate non-linear predictors to match in-sample
    ## range
    # Predictors to truncate
    gam_cols <- c(object_temp$best$vars_index, object_temp$best$vars_s)
    trunc_indx <- !(gam_cols %in% exclude.trunc)
    trunc_cols <- gam_cols[trunc_indx]
    # Truncate
    if(length(trunc_cols != 0)){
      data_temp <- truncate_vars(range.object = object_temp$best$vars_range,
                                 data = data_temp,
                                 cols.trunc = trunc_cols)
    }
    list_index <- object_temp$best$alpha
    numInd <- NCOL(list_index)
    alpha <- vector(mode = "list", length = numInd)
    for(w in 1:numInd){
      alpha[[w]] <- list_index[ , w]
    }
    alpha <- unlist(alpha)
    names(alpha) <- NULL
    if(all(alpha == 0)){
      data_list[[NROW(newdata)]] <- data_temp
    }else{
      X_test <- as.matrix(data_temp[ , object_temp$best$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = numInd)
      for(a in 1:numInd){
        ind[[a]] <- as.numeric(X_test %*% as.matrix(list_index[ , a], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list[[NROW(newdata)]] <- dplyr::bind_cols(data_temp, dat)
    }
    # Prediction
    pred <- predict(object_temp$best$gam, data_list[[NROW(newdata)]], type = "response")
    predictions[[NROW(newdata)]] <- pred
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    data_list <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key_n }}] == object$key[i], ]
      object_temp <- object$fit[[i]]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      gam_cols <- c(object_temp$best$vars_index, object_temp$best$vars_s)
      trunc_indx <- !(gam_cols %in% exclude.trunc)
      trunc_cols <- gam_cols[trunc_indx]
      # Truncate
      if(length(trunc_cols != 0)){
        newdata_cat <- truncate_vars(range.object = object_temp$best$vars_range,
                                     data = newdata_cat,
                                     cols.trunc = trunc_cols)
      }
      list_index <- object$fit[[i]]$best$alpha
      numInd <- NCOL(list_index)
      alpha <- vector(mode = "list", length = numInd)
      for(z in 1:numInd){
        alpha[[z]] <- list_index[ , z]
      }
      alpha <- unlist(alpha)
      names(alpha) <- NULL
      if(all(alpha == 0)){
        data_list[[i]] <- newdata_cat
      }else{
        X_test <- as.matrix(newdata_cat[ , object$fit[[i]]$best$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = numInd)
        for(k in 1:length(ind)){
          ind[[k]] <- as.numeric(X_test %*% as.matrix(list_index[ , k], ncol = 1))
        }
        names(ind) <- colnames(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[i]] <- dplyr::bind_cols(newdata_cat, dat)
      }
      # Prediction
      predictions[[i]] <- predict(object$fit[[i]]$best$gam, data_list[[i]],
                                     type = "response")
    }
  }
  newdata1 <- dplyr::bind_rows(data_list)
  pred <- unlist(predictions)
  pred_F <- newdata1 |>
    dplyr::mutate(.predict = pred) |>
    tsibble::as_tsibble(index = {{ index_n }}, key = {{ key_n }})
  return(pred_F)
}


#' Obtaining forecasts on a test set from a \code{smimodelFit}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{smimodelFit} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble} with forecasts on test set.
#'
#' @method predict smimodelFit
#' @export
predict.smimodelFit <- function(object, newdata, exclude.trunc = NULL, 
                                recursive = FALSE,
                                recursive_colRange = NULL, ...) {
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  list_index <- object$alpha
  num_ind <- NCOL(list_index)
  alpha <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    alpha[[i]] <- list_index[ , i]
  }
  alpha <- unlist(alpha)
  names(alpha) <- NULL
  # Predictors to truncate
  gam_cols <- c(object$vars_index, object$vars_s)
  trunc_indx <- !(gam_cols %in% exclude.trunc)
  trunc_cols <- gam_cols[trunc_indx]
  if(all(alpha == 0)){
    if(recursive == TRUE){
      # Prepare newdata for recursive forecasting
      newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
      # Recursive forecasting
      predictions =  vector(mode = "list", length = NROW(newdata))
      for(m in 1:(NROW(newdata) - 1)){
        data_temp = newdata[m, ]
        if(length(trunc_cols != 0)){
          ## Avoid extrapolation; truncate non-linear predictors to match 
          ## in-sample range
          data_temp <- truncate_vars(range.object = object$vars_range,
                                     data = data_temp,
                                     cols.trunc = trunc_cols)
        }
        pred <- predict(object$gam, data_temp, type = "response")
        predictions[[m]] <- pred
        x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
        y_seq = recursive_colRange
        for(l in 1:length(recursive_colRange)){
          if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
            newdata[x_seq[l], y_seq[l]] = pred
          }
        }
      }
      data_temp = newdata[NROW(newdata), ]
      if(length(trunc_cols != 0)){
        ## Avoid extrapolation; truncate non-linear predictors to match in-sample
        data_temp <- truncate_vars(range.object = object$vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      predictions[[NROW(newdata)]] = predict(object$gam, data_temp, type = "response")
      pred <- unlist(predictions)
      pred_F <- newdata |>
        dplyr::mutate(.predict = pred)
    }else if(recursive == FALSE){
      # Index
      index_data <- index(newdata)
      # Convert to a tibble
      newdata <- newdata |>
        tibble::as_tibble() |>
        dplyr::arrange({{index_data}})
      if(length(trunc_cols != 0)){
        newdata <- truncate_vars(range.object = object$vars_range,
                                 data = newdata,
                                 cols.trunc = trunc_cols)
      }
      pred <- predict(object$gam, newdata, type = "response")
      pred_F <- newdata |>
        dplyr::mutate(.predict = pred)
    }
  }else{
    if(recursive == TRUE){
      # Prepare newdata for recursive forecasting
      newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
      predictions =  vector(mode = "list", length = NROW(newdata))
      data_list <- vector(mode = "list", length = NROW(newdata))
      for(m in 1:(NROW(newdata) - 1)){
        data_temp = newdata[m, ]
        if(length(trunc_cols != 0)){
          ## Avoid extrapolation; truncate non-linear predictors to match 
          ## in-sample range
          data_temp <- truncate_vars(range.object = object$vars_range,
                                     data = data_temp,
                                     cols.trunc = trunc_cols)
        }
        X_test <- as.matrix(data_temp[ , object$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
        }
        names(ind) <- colnames(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[m]] <- dplyr::bind_cols(data_temp, dat)
        pred <- predict(object$gam, data_list[[m]], type = "response")
        predictions[[m]] <- pred
        x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
        y_seq = recursive_colRange
        for(l in 1:length(recursive_colRange)){
          if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
            newdata[x_seq[l], y_seq[l]] = pred
          }
        }
      }
      data_temp = newdata[NROW(newdata), ]
      if(length(trunc_cols != 0)){
        ## Avoid extrapolation; truncate non-linear predictors to match 
        ## in-sample range
        data_temp <- truncate_vars(range.object = object$vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      X_test <- as.matrix(data_temp[ , object$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = num_ind)
      for(i in 1:num_ind){
        ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list[[NROW(newdata)]] <- dplyr::bind_cols(data_temp, dat)
      predictions[[NROW(newdata)]] = predict(object$gam, data_list[[NROW(newdata)]], type = "response")
      newdata1 <- dplyr::bind_rows(data_list)
      pred <- unlist(predictions)
      pred_F <- newdata1 |>
        dplyr::mutate(.predict = pred)
    }else if(recursive == FALSE){
      # Index
      index_data <- index(newdata)
      # Convert to a tibble
      newdata <- newdata |>
        tibble::as_tibble() |>
        dplyr::arrange({{index_data}})
      if(length(trunc_cols != 0)){
        ## Avoid extrapolation; truncate non-linear predictors to match 
        ## in-sample range
        newdata <- truncate_vars(range.object = object$vars_range,
                                 data = newdata,
                                 cols.trunc = trunc_cols)
      }
      X_test <- as.matrix(newdata[ , object$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = num_ind)
      for(i in 1:num_ind){
        ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list <- dplyr::bind_cols(newdata, dat)
      pred <- predict(object$gam, data_list, type = "response")
      pred_F <- data_list |>
        dplyr::mutate(.predict = pred)
    }
  }
  return(pred_F)
}

#' Obtaining forecasts on a test set from a fitted \code{backward}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{backward} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tsibble} with forecasts on test set.
#'
#' @method predict backward
#'
#' @export
predict.backward <- function(object, newdata, exclude.trunc = NULL,
                             recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key_n <- key(newdata)[[1]]
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key_n }}][[1]]
      key22_pos = which(object$key == key22)
      object_temp <- object$fit[[key22_pos]]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      gam_cols <- colnames(object_temp$model)
      remove_temp <- as.character(attributes(object_temp$pterms)$variables)
      no_trunc_cols <- unique(c(as.character(index_n), as.character(key_n),
                                remove_temp[2:length(remove_temp)], 
                                exclude.trunc))
      trunc_indx <- !(gam_cols %in% no_trunc_cols)
      trunc_cols <- gam_cols[trunc_indx]
      # In-sample range 
      vars_original <- object_temp$model[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      # Truncate
      if(length(trunc_cols) != 0){
        data_temp <- truncate_vars(range.object = vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred <- predict(object_temp, data_temp, type = "response")
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    key22 = data_temp[ , {{ key_n }}][[1]]
    key22_pos = which(object$key == key22)
    object_temp <- object$fit[[key22_pos]]
    ## Avoid extrapolation; truncate non-linear predictors to match in-sample
    ## range
    # Predictors to truncate
    gam_cols <- colnames(object_temp$model)
    remove_temp <- as.character(attributes(object_temp$pterms)$variables)
    no_trunc_cols <- unique(c(as.character(index_n), as.character(key_n),
                              remove_temp[2:length(remove_temp)], 
                              exclude.trunc))
    trunc_indx <- !(gam_cols %in% no_trunc_cols)
    trunc_cols <- gam_cols[trunc_indx]
    # In-sample range 
    vars_original <- object_temp$model[ , trunc_cols]
    vars_range <- apply(vars_original, 2, range)
    # Truncate
    if(length(trunc_cols) != 0){
      data_temp <- truncate_vars(range.object = vars_range,
                                 data = data_temp,
                                 cols.trunc = trunc_cols)
    }
    predictions[[NROW(newdata)]] = predict(object_temp, data_temp, type = "response")
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key_n }}] == object$key[i], ]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      gam_cols <- colnames(object$fit[[i]]$model)
      remove_temp <- as.character(attributes(object$fit[[i]]$pterms)$variables)
      no_trunc_cols <- unique(c(as.character(index_n), as.character(key_n),
                                remove_temp[2:length(remove_temp)], 
                                exclude.trunc))
      trunc_indx <- !(gam_cols %in% no_trunc_cols)
      trunc_cols <- gam_cols[trunc_indx]
      # In-sample range 
      vars_original <- object$fit[[i]]$model[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      if(length(trunc_cols) != 0){
        newdata_cat <- truncate_vars(range.object = vars_range,
                                     data = newdata_cat,
                                     cols.trunc = trunc_cols)
      }
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
  }
  pred <- unlist(predictions)
  pred_F <- newdata |>
    mutate(.predict = pred) |>
    as_tsibble(index = {{ index_n }}, key = {{ key_n }})
  return(pred_F)
}


#' Obtaining forecasts on a test set from a fitted \code{pprFit}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{pprFit} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tsibble} with forecasts on test set.
#'
#' @method predict pprFit
#'
#' @export
predict.pprFit <- function(object, newdata, exclude.trunc = NULL,
                           recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key_n <- key(newdata)[[1]]
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      col_retain <- !(colnames(data_temp) %in% c("indexLag", "indexDiff", "row_idx", "grp"))
      if(any(is.na(data_temp[ , col_retain]))){
        pred <- NA
      }else{
        key22 = data_temp[ , {{ key_n }}][[1]]
        key22_pos = which(object$key == key22)
        object_temp <- object$fit[[key22_pos]]
        ## Avoid extrapolation; truncate non-linear predictors to match 
        ## in-sample range
        # Predictors to truncate
        trunc_indx <- !(object_temp$xnames %in% exclude.trunc)
        trunc_cols <- object_temp$xnames[trunc_indx]
        # In-sample range 
        vars_original <- object_temp$model[ , trunc_cols]
        vars_range <- apply(vars_original, 2, range)
        if(length(trunc_cols) != 0){
          data_temp <- truncate_vars(range.object = vars_range,
                                     data = data_temp,
                                     cols.trunc = trunc_cols)
        }
        pred <- predict(object_temp, data_temp, type = "response")
      }
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    if(any(is.na(data_temp[ , 1:(NCOL(data_temp) - 2)]))){
      pred <- NA
    }else{
      key22 = data_temp[ , {{ key_n }}][[1]]
      key22_pos = which(object$key == key22)
      object_temp <- object$fit[[key22_pos]]
      ## Avoid extrapolation; truncate non-linear predictors to match 
      ## in-sample range
      # Predictors to truncate
      trunc_indx <- !(object_temp$xnames %in% exclude.trunc)
      trunc_cols <- object_temp$xnames[trunc_indx]
      # In-sample range 
      vars_original <- object_temp$model[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      if(length(trunc_cols) != 0){
        data_temp <- truncate_vars(range.object = vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred <- predict(object_temp, data_temp, type = "response")
    }
    predictions[[NROW(newdata)]] <- pred
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key_n }}] == object$key[i], ]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      trunc_indx <- !(object$fit[[i]]$xnames %in% exclude.trunc)
      trunc_cols <- object$fit[[i]]$xnames[trunc_indx]
      # In-sample range 
      vars_original <- object$fit[[i]]$model[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      if(length(trunc_cols) != 0){
        newdata_cat <- truncate_vars(range.object = vars_range,
                                     data = newdata_cat,
                                     cols.trunc = trunc_cols)
      }
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
  }
  pred <- unlist(predictions)
  pred_F <- newdata |>
    mutate(.predict = pred) |>
    as_tsibble(index = {{ index_n }}, key = {{ key_n }})
  return(pred_F)
}


#' Obtaining forecasts on a test set from a fitted \code{gaimFit}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{gaimFit} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tsibble} with forecasts on test set.
#'
#' @method predict gaimFit
#'
#' @export
predict.gaimFit <- function(object, newdata, exclude.trunc = NULL,
                            recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key_n <- key(newdata)[[1]]
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      #print(paste0("This is newdata", m))
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key_n }}][[1]]
      key22_pos = which(object$key == key22)
      object_temp <- object$fit[[key22_pos]]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      cgaim_cols <- c(object_temp$vars_index, object_temp$vars_s)
      trunc_indx <- !(cgaim_cols %in% exclude.trunc)
      trunc_cols <- cgaim_cols[trunc_indx]
      # In-sample range 
      vars_original <- bind_cols(object_temp$x, object_temp$sm_mod$Xcov)[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      # Truncate
      if(length(trunc_cols) != 0){
        data_temp <- truncate_vars(range.object = vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred <- predict(object_temp, data_temp, type = "response")
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    key22 = data_temp[ , {{ key_n }}][[1]]
    key22_pos = which(object$key == key22)
    object_temp <- object$fit[[key22_pos]]
    ## Avoid extrapolation; truncate non-linear predictors to match in-sample
    ## range
    # Predictors to truncate
    cgaim_cols <- c(object_temp$vars_index, object_temp$vars_s)
    trunc_indx <- !(cgaim_cols %in% exclude.trunc)
    trunc_cols <- cgaim_cols[trunc_indx]
    # In-sample range 
    vars_original <- bind_cols(object_temp$x, object_temp$sm_mod$Xcov)[ , trunc_cols]
    vars_range <- apply(vars_original, 2, range)
    # Truncate
    if(length(trunc_cols) != 0){
      data_temp <- truncate_vars(range.object = vars_range,
                                 data = data_temp,
                                 cols.trunc = trunc_cols)
    }
    predictions[[NROW(newdata)]] = predict(object_temp, data_temp, type = "response")
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key_n }}] == object$key[i], ]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      cgaim_cols <- c(object$fit[[i]]$vars_index, object$fit[[i]]$vars_s)
      trunc_indx <- !(cgaim_cols %in% exclude.trunc)
      trunc_cols <- cgaim_cols[trunc_indx]
      # In-sample range 
      vars_original <- bind_cols(object$fit[[i]]$x, object$fit[[i]]$sm_mod$Xcov)[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      # Truncate
      if(length(trunc_cols) != 0){
        newdata_cat <- truncate_vars(range.object = vars_range,
                                     data = newdata_cat,
                                     cols.trunc = trunc_cols)
      }
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
  }
  pred <- unlist(predictions)
  pred_F <- newdata |>
    mutate(.predict = pred) |>
    as_tsibble(index = {{ index_n }}, key = {{ key_n }})
  return(pred_F)
}


#' Obtaining forecasts on a test set from a fitted \code{lmFit}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{lmFit} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tsibble} with forecasts on test set.
#'
#' @method predict lmFit
#'
#' @export
predict.lmFit <- function(object, newdata,
                           recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key11 <- key(newdata)[[1]]
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      pred <- predict(object$fit[[key22_pos]], data_temp, type = "response")
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    key22 = data_temp[ , {{ key11 }}][[1]]
    key22_pos = which(object$key == key22)
    predictions[[NROW(newdata)]] = predict(object$fit[[key22_pos]], data_temp,
                                           type = "response")
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key11 }}] == object$key[i], ]
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
  }
  pred <- unlist(predictions)
  pred_F <- newdata |>
    mutate(.predict = pred) |>
    as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  return(pred_F)
}


#' Obtaining forecasts on a test set from a fitted \code{gamFit}
#'
#' Gives forecasts on a test set.
#'
#' @param object A \code{gamFit} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tsibble} with forecasts on test set.
#'
#' @method predict gamFit
#'
#' @export
predict.gamFit <- function(object, newdata, exclude.trunc = NULL,
                          recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata |>
      mutate(dummy_key = rep(1, NROW(newdata))) |>
      as_tsibble(index = index_n, key = dummy_key)
  }
  key_n <- key(newdata)[[1]]
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key_n }}][[1]]
      key22_pos = which(object$key == key22)
      object_temp <- object$fit[[key22_pos]]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      gam_cols <- colnames(object_temp$model)
      remove_temp <- as.character(attributes(object_temp$pterms)$variables)
      no_trunc_cols <- unique(c(as.character(index_n), as.character(key_n),
                                remove_temp[2:length(remove_temp)], 
                                exclude.trunc))
      trunc_indx <- !(gam_cols %in% no_trunc_cols)
      trunc_cols <- gam_cols[trunc_indx]
      # In-sample range
      vars_original <- object_temp$model[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      if(length(trunc_cols) != 0){
        data_temp <- truncate_vars(range.object = vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred <- predict(object_temp, data_temp, type = "response")
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    key22 = data_temp[ , {{ key_n }}][[1]]
    key22_pos = which(object$key == key22)
    object_temp <- object$fit[[key22_pos]]
    ## Avoid extrapolation; truncate non-linear predictors to match in-sample
    ## range
    # Predictors to truncate
    gam_cols <- colnames(object_temp$model)
    remove_temp <- as.character(attributes(object_temp$pterms)$variables)
    no_trunc_cols <- unique(c(as.character(index_n), as.character(key_n),
                              remove_temp[2:length(remove_temp)], 
                              exclude.trunc))
    trunc_indx <- !(gam_cols %in% no_trunc_cols)
    trunc_cols <- gam_cols[trunc_indx]
    # In-sample range
    vars_original <- object_temp$model[ , trunc_cols]
    vars_range <- apply(vars_original, 2, range)
    if(length(trunc_cols) != 0){
      data_temp <- truncate_vars(range.object = vars_range,
                                 data = data_temp,
                                 cols.trunc = trunc_cols)
    }
    predictions[[NROW(newdata)]] = predict(object_temp, data_temp, type = "response")
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key_n }}] == object$key[i], ]
      ## Avoid extrapolation; truncate non-linear predictors to match in-sample
      ## range
      # Predictors to truncate
      gam_cols <- colnames(object$fit[[i]]$model)
      remove_temp <- as.character(attributes(object$fit[[i]]$pterms)$variables)
      no_trunc_cols <- unique(c(as.character(index_n), as.character(key_n),
                                remove_temp[2:length(remove_temp)], 
                                exclude.trunc))
      trunc_indx <- !(gam_cols %in% no_trunc_cols)
      trunc_cols <- gam_cols[trunc_indx]
      # In-sample range
      vars_original <- object$fit[[i]]$model[ , trunc_cols]
      vars_range <- apply(vars_original, 2, range)
      if(length(trunc_cols) != 0){
        newdata_cat <- truncate_vars(range.object = vars_range,
                                     data = newdata_cat,
                                     cols.trunc = trunc_cols)
      }
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
  }
  pred <- unlist(predictions)
  pred_F <- newdata |>
    mutate(.predict = pred) |>
    as_tsibble(index = {{ index_n }}, key = {{ key_n }})
  return(pred_F)
}


#' Obtaining recursive forecasts on a test set from a fitted \code{mgcv::gam}
#'
#' Gives recursive forecasts on a test set.
#'
#' @param object A \code{gam} object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a \code{tsibble}).
#' @param exclude.trunc The names of the predictor variables that should not be
#'   truncated for stable predictions as a character string. (Since the
#'   nonlinear functions are estimated using splines, extrapolation is not
#'   desirable. Hence, if any predictor variable in `newdata` that is treated
#'   non-linearly in the estimated model, will be truncated to be in the
#'   in-sample range before obtaining predictions. If any variables are listed
#'   here will be excluded from such truncation.)
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   \code{FALSE}).
#' @param recursive_colRange If \code{recursive = TRUE}, the range of column
#'   numbers in \code{newdata} to be filled with forecasts.
#'   Recursive/autoregressive forecasting is required when the lags of the
#'   response variable itself are used as predictor variables into the model.
#'   Make sure such lagged variables are positioned together in increasing lag
#'   order (i.e. \code{lag_1, lag_2, ..., lag_m}, \code{lag_m =} maximum lag
#'   used) in \code{newdata}, with no break in the lagged variable sequence even
#'   if some of the intermediate lags are not used as predictors.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble} with forecasts on test set.
predict_gam <- function(object, newdata, exclude.trunc = NULL, 
                          recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  # Predictors to truncate
  gam_cols <- colnames(object$model)
  remove_temp <- as.character(attributes(object$pterms)$variables)
  no_trunc_cols <- unique(c(as.character(index_n),
                            remove_temp[2:length(remove_temp)], 
                            exclude.trunc))
  trunc_indx <- !(gam_cols %in% no_trunc_cols)
  trunc_cols <- gam_cols[trunc_indx]
  # In-sample range 
  vars_original <- object$model[ , trunc_cols]
  vars_range <- apply(vars_original, 2, range)
  if(recursive == TRUE){
    # Prepare newdata for recursive forecasting
    newdata <- prep_newdata(newdata = newdata, recursive_colRange = recursive_colRange)
    # Recursive forecasting
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      if(length(trunc_cols) != 0){
        ## Avoid extrapolation; truncate non-linear predictors to match 
        ## in-sample range
        data_temp <- truncate_vars(range.object = vars_range,
                                   data = data_temp,
                                   cols.trunc = trunc_cols)
      }
      pred <- predict(object, data_temp, type = "response")
      predictions[[m]] <- pred
      x_seq = seq((m+1), (m+((max(recursive_colRange) - min(recursive_colRange)) + 1)))
      y_seq = recursive_colRange
      for(l in 1:length(recursive_colRange)){
        if((x_seq[l] <= NROW(newdata)) & is.na(newdata[x_seq[l], y_seq[l]])){
          newdata[x_seq[l], y_seq[l]] = pred
        }
      }
    }
    data_temp = newdata[NROW(newdata), ]
    if(length(trunc_cols) != 0){
      ## Avoid extrapolation; truncate non-linear predictors to match 
      ## in-sample range
      data_temp <- truncate_vars(range.object = vars_range,
                                 data = data_temp,
                                 cols.trunc = trunc_cols)
    }
    predictions[[NROW(newdata)]] = predict(object, data_temp, type = "response")
    pred <- unlist(predictions)
  }else if(recursive == FALSE){
    # Index
    index_data <- index(newdata)
    # Convert to a tibble
    newdata <- newdata |>
      tibble::as_tibble() |>
      dplyr::arrange({{index_data}})
    if(length(trunc_cols) != 0){
      ## Avoid extrapolation; truncate non-linear predictors to match 
      ## in-sample range
      newdata <- truncate_vars(range.object = vars_range,
                               data = newdata,
                               cols.trunc = trunc_cols)
    }
    pred <- predict(object, newdata, type = "response")
  }
  pred_F <- newdata |>
    dplyr::mutate(.predict = pred)
  return(pred_F)
}


#' Truncating predictors to be in the in-sample range
#'
#' Truncates predictors to be in the in-sample range to avoid spline 
#' extrapolation.
#'
#' @param range.object A matrix containing range of each predictor variable. 
#' Should be a matrix with two rows for min and max, and the columns should 
#' correspond to variables.
#' @param data Out-of-sample data set of which variables should be truncated.
#' @param cols.trunc Column names of the variables to be truncated.
truncate_vars <- function(range.object, data, cols.trunc){
  for(i in 1:length(cols.trunc)){
    # In-sample range
    insample_range <- range.object[ , cols.trunc[i]]
    # Truncate
    if(NROW(data) == 1){
      if(data[ , cols.trunc[i]] < insample_range[1]){
        # Less than minimum
        data[ , cols.trunc[i]] <- insample_range[1]
      }else if(data[ , cols.trunc[i]] > insample_range[2]){
        # Greater than maximum
        data[ , cols.trunc[i]] <- insample_range[2]
      }
    }else{
      # Less than minimum
      smaller_indx <- which(data[ , cols.trunc[i]] < insample_range[1])
      # Greater than maximum
      larger_indx <- which(data[ , cols.trunc[i]] > insample_range[2])
      # Truncate
      data[smaller_indx, cols.trunc[i]] <- insample_range[1]
      data[larger_indx, cols.trunc[i]] <- insample_range[2]
    }
  }
  return(data)
}


#' Predictions from a fitted CGAIM object - copied from cgaim:::predict.cgaim()
#' and modified
#'
#' Uses a fitted \code{cgaim} object and computes prediction for the observed
#' data or new data. Predicts the response, indices or ridge functions values at
#' the provided data.
#'
#' @param object A \code{gaim} object.
#' @param newdata A list or data.frame containing the new data to predict. If
#'   missing, fitted values from the model are returned.
#' @param type A character indicating the type of prediction to return.
#'   \code{type = "response"} returns the predicted response. \code{type =
#'   "terms"}, returns ridge and smooth functions evaluated at index predicted
#'   for \code{newdata}. \code{type = "scterms"} is the same, except that terms
#'   are postmultiplied by their scaling coefficients beta. \code{type =
#'   "indices"} returns predicted indices values.
#' @param select A numeric or character vector indicating terms to return for
#'   all types except \code{"response"}.
#' @param na.action A function indicating how to treat NAs. See
#'   \code{\link[stats]{na.fail}}.
#' @param ... For compatibility with the default \code{predict} method. Unused
#'   at the moment.
#'
#' @details \code{type = "terms"} returns the scaled ridge functions, i.e.
#'   before being multiplied by scaling coefficients beta.
#'
#' @return When \code{type = "response"} returns a vector of predicted response.
#'   When \code{type = "terms"} or \code{"scterms"}, returns a matrix of
#'   evaluated ridge and smooth terms. When \code{type = "indices"}, returns a
#'   matrix of evaluated indices.
#'
#' @examples
#' ## Load library
#' library(cgaim)
#' ## Simulate some data
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' mu <- 4 * exp(8 * x1) / (1 + exp(8 * x1)) + exp(x3)
#' y <- mu + rnorm(n)
#' df1 <- data.frame(y, x1, x2, x3, x4)
#'
#' ## Fit an unconstrained the model
#' ans <- cgaim(y ~ g(x1, x2, label = "foo") + g(x3, x4, label = "bar"),
#'   data = df1)
#'
#' ## Get fitted values
#' yhat <- predict(ans)
#'
#' ## Predict on new data
#' newdf <- as.data.frame(matrix(rnorm(100), 25, 4))
#' names(newdf) <- sprintf("x%i", 1:4)
#'
#' # predicted response
#' ypred <- predict(ans, newdf)
#'
#' # Indices
#' indices <- predict(ans, newdata = newdf, type = "indices")
#'
#' # Ridge functions
#' funs <- predict(ans, newdata = newdf, type = "terms")
#'
#' ## Select specific terms
#' ind1 <- predict(ans, newdata = newdf, select = "foo", type = "indices")
#' fun1 <- predict(ans, newdata = newdf, select = "foo", type = "terms")
#'
#' # Plot
#' plot(ans, select = "foo")
#' points(ind1, fun1)
#'
#' ## Scaled terms
#' fun2 <- predict(ans, newdata = newdf, select = "foo", type = "scterms")
#'
#' # Plot
#' plot(ans, select = "foo", yscale = TRUE)
#' points(ind1, fun2)
#'
#' @export
predict.cgaim <- function(object, newdata,
                          type = c("response", "terms", "scterms", "indices"),
                          select = NULL,
                          na.action = "na.pass", ...)
{
  type <- match.arg(type)
  n <- length(object$fitted)
  p <- length(object$beta) - 1
  # Check select
  if (is.null(select)){
    select <- seq_len(p)
  } else {
    if (is.character(select)){
      selmatch <- match(select, colnames(object$gfit))
      nas <- is.na(selmatch)
      if (any(nas)) warning(paste0("Incorrect select removed: ",
                                   paste(select[nas], collapse = ", ")))
      select <- selmatch[!nas]
    } else {
      inc <- select > p
      if (any(inc)) warning(paste0("Incorrect select removed: ",
                                   paste(select[inc], collapse = ", ")))
      select <- select[!inc]
    }
  }
  # Extract terms
  mt <- stats::terms(object)
  if (!missing(newdata)){
    mt <- stats::delete.response(mt)
    gind <- attr(mt, "specials")$g
    mfind <- stats::model.frame(mt[gind], data = newdata,
                                na.action = na.action)
    alphas <- object$alpha
    newindex <- do.call(cbind, Map("%*%", mfind, alphas))
    colnames(newindex) <- names(alphas)
    if (type == "indices"){
      newindex <- newindex[,select, drop = FALSE]
      return(newindex)
    }
    Xterms <- c(data.frame(newindex), newdata[names(object$sm_mod$Xcov)])
    sind <- attr(mt, "specials")$s
    objx <- cbind(object$indexfit, object$sm_mod$Xcov)
    gterms <- matrix(0, nrow(mfind), p)
    for (j in 1:p){
      if (is.numeric(Xterms[[j]])){
        # Interpolation
        inter_fun <- ifelse(j %in% c(gind, sind), "spline", "approx")
        if(inter_fun == "approx"){
          ## Fix numerical precision issues if any
          # In-sample data range
          objx_range <- range(objx[,j])
          # Tolerance for absolute difference between two numbers
          tolerance <- sqrt(.Machine$double.eps)
          # Difference with minimum
          Xterm_diff_min <- abs(Xterms[[j]] - objx_range[1])
          # Less than minimum
          smaller_indx <- which(Xterm_diff_min < tolerance)
          # Set to minimum
          Xterms[[j]][smaller_indx] <- objx_range[1]
          # Difference with maximum
          Xterm_diff_max <- abs(Xterms[[j]] - objx_range[2])
          # Less than minimum
          larger_indx <- which(Xterm_diff_max < tolerance)
          # Set to maximum
          Xterms[[j]][larger_indx] <- objx_range[2]
        }
        gterms[,j] <- suppressWarnings(do.call(inter_fun,
                                               list(x = objx[,j],
                                                    y = object$gfit[,j],
                                                    xout = Xterms[[j]]))$y)
      } else {
        faccoefs <- unique(object$gfit[,j])
        names(faccoefs) <- unique(objx[,j])
        gterms[,j] <- faccoefs[Xterms[[j]]]
      }
    }
    colnames(gterms) <- colnames(object$gfit)
    if (type == "terms"){
      gterms <- gterms[,select, drop = FALSE]
      return(gterms)
    }
    betas <- object$beta
    if(type == "scterms"){
      scgterms <- gterms * matrix(betas[colnames(gterms)], nrow(gterms),
                                  ncol(gterms), byrow = TRUE)
      return(scgterms[, select, drop = FALSE])
    }
    yhat <- cbind(1, gterms) %*% betas
    return(yhat)
  } else {
    out <- switch(type,
                  response = object$fitted,
                  terms = object$gfit[, select, drop = FALSE],
                  scterms = (object$gfit * matrix(object$beta[colnames(object$gfit)],
                                                  n, p,
                                                  byrow = TRUE))[, select, drop = FALSE],
                  indices = object$indexfit[, select, drop = FALSE]
    )
    return(out)
  }
}