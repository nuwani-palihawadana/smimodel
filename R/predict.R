#' Obtaining forecasts on a test set from a fitted `smimodel`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `smimodel` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `newdata` to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @method predict smimodel
#'
#' @export
predict.smimodel <- function(object, newdata, recursive = FALSE, 
                             recursive_colRange = NULL, ...) {
  if (!tsibble::is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata %>%
      dplyr::mutate(dummy_key = rep(1, NROW(newdata))) %>%
      tsibble::as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  predict_fn <- mgcv::predict.gam
  if(recursive == TRUE){
    newdata <- newdata %>%
      tibble::as_tibble() %>%
      dplyr::arrange({{index_n}})
    predictions =  vector(mode = "list", length = NROW(newdata))
    data_list <- vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      list_index <- object$fit[[key22_pos]]$best$alpha
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
        X_test <- as.matrix(newdata[m, object$fit[[key22_pos]]$best$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = numInd)
        for(i in 1:numInd){
          ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
        }
        names(ind) <- colnames(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[m]] <- dplyr::bind_cols(data_temp, dat)
      }
      pred <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list[[m]], type = "response")
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
    list_index <- object$fit[[key22_pos]]$best$alpha
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
      X_test <- as.matrix(newdata[NROW(newdata), object$fit[[key22_pos]]$best$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = numInd)
      for(a in 1:numInd){
        ind[[a]] <- as.numeric(X_test %*% as.matrix(list_index[ , a], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list[[NROW(newdata)]] <- dplyr::bind_cols(data_temp, dat)
    }
    pred <- predict_fn(object$fit[[key22_pos]]$best$gam, data_list[[NROW(newdata)]], type = "response")
    predictions[[NROW(newdata)]] <- pred
    newdata1 <- dplyr::bind_rows(data_list)
    pred <- unlist(predictions)
    pred_F <- newdata1 %>% 
      dplyr::mutate(.predict = pred) %>%
      tsibble::as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    data_list <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key11 }}] == object$key[i], ]
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
      predictions[[i]] <- predict_fn(object$fit[[i]]$best$gam, data_list[[i]], 
                                     type = "response")
    }
    newdata1 <- dplyr::bind_rows(data_list)
    pred <- unlist(predictions)
    pred_F <- newdata1 %>% 
      dplyr::mutate(.predict = pred) %>%
      tsibble::as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }
  return(pred_F)
}



#' Obtaining forecasts on a test set from a `smimodelFit`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `smimodelFit` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tibble`).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `newdata` to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @method predict smimodelFit
#'
#' @export
predict.smimodelFit <- function(object, newdata, recursive = FALSE, 
                                recursive_colRange = NULL, ...) {
  if (!is_tibble(newdata)) stop("newdata is not a tibble.")
  predict_fn <- mgcv::predict.gam
  list_index <- object$alpha
  num_ind <- NCOL(list_index)
  alpha <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    alpha[[i]] <- list_index[ , i]
  }
  alpha <- unlist(alpha)
  names(alpha) <- NULL
  if(all(alpha == 0)){
    if(recursive == TRUE){
      predictions =  vector(mode = "list", length = NROW(newdata))
      for(m in 1:(NROW(newdata) - 1)){
        data_temp = newdata[m, ]
        pred <- predict_fn(object$gam, data_temp, type = "response")
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
      predictions[[NROW(newdata)]] = predict_fn(object$gam, data_temp, type = "response")
      pred <- unlist(predictions)
      pred_F <- newdata %>% 
        dplyr::mutate(.predict = pred) 
    }else if(recursive == FALSE){
      pred <- predict_fn(object$gam, newdata, type = "response")
      pred_F <- newdata %>% 
        dplyr::mutate(.predict = pred)
    }
  }else{
    if(recursive == TRUE){
      predictions =  vector(mode = "list", length = NROW(newdata))
      data_list <- vector(mode = "list", length = NROW(newdata))
      for(m in 1:(NROW(newdata) - 1)){
        data_temp = newdata[m, ]
        X_test <- as.matrix(newdata[m, object$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = num_ind)
        for(i in 1:num_ind){
          ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
        }
        names(ind) <- colnames(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[m]] <- dplyr::bind_cols(data_temp, dat)
        pred <- predict_fn(object$gam, data_list[[m]], type = "response")
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
      X_test <- as.matrix(newdata[NROW(newdata), object$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = num_ind)
      for(i in 1:num_ind){
        ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list[[NROW(newdata)]] <- dplyr::bind_cols(data_temp, dat)
      predictions[[NROW(newdata)]] = predict_fn(object$gam, data_list[[NROW(newdata)]], type = "response")
      newdata1 <- dplyr::bind_rows(data_list)
      pred <- unlist(predictions)
      pred_F <- newdata1 %>% 
        dplyr::mutate(.predict = pred) 
    }else if(recursive == FALSE){
      X_test <- as.matrix(newdata[ , object$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = num_ind)
      for(i in 1:num_ind){
        ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[ , i], ncol = 1))
      }
      names(ind) <- colnames(list_index)
      dat <- tibble::as_tibble(ind)
      data_list <- dplyr::bind_cols(newdata, dat)
      pred <- predict_fn(object$gam, data_list, type = "response")
      pred_F <- data_list %>% 
        dplyr::mutate(.predict = pred)
    }
  }
  return(pred_F)
}



#' Obtaining forecasts on a test set from a fitted `backward`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `backward` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `newdata` to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @method predict backward
#'
#' @export
predict.backward <- function(object, newdata, 
                             recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata %>%
      mutate(dummy_key = rep(1, NROW(newdata))) %>%
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  predict_fn <- mgcv::predict.gam
  key11 <- key(newdata)[[1]]
  if(recursive == TRUE){
    newdata <- newdata %>%
      as_tibble() %>%
      arrange({{index_n}})
    predictions =  vector(mode = "list", length = NROW(newdata))
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      pred <- predict_fn(object$fit[[key22_pos]], data_temp, type = "response")
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
    predictions[[NROW(newdata)]] = predict_fn(object$fit[[key22_pos]], data_temp, 
                                              type = "response")
    pred <- unlist(predictions)
    pred_F <- newdata %>% 
      mutate(.predict = pred) %>%
      as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key11 }}] == object$key[i], ]
      predictions[[i]] <- predict_fn(object$fit[[i]], newdata_cat, type = "response")
    }
    pred <- unlist(predictions)
    pred_F <- newdata %>% 
      mutate(.predict = pred) %>%
      as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }
  return(pred_F)
}



#' Obtaining forecasts on a test set from a fitted `pprFit`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `pprFit` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `newdata` to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @method predict pprFit
#'
#' @export
predict.pprFit <- function(object, newdata, 
                           recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata %>%
      mutate(dummy_key = rep(1, NROW(newdata))) %>%
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  if(recursive == TRUE){
    newdata <- newdata %>%
      as_tibble() %>%
      arrange({{index_n}})
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
    pred <- unlist(predictions)
    pred_F <- newdata %>% 
      mutate(.predict = pred) %>%
      as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key11 }}] == object$key[i], ]
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
    pred <- unlist(predictions)
    pred_F <- newdata %>% 
      mutate(.predict = pred) %>%
      as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }
  return(pred_F)
}



#' Obtaining forecasts on a test set from a fitted `gaimFit`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `gaimFit` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `newdata` to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @method predict gaimFit
#'
#' @export
predict.gaimFit <- function(object, newdata, 
                            recursive = FALSE, recursive_colRange = NULL, ...){
  if (!is_tsibble(newdata)) stop("newdata is not a tsibble.")
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata %>%
      mutate(dummy_key = rep(1, NROW(newdata))) %>%
      as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  key11 <- key(newdata)[[1]]
  if(recursive == TRUE){
    newdata <- newdata %>%
      as_tibble() %>%
      arrange({{index_n}})
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
    pred <- unlist(predictions)
    pred_F <- newdata %>% 
      mutate(.predict = pred) %>%
      as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }else if(recursive == FALSE){
    predictions <- vector(mode = "list", length = NROW(object))
    for (i in 1:NROW(object)) {
      newdata_cat <- newdata[newdata[{{ key11 }}] == object$key[i], ]
      predictions[[i]] <- predict(object$fit[[i]], newdata_cat, type = "response")
    }
    pred <- unlist(predictions)
    pred_F <- newdata %>% 
      mutate(.predict = pred) %>%
      as_tsibble(index = {{ index_n }}, key = {{ key11 }})
  }
  return(pred_F)
}