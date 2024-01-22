#' Obtaining forecasts on a test set from a fitted `smimodel`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `smimodel` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tibble`).
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
  if (!is_tibble(newdata)) stop("newdata is not a tibble.")
  predict_fn <- mgcv::predict.gam
  list_index <- object[1:(length(object)-4)]
  alpha <- vector(mode = "list", length = length(list_index))
  for(i in 1:length(list_index)){
    alpha[[i]] <- list_index[[i]]$coefficients
  }
  alpha <- unlist(alpha)
  if(all(alpha == 0)){
    if(recursive == TRUE){
      predictions =  vector(mode = "list", length = NROW(newdata))
      #pb <- lazybar::lazyProgressBar(NROW(newdata) - 1)
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
        #pb$tick()$print()
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
      #pb <- lazybar::lazyProgressBar(NROW(newdata) - 1)
      for(m in 1:(NROW(newdata) - 1)){
        data_temp = newdata[m, ]
        X_test <- as.matrix(newdata[m, object$vars_index])
        # if("intercept" %in% names(list_index[[1]])){
        #   X_test <- X_test %>%
        #     tibble::as_tibble() %>%
        #     dplyr::mutate(intercept = 1)
        #   # Calculating indices
        #   ind <- vector(mode = "list", length = length(list_index))
        #   for(i in 1:length(list_index)){
        #     ind[[i]] <- as.numeric(as.matrix(X_test[ , object$vars_index]) %*%
        #                              as.matrix(list_index[[i]]$coefficients, ncol = 1)) +
        #       as.numeric(as.matrix(X_test[, "intercept"]) %*%
        #                    as.matrix(list_index[[i]]$intercept, ncol = 1))
        #   }
        # }else{
        #   # Calculating indices
        #   ind <- vector(mode = "list", length = length(list_index))
        #   for(i in 1:length(list_index)){
        #     ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
        #   }
        # }
        
        # Calculating indices
        ind <- vector(mode = "list", length = length(list_index))
        for(i in 1:length(list_index)){
          ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
        }
        names(ind) <- names(list_index)
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
        #pb$tick()$print()
      }
      data_temp = newdata[NROW(newdata), ]
      X_test <- as.matrix(newdata[NROW(newdata), object$vars_index])
      # if("intercept" %in% names(list_index[[1]])){
      #   X_test <- X_test %>%
      #     tibble::as_tibble() %>%
      #     dplyr::mutate(intercept = 1)
      #   # Calculating indices
      #   ind <- vector(mode = "list", length = length(list_index))
      #   for(i in 1:length(list_index)){
      #     ind[[i]] <- as.numeric(as.matrix(X_test[ , object$vars_index]) %*%
      #                              as.matrix(list_index[[i]]$coefficients, ncol = 1)) +
      #       as.numeric(as.matrix(X_test[, "intercept"]) %*%
      #                    as.matrix(list_index[[i]]$intercept, ncol = 1))
      #   }
      # }else{
      #   # Calculating indices
      #   ind <- vector(mode = "list", length = length(list_index))
      #   for(i in 1:length(list_index)){
      #     ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
      #   }
      # }
      
      # Calculating indices
      ind <- vector(mode = "list", length = length(list_index))
      for(i in 1:length(list_index)){
        ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
      }
      names(ind) <- names(list_index)
      dat <- tibble::as_tibble(ind)
      data_list[[NROW(newdata)]] <- dplyr::bind_cols(data_temp, dat)
      predictions[[NROW(newdata)]] = predict_fn(object$gam, data_list[[NROW(newdata)]], type = "response")
      newdata1 <- dplyr::bind_rows(data_list)
      pred <- unlist(predictions)
      pred_F <- newdata1 %>% 
        dplyr::mutate(.predict = pred) 
    }else if(recursive == FALSE){
      X_test <- as.matrix(newdata[ , object$vars_index])
      # if("intercept" %in% names(list_index[[1]])){
      #   X_test <- X_test %>%
      #     tibble::as_tibble() %>%
      #     dplyr::mutate(intercept = 1)
      #   # Calculating indices
      #   ind <- vector(mode = "list", length = length(list_index))
      #   for(i in 1:length(list_index)){
      #     ind[[i]] <- as.numeric(as.matrix(X_test[ , object$vars_index]) %*%
      #                              as.matrix(list_index[[i]]$coefficients, ncol = 1)) +
      #       as.numeric(as.matrix(X_test[, "intercept"]) %*%
      #                    as.matrix(list_index[[i]]$intercept, ncol = 1))
      #   }
      # }else{
      #   # Calculating indices
      #   ind <- vector(mode = "list", length = length(list_index))
      #   for(i in 1:length(list_index)){
      #     ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
      #   }
      # }
      
      # Calculating indices
      ind <- vector(mode = "list", length = length(list_index))
      for(i in 1:length(list_index)){
        ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
      }
      names(ind) <- names(list_index)
      dat <- tibble::as_tibble(ind)
      data_list <- dplyr::bind_cols(newdata, dat)
      pred <- predict_fn(object$gam, data_list, type = "response")
      pred_F <- data_list %>% 
        dplyr::mutate(.predict = pred)
    }
  }
  return(pred_F)
}
