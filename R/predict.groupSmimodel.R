#' Obtaining forecasts on a test set from a fitted `groupSmimodel`
#'
#' Gives forecasts on a test set.
#'
#' @param object A `groupSmimodel` object.
#' @param newdata The set of new data on for which the forecasts are required
#'   (i.e. test set; should be a `tsibble`).
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class `tsibble`.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an integer. If `neighbour = x`, `x`
#'   number of keys before the key of interest and `x` number of keys after the
#'   key of interest are grouped together for model fitting. The default is `0`
#'   (i.e. no neighbours are considered for model fitting).
#' @param recursive Whether to obtain recursive forecasts or not (default -
#'   FALSE).
#' @param recursive_colRange If `recursive = TRUE`, The range of column numbers
#'   in `newdata` to be filled with forecasts.
#' @param ... Other arguments not currently used.
#'
#' @method predict groupSmimodel
#'
#' @export
predict.groupSmimodel <- function(object, newdata, data, neighbour = 0,
                                    recursive = FALSE, recursive_colRange = NULL, ...) {
  if (!tsibble::is_tsibble(newdata)) stop("newdata is not a tsibble.")
  if (!tsibble::is_tsibble(data)) stop("data is not a tsibble.")
  data1 <- data
  data_index <- index(data1)
  data_key <- key(data1)
  if (length(key(data1)) == 0) {
    data1 <- data1 %>%
      dplyr::mutate(dummy_key = rep(1, NROW(data1))) %>%
      tsibble::as_tsibble(index = data_index, key = dummy_key)
    data_key <- key(data1)
  }
  key11 <- key(data1)[[1]]
  key_unique <- unique(as.character(sort(dplyr::pull((data1[, {{ key11 }}])[, 1]))))
  key_num <- seq_along(key_unique)
  ref <- data.frame(key_unique, key_num)
  data1 <- data1 %>%
    dplyr::mutate(
      num_key = as.numeric(factor(as.character({{ key11 }}), levels = key_unique))
    )
  index_n <- index(newdata)
  key_n <- key(newdata)
  if (length(key(newdata)) == 0) {
    newdata <- newdata %>%
      dplyr::mutate(dummy_key = rep(1, NROW(newdata))) %>%
      tsibble::as_tsibble(index = index_n, key = dummy_key)
    key_n <- key(newdata)
  }
  predict_fn <- mgcv::predict.gam
  yvar <- object$fit[[1]]$var_y
  gam <- vector(mode = "list", length = NROW(object))
  for(j in 1:NROW(object)){
    df_cat <- data1 %>%
      dplyr::filter((abs(num_key - ref$key_num[j]) <= neighbour) |
                      (abs(num_key - ref$key_num[j] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[j] - NROW(ref)) <= neighbour)) %>%
      tibble::as_tibble() %>%
      dplyr::arrange({{data_index}})
    gam[[j]] <- make_gam(x = object$fit[[j]], data = df_cat)
  }
  gam_list <- list(object$key, gam)
  ref_gam <- tibble::as_tibble(
    x = gam_list, .rows = length(gam_list[[1]]),
    .name_repair = ~ vctrs::vec_as_names(..., repair = "universal", quiet = TRUE)
  )
  ref_gam <- ref_gam %>%
    dplyr::rename(key = ...1) %>%
    dplyr::rename(gam = ...2)
  if(recursive == TRUE){
    newdata <- newdata %>%
      tibble::as_tibble() %>%
      dplyr::arrange({{index_n}})
    predictions =  vector(mode = "list", length = NROW(newdata))
    data_list <- vector(mode = "list", length = NROW(newdata))
    #pb <- lazybar::lazyProgressBar(NROW(newdata) - 1)
    for(m in 1:(NROW(newdata) - 1)){
      data_temp = newdata[m, ]
      key22 = data_temp[ , {{ key11 }}][[1]]
      key22_pos = which(object$key == key22)
      list_index <- object$fit[[key22_pos]][1:(length(object$fit[[key22_pos]])-3)]
      alpha <- vector(mode = "list", length = length(list_index))
      for(b in 1:length(list_index)){
        alpha[[b]] <- list_index[[b]]$coefficients
      }
      alpha <- unlist(alpha)
      if(all(alpha == 0)){
        data_list[[m]] <- data_temp
      }else{
        X_test <- as.matrix(newdata[m, object$fit[[key22_pos]]$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = length(list_index))
        for(i in 1:length(list_index)){
          ind[[i]] <- as.numeric(X_test %*% as.matrix(list_index[[i]]$coefficients, ncol = 1))
        }
        names(ind) <- names(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[m]] <- dplyr::bind_cols(data_temp, dat)
      }
      pred <- predict_fn(ref_gam$gam[[key22_pos]], data_list[[m]], type = "response")
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
    key22 = data_temp[ , {{ key11 }}][[1]]
    key22_pos = which(object$key == key22)
    list_index <- object$fit[[key22_pos]][1:(length(object$fit[[key22_pos]])-3)]
    alpha <- vector(mode = "list", length = length(list_index))
    for(w in 1:length(list_index)){
      alpha[[w]] <- list_index[[w]]$coefficients
    }
    alpha <- unlist(alpha)
    if(all(alpha == 0)){
      data_list[[NROW(newdata)]] <- data_temp
    }else{
      X_test <- as.matrix(newdata[NROW(newdata), object$fit[[key22_pos]]$vars_index])
      # Calculating indices
      ind <- vector(mode = "list", length = length(list_index))
      for(a in 1:length(list_index)){
        ind[[a]] <- as.numeric(X_test %*% as.matrix(list_index[[a]]$coefficients, ncol = 1))
      }
      names(ind) <- names(list_index)
      dat <- tibble::as_tibble(ind)
      data_list[[NROW(newdata)]] <- dplyr::bind_cols(data_temp, dat)
    }
    pred <- predict_fn(ref_gam$gam[[key22_pos]], data_list[[NROW(newdata)]], type = "response")
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
      list_index <- object$fit[[i]][1:(length(object$fit[[i]])-3)]
      alpha <- vector(mode = "list", length = length(list_index))
      for(z in 1:length(list_index)){
        alpha[[z]] <- list_index[[z]]$coefficients
      }
      alpha <- unlist(alpha)
      if(all(alpha == 0)){
        data_list[[i]] <- newdata_cat
      }else{
        X_test <- as.matrix(newdata_cat[ , object$fit[[i]]$vars_index])
        # Calculating indices
        ind <- vector(mode = "list", length = length(list_index))
        for(k in 1:length(ind)){
          ind[[k]] <- as.numeric(X_test %*% as.matrix(list_index[[k]]$coefficients, ncol = 1))
        }
        names(ind) <- names(list_index)
        dat <- tibble::as_tibble(ind)
        data_list[[i]] <- dplyr::bind_cols(newdata_cat, dat) 
      }
      key_pos <- which(ref_gam$key == object$key[i])
      predictions[[i]] <- predict_fn(ref_gam$gam[[key_pos]], data_list[[i]], 
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