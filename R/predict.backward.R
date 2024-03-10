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
    #pb <- lazybar::lazyProgressBar(NROW(newdata) - 1)
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
      #pb$tick()$print()
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