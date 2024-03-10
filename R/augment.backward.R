#' Augment function for class `backward`
#'
#' Generates residuals and fitted values of a fitted `backward` object.
#'
#' @param x A `backward` object.
#' @param ... Other arguments not currently used.
#'
#' @method augment backward
#'
#' @export
augment.backward <- function(x, ...) {
  model.resid <- vector(mode = "list", length = nrow(x))
  model.fitted <- vector(mode = "list", length = nrow(x))
  time_variable <- vector(mode = "list", length = nrow(x))
  old_group <- vector(mode = "list", length = nrow(x))
  new_group <- vector(mode = "list", length = nrow(x))
  df <- vector(mode = "list", length = nrow(x))
  for (i in seq(NROW(x))) {
    model.resid[[i]] <- x$fit[[i]]$residuals
    model.fitted[[i]] <- x$fit[[i]]$fitted.values
    time_variable[[i]] <- x$fit[[i]]$model[, index(x$fit[[i]]$model)][[1]]
    old_group[[i]] <- x$fit[[i]]$model[, key(x$fit[[i]]$model)[[1]]][[1]]
    new_group[[i]] <- rep(x[[1]][i], nrow(x$fit[[i]]$model))
    df[[i]] <- as_tibble(data.frame(
      Index = time_variable[[i]], New = new_group[[i]],
      Old = old_group[[i]], .resid = model.resid[[i]],
      .fitted = model.fitted[[i]]
    )) %>%
      dplyr::filter(New == Old)
  }
  mod_res <- bind_rows(df)
  return(mod_res)
}