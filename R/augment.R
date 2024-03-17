#' Augment function for class `smimodel`
#'
#' Generates residuals and fitted values of a fitted `smimodel` object.
#'
#' @param x A `smimodel` object.
#' @param ... Other arguments not currently used.
#'
#'
#' @method augment smimodel
#'
#' @export
augment.smimodel <- function(x, ...) {
  # if (!tsibble::is_tsibble(data)) stop("data is not a tsibble.")
  # yvar <- x$fit[[1]]$best$var_y
  # data1 <- data
  # data_index <- index(data1)
  # data_key <- key(data1)
  # if (length(key(data1)) == 0) {
  #   data1 <- data1 %>%
  #     dplyr::mutate(dummy_key = rep(1, NROW(data1))) %>%
  #     tsibble::as_tsibble(index = data_index, key = dummy_key)
  #   data_key <- key(data1)
  # }
  # key11 <- key(data1)[[1]]
  # key_unique <- unique(as.character(sort(dplyr::pull((data1[, {{ key11 }}])[, 1]))))
  # key_num <- seq_along(key_unique)
  # ref <- data.frame(key_unique, key_num)
  # data1 <- data1 %>%
  #   dplyr::mutate(
  #     num_key = as.numeric(factor(as.character({{ key11 }}), levels = key_unique))
  #   )
  model.resid <- vector(mode = "list", length = NROW(x))
  model.fitted <- vector(mode = "list", length = NROW(x))
  time_variable <- vector(mode = "list", length = NROW(x))
  old_group <- vector(mode = "list", length = NROW(x))
  new_group <- vector(mode = "list", length = NROW(x))
  df <- vector(mode = "list", length = NROW(x))
  for (i in seq(NROW(x))) {
    # df_cat <- data1 %>%
    #   dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
    #                   (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
    #                   (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) %>%
    #   tibble::as_tibble() %>%
    #   dplyr::arrange({{data_index}})
    model.resid[[i]] <- x$fit[[i]]$best$gam$residuals
    model.fitted[[i]] <- x$fit[[i]]$best$gam$fitted.values
    time_variable[[i]] <- x$fit[[i]]$best$gam$model[ , index(x$fit[[i]]$best$gam$model)][[1]]
    #time_variable[[i]] <- df_cat[, data_index][[1]]
    old_group[[i]] <- x$fit[[i]]$best$gam$model[ , key(x$fit[[i]]$best$gam$model)[[1]]][[1]]
    #old_group[[i]] <- df_cat[, key11][[1]]
    new_group[[i]] <- rep(x[[1]][i], nrow(x$fit[[1]]$best$gam$model))
    #new_group[[i]] <- rep(x[[1]][i], nrow(df_cat))
    df[[i]] <- tibble::as_tibble(data.frame(Index = time_variable[[i]], 
                                            New = new_group[[i]],
                                            Old = old_group[[i]],
                                            .resid = model.resid[[i]],
                                            .fitted = model.fitted[[i]])) %>%
      dplyr::filter(New == Old)
  }
  mod_res <- dplyr::bind_rows(df)
  return(mod_res)
}
utils::globalVariables(c("dummy_key", "num_key", "New", "Old"))
#' @export
generics::augment



#' Augment function for class `smimodelFit`
#'
#' Generates residuals and fitted values of a fitted `smimodelFit` object.
#'
#' @param x A `smimodelFit` object.
#' @param ... Other arguments not currently used.
#'
#' @importFrom generics augment
#'
#' @method augment smimodelFit
#'
#' @export
augment.smimodelFit <- function(x, ...) {
  smimodel.resid <- x$gam$residuals
  smimodel.fitted <- x$gam$fitted.values
  df <- tibble::tibble(
    .resid = smimodel.resid,
    .fitted = smimodel.fitted
  )
  return(df)
}



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



#' Augment function for class `pprFit`
#'
#' Generates residuals and fitted values of a fitted `pprFit` object.
#'
#' @param x A `pprFit` object.
#' @param ... Other arguments not currently used.
#'
#' @method augment pprFit
#'
#' @export
augment.pprFit <- function(x, ...) {
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



#' Augment function for class `gaimFit`
#'
#' Generates residuals and fitted values of a fitted `gaimFit` object.
#'
#' @param x A `gaimFit` object.
#' @param ... Other arguments not currently used.
#'
#' @method augment gaimFit
#'
#' @export
augment.gaimFit <- function(x, ...) {
  model.resid <- vector(mode = "list", length = nrow(x))
  model.fitted <- vector(mode = "list", length = nrow(x))
  time_variable <- vector(mode = "list", length = nrow(x))
  old_group <- vector(mode = "list", length = nrow(x))
  new_group <- vector(mode = "list", length = nrow(x))
  df <- vector(mode = "list", length = nrow(x))
  for (i in seq(NROW(x))) {
    model.resid[[i]] <- x$fit[[i]]$residuals
    model.fitted[[i]] <- x$fit[[i]]$fitted
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