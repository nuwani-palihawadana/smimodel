#' Augment function for class \code{smimodel}
#'
#' Generates residuals and fitted values of a fitted \code{smimodel} object.
#'
#' @param x A \code{smimodel} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
#'
#' @method augment smimodel
#'
#' @export
augment.smimodel <- function(x, ...) {
  model.resid <- vector(mode = "list", length = NROW(x))
  model.fitted <- vector(mode = "list", length = NROW(x))
  time_variable <- vector(mode = "list", length = NROW(x))
  old_group <- vector(mode = "list", length = NROW(x))
  new_group <- vector(mode = "list", length = NROW(x))
  df <- vector(mode = "list", length = NROW(x))
  for (i in seq(NROW(x))) {
    model.resid[[i]] <- x$fit[[i]]$best$gam$residuals
    model.fitted[[i]] <- x$fit[[i]]$best$gam$fitted.values
    time_variable[[i]] <- x$fit[[i]]$best$gam$model[ , index(x$fit[[i]]$best$gam$model)][[1]]
    old_group[[i]] <- x$fit[[i]]$best$gam$model[ , key(x$fit[[i]]$best$gam$model)[[1]]][[1]]
    new_group[[i]] <- rep(x[[1]][i], nrow(x$fit[[1]]$best$gam$model))
    df[[i]] <- tibble::as_tibble(data.frame(Index = time_variable[[i]], 
                                            New = new_group[[i]],
                                            Old = old_group[[i]],
                                            .resid = model.resid[[i]],
                                            .fitted = model.fitted[[i]])) |>
      dplyr::filter(New == Old) |> 
      select(Index, .resid, .fitted)
  }
  mod_res <- dplyr::bind_rows(df)
  return(mod_res)
}
utils::globalVariables(c("Index", ".resid", ".fitted", "New", "Old"))


#' Augment function for class \code{smimodelFit}
#'
#' Generates residuals and fitted values of a fitted \code{smimodelFit} object.
#'
#' @param x A \code{smimodelFit} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
#'
#' @method augment smimodelFit
augment.smimodelFit <- function(x, ...) {
  smimodel.resid <- x$gam$residuals
  smimodel.fitted <- x$gam$fitted.values
  df <- tibble::tibble(
    .resid = smimodel.resid,
    .fitted = smimodel.fitted
  )
  return(df)
}


#' Augment function for class \code{backward}
#'
#' Generates residuals and fitted values of a fitted \code{backward} object.
#'
#' @param x A \code{backward} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
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
    )) |>
      dplyr::filter(New == Old) |> 
      select(Index, .resid, .fitted)
  }
  mod_res <- bind_rows(df)
  return(mod_res)
}


#' Augment function for class \code{pprFit}
#'
#' Generates residuals and fitted values of a fitted \code{pprFit} object.
#'
#' @param x A \code{pprFit} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
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
    )) |>
      dplyr::filter(New == Old) |> 
      select(Index, .resid, .fitted)
  }
  mod_res <- bind_rows(df)
  return(mod_res)
}


#' Augment function for class \code{gaimFit}
#'
#' Generates residuals and fitted values of a fitted \code{gaimFit} object.
#'
#' @param x A \code{gaimFit} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
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
    )) |>
      dplyr::filter(New == Old) |> 
      select(Index, .resid, .fitted)
  }
  mod_res <- bind_rows(df)
  return(mod_res)
}


#' Augment function for class \code{lmFit}
#'
#' Generates residuals and fitted values of a fitted \code{lmFit} object.
#'
#' @param x A \code{lmFit} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
#'
#' @method augment lmFit
#'
#' @export
augment.lmFit <- function(x, ...) {
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
    )) |>
      dplyr::filter(New == Old) |> 
      select(Index, .resid, .fitted)
  }
  mod_res <- bind_rows(df)
  return(mod_res)
}


#' Augment function for class \code{gamFit}
#'
#' Generates residuals and fitted values of a fitted \code{gamFit} object.
#'
#' @param x A \code{gamFit} object.
#' @param ... Other arguments not currently used.
#' @return A \code{tibble}.
#'
#' @method augment gamFit
#'
#' @export
augment.gamFit <- function(x, ...) {
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
    )) |>
      dplyr::filter(New == Old) |> 
      select(Index, .resid, .fitted)
  }
  mod_res <- bind_rows(df)
  return(mod_res)
}