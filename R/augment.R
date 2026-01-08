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
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' smimodel_ppr <- model_smimodel(data = sim_data,
#'                                yvar = "y",
#'                                index.vars = index.vars,
#'                                initialise = "ppr")
#'
#' # Obtain residuals and fitted values
#' augment(smimodel_ppr)
#' }
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
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1205
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Training set
#' sim_train <- sim_data[1:1000, ]
#' # Validation set
#' sim_val <- sim_data[1001:1200, ]
#'
#' # Predictors taken as non-linear variables
#' s.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' backwardModel <- model_backward(data = sim_train,
#'                                 val.data = sim_val,
#'                                 yvar = "y",
#'                                 s.vars = s.vars)
#' # Obtain residuals and fitted values
#' augment(backwardModel)
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
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' pprModel <- model_ppr(data = sim_data,
#'                       yvar = "y",
#'                       index.vars = index.vars)
#'
#' # Obtain residuals and fitted values
#' augment(pprModel)
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
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Predictors taken as index variables
#' index.vars <- colnames(sim_data)[3:7]
#'
#' # Assign group indices for each predictor
#' index.ind = c(rep(1, 3), rep(2, 2))
#'
#' # Predictors taken as non-linear variables not entering indices
#' s.vars = "x_lag_005"
#'
#' # Model fitting
#' gaimModel <- model_gaim(data = sim_data,
#'                         yvar = "y",
#'                         index.vars = index.vars,
#'                         index.ind = index.ind,
#'                         s.vars = s.vars)
#' # Obtain residuals and fitted values
#' augment(gaimModel)
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
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Predictor variables
#' linear.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' lmModel <- model_lm(data = sim_data,
#'                     yvar = "y",
#'                     linear.vars = linear.vars)
#' # Obtain residuals and fitted values
#' augment(lmModel)
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
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Predictors taken as non-linear variables
#' s.vars <- colnames(sim_data)[3:6]
#'
#' # Predictors taken as linear variables
#' linear.vars <- colnames(sim_data)[7:8]
#'
#' # Model fitting
#' gamModel <- model_gam(data = sim_data,
#'                       yvar = "y",
#'                       s.vars = s.vars,
#'                       linear.vars = linear.vars)
#'
#' # Obtain residuals and fitted values
#' augment(gamModel)
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