#' Function for adding lags of time series variables
#'
#' Generates specified number of lagged variables of the given variable in the
#' form of a tibble.
#'
#' @param variable Variable to be lagged.
#' @param n Number of lags. The default value is \code{n = 10}.
#' @return A \code{tibble}.
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(tidyr)
#' # Adding lagged variables to an existing tibble
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(100)) |>
#'   mutate(x_lag = lag_matrix(x_lag_000, 3)) |>
#'   unpack(x_lag, names_sep = "_")
#' @export
lag_matrix <- function(variable, n = 10) {
  indices <- seq_len(n)
  as_tibble(set_names(
    map(indices, ~ lag(variable, .)),
    sprintf("%03d", indices)
  ))
}