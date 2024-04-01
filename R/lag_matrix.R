#' Function for adding lags of time series variables
#'
#' Generates specified number of lagged variables of the given variable in
#' the form of a tibble.
#'
#' @param variable Variable to be lagged.
#' @param n Number of lags. The default value is `n = 10`.
#' @return A `tibble`.
#'
#' @importFrom purrr set_names map
#' @importFrom dplyr lag
#'
#' @examples
#' # Adding lagged variables to an existing tibble
#' library(dplyr)
#' demand_data <- tsibbledata::vic_elec %>%
#'   mutate(Temperature_lag = lag_matrix(Temperature, 30)) %>%
#'   tidyr::unpack(Temperature_lag, names_sep = "_")
#' @export
lag_matrix <- function(variable, n = 10) {
  indices <- seq_len(n)
  as_tibble(set_names(
    map(indices, ~ lag(variable, .)),
    sprintf("%03d", indices)
  ))
}