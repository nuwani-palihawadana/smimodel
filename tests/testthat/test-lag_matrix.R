# Testing lag_matrix.R

test_that("tests for lag_matrix()", {
  set.seed(123)
  sim_data <- tibble(x_lag_000 = runif(100)) %>%
    mutate(x_lag = lag_matrix(x_lag_000, 3)) %>%
    tidyr::unpack(x_lag, names_sep = "_")

  expect_equal(sim_data$x_lag_001[2:NROW(sim_data)], sim_data$x_lag_000[1:(NROW(sim_data) - 1)])
})