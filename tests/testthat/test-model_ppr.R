# Testing model_ppr.R

test_that("tests for model_ppr()", {
  n = 1005
  set.seed(123)
  sim_data <- tibble(x_lag_000 = runif(n)) %>%
    mutate(
      # Add x_lags
      x_lag = lag_matrix(x_lag_000, 5)) %>%
    tidyr::unpack(x_lag, names_sep = "_") %>%
    mutate(
      # Response variable
      y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
      # Add an index to the data set
      inddd = seq(1, n)) %>%
    tidyr::drop_na() %>%
    select(inddd, y1, starts_with("x_lag")) %>%
    # Make the data set a `tsibble`
    tsibble::as_tsibble(index = inddd)
  # Index variables
  index.vars <- colnames(sim_data)[3:8]
  output1 <- sim_data %>%
    model_ppr(
      yvar = "y1",
      index.vars = index.vars
    )
  
  print(output1)
  
  expect_s3_class(output1, "pprFit")
  expect_true(!is.null(output1$fit[[1]]))
  expect_s3_class(output1$fit[[1]], c("ppr.form", "ppr"))
  expect_s3_class(output1$fit[[1]]$model, "tbl_ts")
})