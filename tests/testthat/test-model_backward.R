# Testing model_backward.R

test_that("tests for model_backward()", {
  n = 1205
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
  # Training set
  sim_train <- sim_data[1:1000, ]
  # Validation set
  sim_val <- sim_data[1001:1200, ]
  # Smooth variables
  s.vars <- colnames(sim_data)[3:8]
  output1 <- model_backward(data = sim_train,
                            val.data = sim_val,
                            yvar = "y1",
                            s.vars = s.vars)
  print(output1)
  
  expect_s3_class(output1, "backward")
  expect_true(!is.null(output1$fit[[1]]))
  expect_s3_class(output1$fit[[1]], "gam")
  expect_s3_class(output1$fit[[1]]$model, "tbl_ts")
})