# Testing model_gaim.R

test_that("tests for model_gaim()", {
  n = 1005
  set.seed(123)
  sim_data <- tibble(x_lag_000 = runif(n)) %>%
    mutate(
      # Add x_lags
      x_lag = lag_matrix(x_lag_000, 5)) %>%
    tidyr::unpack(x_lag, names_sep = "_") %>%
    mutate(
      # Response variable
      y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
      # Add an index to the data set
      inddd = seq(1, n)) %>%
    tidyr::drop_na() %>%
    select(inddd, y, starts_with("x_lag")) %>%
    # Make the data set a `tsibble`
    tsibble::as_tsibble(index = inddd)
  # Index variables
  index.vars <- colnames(sim_data)[3:7]
  # Assign group indices for each predictor
  index.ind = c(rep(1, 3), rep(2, 2))
  # Smooth variables not entering indices
  s.vars = "x_lag_005"
  output1 <- sim_data %>%
    model_gaim(
      yvar = "y",
      index.vars = index.vars,
      index.ind = index.ind,
      s.vars = s.vars
    )
  
  print(output1)
  
  expect_s3_class(output1, "gaimFit")
  expect_true(!is.null(output1$fit[[1]]))
  expect_s3_class(output1$fit[[1]], "cgaim")
  expect_s3_class(output1$fit[[1]]$model, "tbl_ts")
})