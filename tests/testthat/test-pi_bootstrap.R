# Testing pi_bootstrap.R

test_that("tests for pi_bootstrap()", {
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
  # Test set
  sim_test <- sim_data[1001:1200, ]
  # Index variables
  index.vars <- colnames(sim_data)[3:8]
  # Model fitting
  model1 <- model_smimodel(data = sim_train,
                           yvar = "y1",
                           index.vars = index.vars,
                           initialise = "ppr")
  # Calculating lower and upper bounds for 95% prediction interval
  num_futures <- 1000
  PI_model1 <- pi_bootstrap(object = model1,
                           newdata = sim_test,
                           num.futures = num_futures)
  
  expect_equal(names(PI_model1), c("bounds", "sample_paths"))
  expect_true(!is.null(PI_model1$bounds))
  expect_true(!is.null(PI_model1$sample_paths))
  expect_equal(dim(PI_model1$bounds)[1], NROW(sim_test))
  expect_equal(c(dim(PI_model1$sample_paths)[1], dim(PI_model1$sample_paths)[2]), 
               c(NROW(sim_test), num_futures))
})