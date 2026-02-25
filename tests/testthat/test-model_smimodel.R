# Testing model_smimodel.R

test_that("tests for model_smimodel()", {
  skip_on_cran()
  # Test code that requires Gurobi
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
  # Predictor variables
  index.vars = c("x_lag_000", "x_lag_001", "x_lag_002", "x_lag_003")
  s.vars = "x_lag_004"
  linear.vars = "x_lag_005"
  output1 <- sim_data %>%
    model_smimodel(
      yvar = "y",
      index.vars = index.vars,
      initialise = "additive",
      s.vars = s.vars,
      linear.vars = linear.vars
    )
  
  print(output1)
  
  expect_equal(names(output1$fit[[1]]), c("initial", "best"))
  expect_equal(names(output1$fit[[1]]$best), c("alpha", "derivatives", 
                                               "var_y", "vars_index", 
                                               "vars_s", "vars_linear",
                                               "vars_range",
                                               "neighbour", "gam",
                                               "lambda0", "lambda2",
                                               "M", "max.iter",
                                               "tol", "tolCoefs", "TimeLimit", 
                                               "MIPGap", "Nonconvex"))
  expect_true(!is.null(output1$fit[[1]]$best$alpha))
  expect_true(!is.null(output1$fit[[1]]$best$gam))
  expect_equal(output1$fit[[1]]$best$vars_s, s.vars)
  expect_equal(output1$fit[[1]]$best$vars_linear, linear.vars)
  expect_s3_class(output1$fit[[1]]$best$gam$model, "tbl_ts")
})