# Testing greedy_smimodel.R

test_that("tests for greedy_smimodel()", {
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
    drop_na() %>%
    select(inddd, y1, starts_with("x_lag")) %>%
    # Make the data set a `tsibble`
    as_tsibble(index = inddd)
  # Training set
  sim_train <- sim_data[1:1000, ]
  # Validation set
  sim_val <- sim_data[1001:1200, ]
  # Predictor variables
  index.vars = c("x_lag_000", "x_lag_001", "x_lag_002", "x_lag_003")
  s.vars = "x_lag_004"
  linear.vars = "x_lag_005"
  
  # Model fitting
  output1 <- greedy_smimodel(data = sim_train,
                             val.data = sim_val,
                             yvar = "y1",
                             index.vars = index.vars,
                             initialise = "ppr",
                             s.vars = s.vars,
                             linear.vars = linear.vars)
  
  print(output1)
  
  expect_equal(names(output1$fit[[1]]), c("initial", "best", "best_lambdas", "lambda0_seq", "lambda2_seq", "searched"))
  expect_equal(length(output1$fit[[1]]$best_lambdas), 2)
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