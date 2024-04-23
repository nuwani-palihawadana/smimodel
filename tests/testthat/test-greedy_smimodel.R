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
  # Penalty parameter values to search
  # L0 penalty
  lambda0 = seq(1, 12, by = 1)
  # L2 penalty
  lambda2 = seq(0, 12, by = 1)
  # Full grid
  grid1 <- expand.grid(lambda0, lambda2)
  # Starting point options
  starting <- grid1[c(1, 6, 12, 73, 78, 84, 145, 150, 156), ]
  # L0 penalty
  lambda0_start = as.numeric(unique(unlist(starting[1])))
  # L2 penalty
  lambda2_start = as.numeric(unique(unlist(starting[2])))
  # Model fitting
  output1 <- greedy_smimodel(data = sim_train,
                             val.data = sim_val,
                             yvar = "y1",
                             index.vars = index.vars,
                             initialise = "additive",
                             s.vars = s.vars,
                             linear.vars = linear.vars,
                             lambda0_seq = lambda0,
                             lambda2_seq = lambda2,
                             lambda_step = 1,
                             lambda0_start_seq = lambda0_start,
                             lambda2_start_seq = lambda2_start)
  
  print(output1)
  
  expect_equal(names(output1$fit[[1]]), c("initial", "best", "best_lambdas"))
  expect_equal(length(output1$fit[[1]]$best_lambdas), 2)
  expect_equal(names(output1$fit[[1]]$best), c("alpha", "derivatives", 
                                               "var_y", "vars_index", 
                                               "vars_s", "vars_linear", 
                                               "neighbour", "gam"))
  expect_true(!is.null(output1$fit[[1]]$best$alpha))
  expect_true(!is.null(output1$fit[[1]]$best$gam))
  expect_equal(output1$fit[[1]]$best$vars_s, s.vars)
  expect_equal(output1$fit[[1]]$best$vars_linear, linear.vars)
  expect_s3_class(output1$fit[[1]]$best$gam$model, "tbl_ts")
})