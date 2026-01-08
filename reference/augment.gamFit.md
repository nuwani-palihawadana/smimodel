# Augment function for class `gamFit`

Generates residuals and fitted values of a fitted `gamFit` object.

## Usage

``` r
# S3 method for class 'gamFit'
augment(x, ...)
```

## Arguments

- x:

  A `gamFit` object.

- ...:

  Other arguments not currently used.

## Value

A `tibble`.

## Examples

``` r
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1005
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(n)) |>
  mutate(
    # Add x_lags
    x_lag = lag_matrix(x_lag_000, 5)) |>
  unpack(x_lag, names_sep = "_") |>
  mutate(
    # Response variable
    y = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
    # Add an index to the data set
    inddd = seq(1, n)) |>
  drop_na() |>
  select(inddd, y, starts_with("x_lag")) |>
  # Make the data set a `tsibble`
  as_tsibble(index = inddd)

# Predictors taken as non-linear variables
s.vars <- colnames(sim_data)[3:6]

# Predictors taken as linear variables
linear.vars <- colnames(sim_data)[7:8]

# Model fitting
gamModel <- model_gam(data = sim_data,
                      yvar = "y",
                      s.vars = s.vars,
                      linear.vars = linear.vars)
#> [1] "model 1"

# Obtain residuals and fitted values
augment(gamModel)
#> # A tibble: 1,000 × 3
#>    Index  .resid .fitted
#>    <int>   <dbl>   <dbl>
#>  1     6 -0.392    0.959
#>  2     7 -0.219    0.985
#>  3     8  0.515    3.10 
#>  4     9 -0.166    1.41 
#>  5    10 -0.0998   1.01 
#>  6    11  0.356    3.23 
#>  7    12 -0.0723   1.92 
#>  8    13 -0.0837   1.38 
#>  9    14  0.238    2.25 
#> 10    15 -0.134    0.283
#> # ℹ 990 more rows
```
