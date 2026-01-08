# Augment function for class `lmFit`

Generates residuals and fitted values of a fitted `lmFit` object.

## Usage

``` r
# S3 method for class 'lmFit'
augment(x, ...)
```

## Arguments

- x:

  A `lmFit` object.

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

# Predictor variables
linear.vars <- colnames(sim_data)[3:8]

# Model fitting
lmModel <- model_lm(data = sim_data,
                    yvar = "y",
                    linear.vars = linear.vars)
#> [1] "model 1"
# Obtain residuals and fitted values
augment(lmModel)
#> # A tibble: 1,000 × 3
#>    Index  .resid .fitted
#>    <int>   <dbl>   <dbl>
#>  1     6 -0.128    0.695
#>  2     7 -0.292    1.06 
#>  3     8  0.613    3.00 
#>  4     9 -0.233    1.48 
#>  5    10 -0.349    1.26 
#>  6    11  0.584    3.00 
#>  7    12 -0.187    2.04 
#>  8    13 -0.318    1.62 
#>  9    14  0.0687   2.42 
#> 10    15 -0.0921   0.242
#> # ℹ 990 more rows
```
