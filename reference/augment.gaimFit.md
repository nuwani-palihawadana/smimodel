# Augment function for class `gaimFit`

Generates residuals and fitted values of a fitted `gaimFit` object.

## Usage

``` r
# S3 method for class 'gaimFit'
augment(x, ...)
```

## Arguments

- x:

  A `gaimFit` object.

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

# Predictors taken as index variables
index.vars <- colnames(sim_data)[3:7]

# Assign group indices for each predictor
index.ind = c(rep(1, 3), rep(2, 2))

# Predictors taken as non-linear variables not entering indices
s.vars = "x_lag_005"

# Model fitting
gaimModel <- model_gaim(data = sim_data,
                        yvar = "y",
                        index.vars = index.vars,
                        index.ind = index.ind,
                        s.vars = s.vars)
#> [1] "model 1"
# Obtain residuals and fitted values
augment(gaimModel)
#> # A tibble: 1,000 × 3
#>    Index  .resid .fitted
#>    <int>   <dbl>   <dbl>
#>  1     6  0.160    0.407
#>  2     7 -0.230    0.995
#>  3     8  0.381    3.24 
#>  4     9 -0.162    1.41 
#>  5    10 -0.0843   0.997
#>  6    11  0.403    3.18 
#>  7    12  0.0542   1.80 
#>  8    13 -0.0269   1.33 
#>  9    14  0.249    2.24 
#> 10    15 -0.0629   0.213
#> # ℹ 990 more rows
```
