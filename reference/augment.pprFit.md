# Augment function for class `pprFit`

Generates residuals and fitted values of a fitted `pprFit` object.

## Usage

``` r
# S3 method for class 'pprFit'
augment(x, ...)
```

## Arguments

- x:

  A `pprFit` object.

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

# Index variables
index.vars <- colnames(sim_data)[3:8]

# Model fitting
pprModel <- model_ppr(data = sim_data,
                      yvar = "y",
                      index.vars = index.vars)
#> [1] "model 1"

# Obtain residuals and fitted values
augment(pprModel)
#> # A tibble: 1,000 × 3
#>    Index   .resid .fitted
#>    <int>    <dbl>   <dbl>
#>  1     6  0.0637    0.503
#>  2     7  0.0388    0.726
#>  3     8 -0.0444    3.66 
#>  4     9  0.0492    1.20 
#>  5    10 -0.0509    0.964
#>  6    11 -0.0312    3.61 
#>  7    12  0.00121   1.85 
#>  8    13 -0.00829   1.31 
#>  9    14  0.00486   2.48 
#> 10    15 -0.108     0.258
#> # ℹ 990 more rows
```
