# Augment function for class `backward`

Generates residuals and fitted values of a fitted `backward` object.

## Usage

``` r
# S3 method for class 'backward'
augment(x, ...)
```

## Arguments

- x:

  A `backward` object.

- ...:

  Other arguments not currently used.

## Value

A `tibble`.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(tibble)
library(tidyr)
library(tsibble)
#> 
#> Attaching package: ‘tsibble’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, union

# Simulate data
n = 1205
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

# Training set
sim_train <- sim_data[1:1000, ]
# Validation set
sim_val <- sim_data[1001:1200, ]

# Predictors taken as non-linear variables
s.vars <- colnames(sim_data)[3:8]

# Model fitting
backwardModel <- model_backward(data = sim_train,
                                val.data = sim_val,
                                yvar = "y",
                                s.vars = s.vars)
#> [1] "Model 1 fitted!"
# Obtain residuals and fitted values
augment(backwardModel)
#> # A tibble: 1,200 × 3
#>    Index  .resid .fitted
#>    <int>   <dbl>   <dbl>
#>  1     6 -0.415    0.996
#>  2     7 -0.176    0.990
#>  3     8  0.700    3.14 
#>  4     9 -0.242    1.41 
#>  5    10  0.0473   1.01 
#>  6    11  0.213    3.24 
#>  7    12 -0.113    1.90 
#>  8    13 -0.102    1.40 
#>  9    14  0.151    2.26 
#> 10    15 -0.0155   0.300
#> # ℹ 1,190 more rows
```
