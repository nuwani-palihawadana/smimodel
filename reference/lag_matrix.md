# Function for adding lags of time series variables

Generates specified number of lagged variables of the given variable in
the form of a tibble.

## Usage

``` r
lag_matrix(variable, n = 10)
```

## Arguments

- variable:

  Variable to be lagged.

- n:

  Number of lags. The default value is `n = 10`.

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
# Adding lagged variables to an existing tibble
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(100)) |>
  mutate(x_lag = lag_matrix(x_lag_000, 3)) |>
  unpack(x_lag, names_sep = "_")
```
