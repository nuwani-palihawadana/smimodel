# Augment function for class `smimodel`

Generates residuals and fitted values of a fitted `smimodel` object.

## Usage

``` r
# S3 method for class 'smimodel'
augment(x, ...)
```

## Arguments

- x:

  A `smimodel` object.

- ...:

  Other arguments not currently used.

## Value

A `tibble`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
library(ROI)
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
smimodel_ppr <- model_smimodel(data = sim_data,
                               yvar = "y",
                               index.vars = index.vars,
                               initialise = "ppr")

# Obtain residuals and fitted values
augment(smimodel_ppr)
} # }
```
