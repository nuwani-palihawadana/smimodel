# Plot estimated smooths from a fitted `smimodel`

Plots the graphs of fitted spline(s). If a set of multiple models are
fitted, plots graphs of fitted spline(s) of a specified model (in
argument `model`) out of the set of multiple models fitted.

## Usage

``` r
# S3 method for class 'smimodel'
autoplot(object, model = 1, ...)
```

## Arguments

- object:

  A `smimodel` object.

- model:

  An `integer` to indicate the smooths of which model (out of the set of
  multiple models fitted) to be plotted.

- ...:

  Other arguments not currently used.

## Value

Plot(s) of fitted spline(s).

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

autoplot(smimodel_ppr)
} # }
```
