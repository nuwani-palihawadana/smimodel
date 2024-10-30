
<!-- README.md is generated from README.Rmd. Please edit that file -->

# smimodel <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The R package *smimodel* provides functions to estimate Sparse Multiple
Index (SMI) Models for nonparametric forecasting/prediction. The SMI
Modelling algorithm simultaneously performs optimal predictor selection
(hence “sparse”) and predictor grouping, enabling parsimonious
forecasting/prediction models in a high-dimensional context.

The package also includes functions to fit some benchmark comparison
methods namely nonparametric additive model with backward elimination,
group-wise additive index model and projection pursuit regression. These
functions are designed to have a similar syntax and output format for
convenient usage and interchangeability.

For detailed information regarding the SMI Modelling algorithm refer the
related working paper at
<https://www.monash.edu/business/ebs/research/publications/ebs/2024/wp16-2024.pdf>.

## Installation

You can install the the development version of smimodel from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("nuwani-palihawadana/smimodel")
```

***Important:***

The estimation of SMI model is based on an iterative procedure that is
developed using mixed integer programming to solve an
$\ell_{0}$-regularised nonlinear least squares optimisation problem with
linear constraints. We use the commercial MIP solver
[***Gurobi***](https://www.gurobi.com/) to solve the mixed integer
quadratic programs related to the proposed SMI Modelling algorithm.

Hence, users of *smimodel* are required to have an active installation
of Gurobi in their local machines. Since Gurobi is a commercial
software, the users need to obtain a license for Gurobi before they can
install/use it. Please refer to the vignette *“Gurobi installation
help”* for installation help.

## Example

This is a basic example, which shows you how to estimate a SMI model,
and obtain predictions on a test set:

``` r
library(smimodel)
library(dplyr)
library(ROI)
library(tibble)
library(tidyr)
library(tsibble)

# Simulate data
n = 1205
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(n)) %>%
  mutate(
    # Add x_lags
    x_lag = lag_matrix(x_lag_000, 5)) %>%
  unpack(x_lag, names_sep = "_") %>%
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
# Test set
sim_test <- sim_data[1001:1200, ]

# Index variables
index.vars <- colnames(sim_data)[3:8]
```

``` r
# Estimate SMI model (with "PPR" initialisation)
smimodel_ppr <- model_smimodel(data = sim_train,
                               yvar = "y1",
                               index.vars = index.vars,
                               initialise = "ppr")

# Fitted SMI model
smimodel_ppr$fit[[1]]$best
```

``` r
# Plot estimated smooths
autoplot(object = smimodel_ppr)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r

# Obtain residuals and fitted values
augment(x = smimodel_ppr)
#> # A tibble: 1,000 × 5
#>    Index New     Old  .resid .fitted
#>    <int> <chr> <dbl>   <dbl>   <dbl>
#>  1     6 1         1  0.0859   0.495
#>  2     7 1         1  0.0766   0.737
#>  3     8 1         1  0.137    3.70 
#>  4     9 1         1  0.0173   1.15 
#>  5    10 1         1  0.118    0.942
#>  6    11 1         1 -0.201    3.65 
#>  7    12 1         1 -0.0471   1.83 
#>  8    13 1         1  0.0111   1.28 
#>  9    14 1         1 -0.0464   2.46 
#> 10    15 1         1  0.0243   0.261
#> # ℹ 990 more rows
```

``` r
# Obtain predictions on the test set
predict(object = smimodel_ppr, newdata = sim_test)$.predict
```

## License

This package is free and open source software, licensed under GPL-3.
