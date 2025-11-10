# Point estimate accuracy measures

Point estimate accuracy measures

## Usage

``` r
MAE(residuals, na.rm = TRUE, ...)

MSE(residuals, na.rm = TRUE, ...)

point_measures
```

## Format

An object of class `list` of length 2.

## Arguments

- residuals:

  A vector of residuals from either the validation or test data.

- na.rm:

  If `TRUE`, remove the missing values before calculating the accuracy
  measure.

- ...:

  Additional arguments for each measure.

## Value

For the individual functions (`MAE`, `MSE`), returns a single numeric
scalar giving the requested accuracy measure.

For the exported object `point_measures`, returns a named list of
functions that can be supplied to higher-level accuracy routines.

## Examples

``` r
set.seed(123)
ytrain <- rnorm(100)
ytest  <- rnorm(30)
yhat   <- ytest + rnorm(30, sd = 0.3)
resid   <- ytest - yhat

MAE(resid)
#> [1] 0.2840219
MSE(resid)
#> [1] 0.1116308
```
