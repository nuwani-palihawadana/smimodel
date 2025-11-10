# Create lags or leads of a matrix

This is a wrapper for the function
[`conformalForecast::lagmatrix`](https://xqnwang.github.io/conformalForecast/reference/lagmatrix.html).Find
a shifted version of a matrix, adjusting the time base backward (lagged)
or forward (leading) by a specified number of observations for each
column.

## Usage

``` r
leadlagMat(x, lag)
```

## Arguments

- x:

  A matrix or multivariate time series.

- lag:

  A vector of lags (positive values) or leads (negative values) with a
  length equal to the number of columns of `x`.

## Value

A matrix with the same class and size as `x`.

## Examples

``` r
x <- matrix(rnorm(20), nrow = 5, ncol = 4)

# Create lags of a matrix
leadlagMat(x, c(0, 1, 2, 3))
#>             [,1]       [,2]       [,3]       [,4]
#> [1,] -0.02854676         NA         NA         NA
#> [2,] -0.04287046 -1.5487528         NA         NA
#> [3,]  1.36860228  0.5846137 -0.5023235         NA
#> [4,] -0.22577099  0.1238542 -0.3332074 0.44820978
#> [5,]  1.51647060  0.2159416 -1.0185754 0.05300423
#> attr(,"class")
#> [1] "matrix" "array" 

# Create leads of a matrix
leadlagMat(x, c(0, -1, -2, -3))
#>             [,1]      [,2]       [,3]       [,4]
#> [1,] -0.02854676 0.5846137 -1.0185754  2.0500847
#> [2,] -0.04287046 0.1238542 -1.0717912 -0.4910312
#> [3,]  1.36860228 0.2159416  0.3035286         NA
#> [4,] -0.22577099 0.3796395         NA         NA
#> [5,]  1.51647060        NA         NA         NA
#> attr(,"class")
#> [1] "matrix" "array" 
```
