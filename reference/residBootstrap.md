# Generate multiple single season block bootstrap series

Generates multiple replications of single season block bootstrap series.

## Usage

``` r
residBootstrap(x, season.period = 1, m = 1, num.bootstrap = 1000)
```

## Arguments

- x:

  A series of residuals from which bootstrap series to be generated.

- season.period:

  Length of the seasonal period.

- m:

  Multiplier. (Block size = `season.period * m`)

- num.bootstrap:

  Number of bootstrap series to be generated.

## Value

A matrix of bootstrapped series.
