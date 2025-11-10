# Single season block bootstrap

Generates a single replication of single season block bootstrap series.

## Usage

``` r
seasonBootstrap(x, season.period = 1, m = 1)
```

## Arguments

- x:

  A series of residuals from which bootstrap series to be generated.

- season.period:

  Length of the seasonal period.

- m:

  Multiplier. (Block size = `season.period * m`)

## Value

A `numeric` vector.
