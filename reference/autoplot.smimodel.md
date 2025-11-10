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
