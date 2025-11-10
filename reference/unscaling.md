# Unscale a fitted `smimodel`

Transforms back the index coefficients to suit original-scale index
variables if the same were scaled when estimating the `smimodel`
(happens in `initialise = "ppr"` in
[`model_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/model_smimodel.md)
or
[`greedy_smimodel`](https://nuwani-palihawadana.github.io/smimodel/reference/greedy_smimodel.md)).
Users are not expected to directly use this function; usually called
within
[`smimodel.fit`](https://nuwani-palihawadana.github.io/smimodel/reference/smimodel.fit.md).

## Usage

``` r
unscaling(object, scaledInfo)
```

## Arguments

- object:

  A `smimodel` object.

- scaledInfo:

  The list returned from a call of the function
  [`scaling`](https://nuwani-palihawadana.github.io/smimodel/reference/scaling.md).

## Value

A `smimodel` object.
