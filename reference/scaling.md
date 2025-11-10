# Scale data

Scales the columns of the `data` corresponding to `index.vars`.

## Usage

``` r
scaling(data, index.vars)
```

## Arguments

- data:

  Training data set on which models will be trained. Should be a
  `tibble`.

- index.vars:

  A character vector of names of the predictor variables for which
  indices should be estimated.

## Value

A list containing the following components:

- scaled_data:

  The scaled data set of class `tibble`.

- scaled_info:

  A named `numeric` vector of standard deviations of `index.vars` that
  were used to scale the corresponding columns of `data`.
