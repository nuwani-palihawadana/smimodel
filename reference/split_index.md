# Splitting predictors into multiple indices

Splits a given number of predictors into a given number of indices.

## Usage

``` r
split_index(num_pred, num_ind)
```

## Arguments

- num_pred:

  Number of predictors.

- num_ind:

  Number of indices.

## Value

A list containing the following components:

- index:

  An `integer` vector that assigns group indices for each predictor.

- index_positions:

  A list of length = `num_ind` that indicates which predictors belong to
  which index.
