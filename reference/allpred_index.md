# Constructing index coefficient vectors with all predictors in each index

Constructs vectors of coefficients for each index including a
coefficient for all the predictors that are entering indices. i.e. if a
coefficient is not provided for a particular predictor in a particular
index, the function will replace the missing coefficient with a zero.

## Usage

``` r
allpred_index(num_pred, num_ind, ind_pos, alpha)
```

## Arguments

- num_pred:

  Number of predictors.

- num_ind:

  Number of indices.

- ind_pos:

  A list of length = `num_ind` that indicates which predictors belong to
  which index.

- alpha:

  A vector of index coefficients.

## Value

A list containing the following components:

- alpha_init_new:

  A `numeric` vector of index coefficients.

- index:

  An `integer` vector that assigns group indices for each predictor.

- index_positions:

  A list of length = `num_ind` that indicates which predictors belong to
  which index.
