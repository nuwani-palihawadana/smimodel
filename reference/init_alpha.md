# Initialising index coefficients

Initialises index coefficient vector through linear regression or
penalised linear regression.

## Usage

``` r
init_alpha(
  Y,
  X,
  index.ind,
  init.type = "penalisedReg",
  lambda0 = 1,
  lambda2 = 1,
  M = 10
)
```

## Arguments

- Y:

  Column matrix of response.

- X:

  Matrix of predictors entering indices.

- index.ind:

  An `integer` vector that assigns group index for each predictor.

- init.type:

  Type of initialisation for index coefficients. (`"penalisedReg"` -
  Penalised linear regression; `"reg"` - Linear regression)

- lambda0:

  If `init.type = "penalisedReg"`, penalty parameter for L0 penalty.

- lambda2:

  If `init.type = "penalisedReg"`, penalty parameter for L2 penalty.

- M:

  If `init.type = "penalisedReg"`, the big-M value to be used in the
  MIP.

## Value

A list containing the following components:

- alpha_init:

  Normalised vector of index coefficients.

- alpha_nonNormalised:

  Non-normalised (i.e. prior to normalising) vector of index
  coefficients.
