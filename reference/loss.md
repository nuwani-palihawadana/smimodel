# Calculating the loss of the MIP used to estimate a SMI model

Calculates the value of the objective function (loss function) of the
mixed integer program used to estimate a SMI model.

## Usage

``` r
loss(Y, Yhat, alpha, lambda0, lambda2)
```

## Arguments

- Y:

  Column matrix of response.

- Yhat:

  Predicted value of the response.

- alpha:

  Vector of index coefficients.

- lambda0:

  Penalty parameter for L0 penalty.

- lambda2:

  Penalty parameter for L2 penalty.

## Value

A `numeric`.
