# Truncating predictors to be in the in-sample range

Truncates predictors to be in the in-sample range to avoid spline
extrapolation.

## Usage

``` r
truncate_vars(range.object, data, cols.trunc)
```

## Arguments

- range.object:

  A matrix containing range of each predictor variable. Should be a
  matrix with two rows for min and max, and the columns should
  correspond to variables.

- data:

  Out-of-sample data set of which variables should be truncated.

- cols.trunc:

  Column names of the variables to be truncated.
