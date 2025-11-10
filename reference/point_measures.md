# Point estimate accuracy measures

Point estimate accuracy measures

## Usage

``` r
MSE(residuals, na.rm = TRUE, ...)

MAE(residuals, na.rm = TRUE, ...)

point_measures
```

## Format

An object of class `list` of length 2.

## Arguments

- residuals:

  A vector of residuals from either the validation or test data.

- na.rm:

  If `TRUE`, remove the missing values before calculating the accuracy
  measure.

- ...:

  Additional arguments for each measure.
