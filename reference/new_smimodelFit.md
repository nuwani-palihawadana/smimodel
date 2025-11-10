# Constructor function for the class `smimodelFit`

Constructs an object of class `smimodelFit` using the information passed
to arguments.

## Usage

``` r
new_smimodelFit(
  data,
  yvar,
  neighbour = 0,
  family = gaussian(),
  index.vars,
  initialise = c("additive", "linear", "userInput"),
  index.ind = NULL,
  index.coefs = NULL,
  s.vars = NULL,
  linear.vars = NULL
)
```

## Arguments

- data:

  Training data set on which models will be trained. Must be a data set
  of class `tsibble`.(Make sure there are no additional date or time
  related variables except for the `index` of the `tsibble`).

- yvar:

  Name of the response variable as a character string.

- neighbour:

  `neighbour` argument passed from the outer function.

- family:

  A description of the error distribution and link function to be used
  in the model (see [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`family`](https://rdrr.io/r/stats/family.html)).

- index.vars:

  A `character` vector of names of the predictor variables for which
  indices should be estimated.

- initialise:

  The model structure with which the estimation process should be
  initialised. The default is "additive", where the initial model will
  be a nonparametric additive model. The other options are "linear" -
  linear regression model (i.e. a special case single-index model, where
  the initial values of the index coefficients are obtained through a
  linear regression), and "userInput" - user specifies the initial model
  structure (i.e. the number of indices and the placement of index
  variables among indices) and the initial index coefficients through
  `index.ind` and `index.coefs` arguments respectively.

- index.ind:

  If `initialise = "userInput"`: an `integer` vector that assigns group
  index for each predictor in `index.vars`.

- index.coefs:

  If `initialise = "userInput"`: a `numeric` vector of index
  coefficients.

- s.vars:

  A `character` vector of names of the predictor variables for which
  splines should be fitted individually (rather than considering as part
  of an index).

- linear.vars:

  A `character` vector of names of the predictor variables that should
  be included linearly into the model.

## Value

A list of initial model information. For descriptions of the list
elements see
[`make_smimodelFit`](https://nuwani-palihawadana.github.io/smimodel/reference/make_smimodelFit.md)).
