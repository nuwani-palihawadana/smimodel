% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_ppr.R
\name{model_ppr}
\alias{model_ppr}
\title{Projection Pursuit Regression (PPR) models}
\usage{
model_ppr(data, yvar, neighbour = 0, index.vars, num_ind = 5, ...)
}
\arguments{
\item{data}{Training data set on which models will be trained. Must be a data
set of class \code{tsibble}.(Make sure there are no additional
date/time/date-time/yearmonth/POSIXct/POSIXt variables except for the
\code{index} of the \code{tsibble}). If multiple models are fitted, the grouping
variable should be the key of the \code{tsibble}.}

\item{yvar}{Name of the response variable as a character string.}

\item{neighbour}{If multiple models are fitted: Number of neighbours of each
key (i.e. grouping variable) to be considered in model fitting to handle
smoothing over the key. Should be an integer. If \code{neighbour = x}, \code{x}
number of keys before the key of interest and \code{x} number of keys after the
key of interest are grouped together for model fitting. The default is \code{0}
(i.e. no neighbours are considered for model fitting).}

\item{index.vars}{A character vector of names of the predictor variables for
which indices should be estimated.}

\item{num_ind}{An integer that specifies the number of indices to be used in
the model(s). (Corresponds to \code{nterms} in \code{stats::ppr()}.)}

\item{...}{Other arguments not currently used. (For more information on other
arguments that can be passed, refer \code{stats::ppr()}.)}
}
\description{
A wrapper for \code{stats::ppr()} enabling multiple PPR models based on a grouping
variable.
}
\examples{
library(dplyr)
library(tibble)
library(tidyr)
library(tsibble)
n = 1005
set.seed(123)
sim_data <- tibble(x_lag_000 = runif(n)) |>
  mutate(
    # Add x_lags
    x_lag = lag_matrix(x_lag_000, 5)) |>
  unpack(x_lag, names_sep = "_") |>
  mutate(
    # Response variable
    y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
    # Add an index to the data set
    inddd = seq(1, n)) |>
  drop_na() |>
  select(inddd, y1, starts_with("x_lag")) |>
  # Make the data set a `tsibble`
  as_tsibble(index = inddd)
# Index variables
index.vars <- colnames(sim_data)[3:8]
# Model fitting
pprModel <- model_ppr(data = sim_data,
                      yvar = "y1",
                      index.vars = index.vars)
pprModel$fit[[1]]
}
