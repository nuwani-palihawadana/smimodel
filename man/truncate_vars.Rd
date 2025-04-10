% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{truncate_vars}
\alias{truncate_vars}
\title{Truncating predictors to be in the in-sample range}
\usage{
truncate_vars(range.object, data, cols.trunc)
}
\arguments{
\item{range.object}{A matrix containing range of each predictor variable.
Should be a matrix with two rows for min and max, and the columns should
correspond to variables.}

\item{data}{Out-of-sample data set of which variables should be truncated.}

\item{cols.trunc}{Column names of the variables to be truncated.}
}
\description{
Truncates predictors to be in the in-sample range to avoid spline
extrapolation.
}
