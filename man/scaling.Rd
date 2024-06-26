% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling.R
\name{scaling}
\alias{scaling}
\title{Scale data}
\usage{
scaling(data, index.vars)
}
\arguments{
\item{data}{Training data set on which models will be trained. Should be a
\code{tibble}.}

\item{index.vars}{A character vector of names of the predictor variables for
which indices should be estimated.}
}
\description{
Standardises the columns of the \code{data} corresponding to \code{index.vars}.
}
