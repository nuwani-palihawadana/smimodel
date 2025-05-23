% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_smimodel.R
\name{inner_update}
\alias{inner_update}
\title{Updating index coefficients and non-linear functions iteratively}
\usage{
inner_update(
  x,
  data,
  yvar,
  family = gaussian(),
  index.vars,
  s.vars,
  linear.vars,
  num_ind,
  dgz,
  alpha_old,
  lambda0 = 1,
  lambda2 = 1,
  M = 10,
  max.iter = 50,
  tol = 0.001,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE
)
}
\arguments{
\item{x}{Fitted \code{gam}.}

\item{data}{Training data set on which models will be trained. Should be a
\code{tsibble}.}

\item{yvar}{Name of the response variable as a character string.}

\item{family}{A description of the error distribution and link function to be
used in the model (see \code{\link{glm}} and \code{\link{family}}).}

\item{index.vars}{A \code{character} vector of names of the predictor
variables for which indices should be estimated.}

\item{s.vars}{A \code{character} vector of names of the predictor variables
for which splines should be fitted individually (rather than considering as
part of an index).}

\item{linear.vars}{A \code{character} vector of names of the predictor
variables that should be included linearly into the model.}

\item{num_ind}{Number of indices.}

\item{dgz}{The \code{tibble} of derivatives of the estimated smooths.}

\item{alpha_old}{Current vector of index coefficients.}

\item{lambda0}{Penalty parameter for L0 penalty.}

\item{lambda2}{Penalty parameter for L2 penalty.}

\item{M}{Big-M value to be used in MIP.}

\item{max.iter}{Maximum number of MIP iterations performed to update index
coefficients for a given model.}

\item{tol}{Tolerance for loss.}

\item{TimeLimit}{A limit for the total time (in seconds) expended in a single
MIP iteration.}

\item{MIPGap}{Relative MIP optimality gap.}

\item{NonConvex}{The strategy for handling non-convex quadratic objectives or
non-convex quadratic constraints in Gurobi solver.}

\item{verbose}{The option to print detailed solver output.}
}
\value{
A list containing following elements: \item{best_alpha}{The vector of
best index coefficient estimates.} \item{min_loss}{Minimum value of the
objective function(loss).}
\item{index.ind}{An \code{integer} vector that assigns group index for each
predictor, corresponding to \code{best_alpha}.}
\item{ind_pos}{A list that indicates which predictors belong to which index,
corresponding to \code{best_alpha}.}
\item{X_new}{A matrix of selected predictor variables, corresponding to
\code{best_alpha}.}
}
\description{
Iteratively updates index coefficients and non-linear functions using mixed
integer programming. (A helper function used within
\code{\link{update_smimodelFit}}; users are not expected to directly call
this function.)
}
