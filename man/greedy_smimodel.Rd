% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/greedy_smimodel.R
\name{greedy_smimodel}
\alias{greedy_smimodel}
\title{SMI model estimation through a greedy search for penalty parameters}
\usage{
greedy_smimodel(
  data,
  val.data,
  yvar,
  neighbour = 0,
  family = gaussian(),
  index.vars,
  initialise = c("ppr", "additive", "linear", "multiple", "userInput"),
  num_ind = 5,
  num_models = 5,
  seed = 123,
  index.ind = NULL,
  index.coefs = NULL,
  s.vars = NULL,
  linear.vars = NULL,
  lambda0_seq,
  lambda2_seq,
  lambda_step,
  lambda0_start_seq,
  lambda2_start_seq,
  refit = TRUE,
  M = 10,
  max.iter = 50,
  tol = 0.001,
  tolCoefs = 0.001,
  TimeLimit = Inf,
  MIPGap = 1e-04,
  NonConvex = -1,
  verbose = FALSE,
  parallel = FALSE,
  workers = NULL,
  recursive = FALSE,
  recursive_colRange = NULL
)
}
\arguments{
\item{data}{Training data set on which models will be trained. Should be a
\code{tsibble}.}

\item{val.data}{Validation data set. (The data set on which the penalty
parameter selection will be performed.) Must be a data set of class
\code{tsibble}. (Once the penalty parameter selection is completed, the best
model will be re-fitted for the combined data set \code{data} + \code{val.data}.)}

\item{yvar}{Name of the response variable as a character string.}

\item{neighbour}{If multiple models are fitted: Number of neighbours of each
key (i.e. grouping variable) to be considered in model fitting to handle
smoothing over the key. Should be an integer. If \code{neighbour = x}, \code{x}
number of keys before the key of interest and \code{x} number of keys after the
key of interest are grouped together for model fitting. The default is \code{0}
(i.e. no neighbours are considered for model fitting).}

\item{family}{A description of the error distribution and link function to be
used in the model (see \code{\link{glm}} and \code{\link{family}}).}

\item{index.vars}{A character vector of names of the predictor variables for
which indices should be estimated.}

\item{initialise}{The model structure with which the estimation process
should be initialised. The default is "ppr", where the initial model is
derived from projection pursuit regression. The other options are
"additive" - nonparametric additive model, "linear" - linear regression
model (i.e. a special case single-index model, where the initial values of
the index coefficients are obtained through a linear regression),
"multiple" - multiple models are fitted starting with different initial
models (number of indices = \code{num_ind}; \code{num_models} random instances of the
model (i.e. the predictor assignment to indices and initial index
coefficients are generated randomly) are considered), and the final optimal
model with the lowest loss is returned, and "userInput" - user specifies
the initial model structure (i.e. the number of indices and the placement
of index variables among indices) and the initial index coefficients
through \code{index.ind} and \code{index.coefs} arguments respectively.}

\item{num_ind}{If \verb{initialise = "ppr" or "multiple"}: an integer that
specifies the number of indices to be used in the model(s). The default is
5.}

\item{num_models}{If \code{initialise = "multiple"}: an integer that specifies the
number of starting models to be checked. The default is 5.}

\item{seed}{If \code{initialise = "multiple"}: the seed to be set when generating
random starting points.}

\item{index.ind}{If \code{initialise = "userInput"}: an integer vector that
assigns group index for each predictor in \code{index.vars}.}

\item{index.coefs}{If \code{initialise = "userInput"}: a numeric vector of index
coefficients.}

\item{s.vars}{A character vector of names of the predictor variables for
which splines should be fitted individually (rather than considering as a
part of an index considered in \code{index.vars}).}

\item{linear.vars}{A character vector of names of the predictor variables
that should be included linearly into the model.}

\item{lambda0_seq}{A numeric vector of candidate values for lambda0 (penalty
parameter for L0 penalty).}

\item{lambda2_seq}{A numeric vector of candidate values for lambda2 (penalty
parameter for L2 penalty).}

\item{lambda_step}{Step size between two adjacent values in \code{lambda0_seq} and
\code{lambda2_seq}.}

\item{lambda0_start_seq}{A subset from \code{lambda0_seq} as candidate starting
points for the greedy search.}

\item{lambda2_start_seq}{A subset from \code{lambda2_seq} as candidate starting
points for the greedy search.}

\item{refit}{Whether to refit the model combining training and validation
sets after parameter tuning. If \code{FALSE}, the final model will be estimated
only on the training set.}

\item{M}{Big-M value used in MIP.}

\item{max.iter}{Maximum number of MIP iterations performed to update index
coefficients for a given model.}

\item{tol}{Tolerance for loss.}

\item{tolCoefs}{Tolerance for coefficients.}

\item{TimeLimit}{A limit for the total time (in seconds) expended in a single
MIP iteration.}

\item{MIPGap}{Relative MIP optimality gap.}

\item{NonConvex}{The strategy for handling non-convex quadratic objectives or
non-convex quadratic constraints in Gurobi solver.}

\item{verbose}{The option to print detailed solver output.}

\item{parallel}{The option to use parallel processing in fitting \code{smimodel}s
for different penalty parameter combinations.}

\item{workers}{If \code{parallel = TRUE}: Number of cores to use.}

\item{recursive}{Whether to obtain recursive forecasts or not (default -
FALSE).}

\item{recursive_colRange}{If \code{recursive = TRUE}, The range of column numbers
in \code{val.data} to be filled with forecasts.}
}
\description{
Performs a greedy search over a given grid of penalty parameter combinations
(lambda0, lambda2), and fits SMI model(s) with best (lowest validation set
MSE) penalty parameter combinations. If a grouping variable is used, penalty
parameters are tuned separately for each individual model.
}
\examples{
library(dplyr)
library(ROI)
library(tibble)
library(tidyr)
library(tsibble)
n = 1205
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
# Training set
sim_train <- sim_data[1:1000, ]
# Validation set
sim_val <- sim_data[1001:1200, ]
# Index variables
index.vars <- colnames(sim_data)[3:8]
# Penalty parameter values to search
# L0 penalty
lambda0 = seq(1, 12, by = 1)
# L2 penalty
lambda2 = seq(0, 12, by = 1)
# Full grid
grid1 <- expand.grid(lambda0, lambda2)
# Starting point options
starting <- grid1[c(1, 6, 12, 73, 78, 84, 145, 150, 156), ]
# L0 penalty
lambda0_start = as.numeric(unique(unlist(starting[1])))
# L2 penalty
lambda2_start = as.numeric(unique(unlist(starting[2])))
# Model fitting
smi_greedy <- greedy_smimodel(data = sim_train,
                              val.data = sim_val,
                              yvar = "y1",
                              index.vars = index.vars,
                              initialise = "additive",
                              lambda0_seq = lambda0,
                              lambda2_seq = lambda2,
                              lambda_step = 1,
                              lambda0_start_seq = lambda0_start,
                              lambda2_start_seq = lambda2_start)
smi_greedy$fit[[1]]
}
