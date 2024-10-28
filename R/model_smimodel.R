#' Sparse Multiple Index (SMI) Models
#'
#' Fits nonparametric multiple index model(s), with simultaneous predictor
#' selection (hence "sparse") and predictor grouping. Possible to fit multiple
#' SMI models based on a grouping variable.
#'
#' @param data Training data set on which models will be trained. Must be a data
#'   set of class \code{tsibble}.(Make sure there are no additional date or time
#'   related variables except for the \code{index} of the \code{tsibble}). If
#'   multiple models are fitted, the grouping variable should be the \code{key}
#'   of the \code{tsibble}. If a \code{key} is not specified, a dummy key with
#'   only one level will be created.
#' @param yvar Name of the response variable as a character string.
#' @param neighbour If multiple models are fitted: Number of neighbours of each
#'   key (i.e. grouping variable) to be considered in model fitting to handle
#'   smoothing over the key. Should be an \code{integer}. If \code{neighbour =
#'   x}, \code{x} number of keys before the key of interest and \code{x} number
#'   of keys after the key of interest are grouped together for model fitting.
#'   The default is \code{neighbour = 0} (i.e. no neighbours are considered for
#'   model fitting).
#' @param family A description of the error distribution and link function to be
#'   used in the model (see \code{\link{glm}} and \code{\link{family}}).
#' @param index.vars A \code{character} vector of names of the predictor
#'   variables for which indices should be estimated.
#' @param initialise The model structure with which the estimation process
#'   should be initialised. The default is \code{"ppr"}, where the initial model
#'   is derived from projection pursuit regression. The other options are
#'   \code{"additive"} - nonparametric additive model, \code{"linear"} - linear
#'   regression model (i.e. a special case single-index model, where the initial
#'   values of the index coefficients are obtained through a linear regression),
#'   \code{"multiple"} - multiple models are fitted starting with different
#'   initial models (number of indices = \code{num_ind}; \code{num_models}
#'   random instances of the model (i.e. the predictor assignment to indices and
#'   initial index coefficients are generated randomly) are considered), and the
#'   final optimal model with the lowest loss is returned, and
#'   \code{"userInput"} - user specifies the initial model structure (i.e. the
#'   number of indices and the placement of index variables among indices) and
#'   the initial index coefficients through \code{index.ind} and
#'   \code{index.coefs} arguments respectively.
#' @param num_ind If \code{initialise = "ppr"} or \code{"multiple"}: an
#'   \code{integer} that specifies the number of indices to be used in the
#'   model(s). The default is \code{num_ind = 5}.
#' @param num_models If \code{initialise = "multiple"}: an \code{integer} that
#'   specifies the number of starting models to be checked. The default is
#'   \code{num_models = 5}.
#' @param seed If \code{initialise = "multiple"}: the seed to be set when
#'   generating random starting points.
#' @param index.ind If \code{initialise = "userInput"}: an \code{integer} vector
#'   that assigns group index for each predictor in \code{index.vars}.
#' @param index.coefs If \code{initialise = "userInput"}: a \code{numeric}
#'   vector of index coefficients.
#' @param s.vars A \code{character} vector of names of the predictor variables
#'   for which splines should be fitted individually (rather than considering as
#'   part of an index).
#' @param linear.vars A \code{character} vector of names of the predictor
#'   variables that should be included linearly into the model.
#' @param lambda0 Penalty parameter for L0 penalty.
#' @param lambda2 Penalty parameter for L2 penalty.
#' @param M Big-M value used in MIP.
#' @param max.iter Maximum number of MIP iterations performed to update index
#'   coefficients for a given model.
#' @param tol Tolerance for the objective function value (loss) of MIP.
#' @param tolCoefs Tolerance for coefficients.
#' @param TimeLimit A limit for the total time (in seconds) expended in a single
#'   MIP iteration.
#' @param MIPGap Relative MIP optimality gap.
#' @param NonConvex The strategy for handling non-convex quadratic objectives or
#'   non-convex quadratic constraints in Gurobi solver.
#' @param verbose The option to print detailed solver output.
#' @return  An object of class \code{smimodel}. This is a \code{tibble} with two
#'   columns: \item{key}{The level of the grouping variable (i.e. key of the
#'   training data set).} \item{fit}{Information of the fitted model
#'   corresponding to the \code{key}.}
#'   Each row of the column \code{fit} contains a list with two elements:
#'   \item{initial}{A list of information of the model initialisation. (For
#'   descriptions of the list elements see \code{\link{make_smimodelFit}}).}
#'   \item{best}{A list of information of the final optimised model. (For
#'   descriptions of the list elements see \code{\link{make_smimodelFit}}).}
#'
#' @details Sparse Multiple Index (SMI) model is a semi-parametric model that
#'   can be
#' written as \deqn{y_{i} = \beta_{0} +
#' \sum_{j = 1}^{p}g_{j}(\bm{\alpha}_{j}^{T}\bm{x}_{ij}) +
#' \sum_{k = 1}^{d}f_{k}(w_{ik}) + \bm{\theta}^{T}\bm{u}_{i} + \varepsilon_{i},
#' \quad i = 1, \dots, n,} where \eqn{y_{i}} is the univariate response,
#' \eqn{\beta_{0}} is the model intercept, \eqn{\bm{x}_{ij} \in
#' \mathbb{R}^{l_{j}}}, \eqn{j = 1, \dots, p} are \eqn{p} subsets of predictors
#'   entering indices, \eqn{\bm{\alpha}_{j}} is a vector of index coefficients
#'   corresponding to the index \eqn{h_{ij} = \bm{\alpha}_{j}^{T}\bm{x}_{ij}},
#'   and \eqn{g_{j}} is a smooth nonlinear function (estimated by a penalised
#'   cubic regression spline). The model also allows for predictors that do not
#'   enter any indices, including covariates \eqn{w_{ik}} that relate to the
#'   response through nonlinear functions \eqn{f_{k}}, \eqn{k = 1, \dots, d},
#'   and linear covariates \eqn{\bm{u}_{i}}.
#'
#'   In the model formulation related to this implementation, both the number of
#'   indices \eqn{p} and the predictor grouping among indices are assumed to be
#'   unknown prior to model estimation. Suppose we observe \eqn{y_1,\dots,y_n},
#'   along with a set of potential predictors, \eqn{\bm{x}_1,\dots,\bm{x}_n},
#'   with each vector \eqn{\bm{x}_i} containing \eqn{q} predictors. This
#'   function implements algorithmic variable selection for index variables
#'   (i.e. predictors entering indices) of the SMI model by allowing for zero
#'   index coefficients for predictors. Non-overlapping predictors among indices
#'   are assumed (i.e. no predictor enters more than one index). For algorithmic
#'   details see reference.
#'
#' @references Palihawadana, N.K., Hyndman, R.J. & Wang, X. (2024). Sparse
#'   Multiple Index Models for High-Dimensional Nonparametric Forecasting.
#'   \url{https://www.monash.edu/business/ebs/research/publications/ebs/2024/wp16-2024.pdf}.
#'
#' @seealso \code{\link{greedy_smimodel}}
#'
#' @examples
#' library(dplyr)
#' library(ROI)
#' library(tibble)
#' library(tidyr)
#' library(tsibble)
#'
#' # Simulate data
#' n = 1005
#' set.seed(123)
#' sim_data <- tibble(x_lag_000 = runif(n)) |>
#'   mutate(
#'     # Add x_lags
#'     x_lag = lag_matrix(x_lag_000, 5)) |>
#'   unpack(x_lag, names_sep = "_") |>
#'   mutate(
#'     # Response variable
#'     y1 = (0.9*x_lag_000 + 0.6*x_lag_001 + 0.45*x_lag_003)^3 + rnorm(n, sd = 0.1),
#'     # Add an index to the data set
#'     inddd = seq(1, n)) |>
#'   drop_na() |>
#'   select(inddd, y1, starts_with("x_lag")) |>
#'   # Make the data set a `tsibble`
#'   as_tsibble(index = inddd)
#'
#' # Index variables
#' index.vars <- colnames(sim_data)[3:8]
#'
#' # Model fitting
#' smimodel_ppr <- model_smimodel(data = sim_data,
#'                                yvar = "y1",
#'                                index.vars = index.vars,
#'                                initialise = "ppr")
#'
#' # Best (optimised) fitted model
#' smimodel_ppr$fit[[1]]$best
#'
#' @export
model_smimodel <- function(data, yvar, neighbour = 0, family = gaussian(), 
                           index.vars, 
                           initialise = c("ppr", "additive", "linear", 
                                          "multiple", "userInput"),
                           num_ind = 5, num_models = 5, seed = 123, 
                           index.ind = NULL, index.coefs = NULL, 
                           s.vars = NULL, linear.vars = NULL, 
                           lambda0 = 1, lambda2 = 1, 
                           M = 10, max.iter = 50, tol = 0.001, tolCoefs = 0.001,
                           TimeLimit = Inf, MIPGap = 1e-4, 
                           NonConvex = -1, verbose = FALSE){
  stopifnot(tsibble::is_tsibble(data))
  initialise <- match.arg(initialise)
  data1 <- data
  data_index <- index(data1)
  data_key <- key(data1)
  if (length(key(data1)) == 0) {
    data1 <- data1 |>
      dplyr::mutate(dummy_key = rep(1, NROW(data1))) |>
      tsibble::as_tsibble(index = data_index, key = dummy_key)
    data_key <- key(data1)
  }
  key11 <- key(data1)[[1]]
  key_unique <- unique(as.character(sort(dplyr::pull((data1[, {{ key11 }}])[, 1]))))
  key_num <- seq_along(key_unique)
  ref <- data.frame(key_unique, key_num)
  data1 <- data1 |>
    dplyr::mutate(
      num_key = as.numeric(factor(as.character({{ key11 }}), levels = key_unique))
    )
  smimodels_list <- vector(mode = "list", length = NROW(ref))
  for (i in seq_along(ref$key_num)){
    print(paste0('model ', paste0(i)))
    df_cat <- data1 |>
      dplyr::filter((abs(num_key - ref$key_num[i]) <= neighbour) |
                      (abs(num_key - ref$key_num[i] + NROW(ref)) <= neighbour) |
                      (abs(num_key - ref$key_num[i] - NROW(ref)) <= neighbour)) 
    smimodels_list[[i]] <- smimodel.fit(data = df_cat, yvar = yvar, 
                                        neighbour = neighbour,
                                        family = family,
                                        index.vars = index.vars, 
                                        initialise = initialise, 
                                        num_ind = num_ind, num_models = num_models, 
                                        seed = seed,
                                        index.ind = index.ind, 
                                        index.coefs = index.coefs,
                                        s.vars = s.vars,
                                        linear.vars = linear.vars,
                                        lambda0 = lambda0, lambda2 = lambda2, 
                                        M = M, max.iter = max.iter, 
                                        tol = tol, tolCoefs = tolCoefs,
                                        TimeLimit = TimeLimit, MIPGap = MIPGap,
                                        NonConvex = NonConvex, verbose = verbose)
  }
  data_list <- list(key_unique, smimodels_list)
  models <- tibble::as_tibble(
    x = data_list, .rows = length(data_list[[1]]),
    .name_repair = ~ make.names(names = c("key", "fit"))
  )
  class(models) <- c("smimodel", "tbl_df", "tbl", "data.frame")
  return(models)
}