#' Scale data
#'
#' Scales the columns of the \code{data} corresponding to \code{index.vars}.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   \code{tibble}.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.
#' @return A list containing the following components: \item{scaled_data}{The
#'   scaled data set of class \code{tibble}.} \item{scaled_info}{A named numeric
#'   vector of standard deviations of \code{index.vars} that were used to scale
#'   the corresponding columns of \code{data}.}

scaling <- function(data, index.vars){
  scaleData <- scale(data[ , index.vars], center = FALSE, 
                     scale = apply(data[ , index.vars], 2, sd, na.rm = TRUE))
  scaleInfo <- attributes(scaleData)$`scaled:scale`
  data <- data |>
    dplyr::select(-{{index.vars}})
  data <- dplyr::bind_cols(data, scaleData)
  output <- list("scaled_data" = data, "scaled_info" = scaleInfo)
  return(output)
}