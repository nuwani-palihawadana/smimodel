#' Scale data
#'
#' Standardises the columns of the `data` corresponding to `index.vars`.
#'
#' @param data Training data set on which models will be trained. Should be a
#'   `tibble`.
#' @param index.vars A character vector of names of the predictor variables for
#'   which indices should be estimated.

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