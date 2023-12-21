#' Splitting predictors into multiple indices
#'
#' Splits a given number of predictors into a given number of indices.
#'
#' @param num_pred Number of predictors.
#' @param num_ind Number of indices.
#'
#' @export
split_index <- function(num_pred, num_ind){
  split_list <- numeric()
  rem_pred <- num_pred
  rem_ind <- num_ind
  count = 1
  while(count <= num_ind){
    split_num <- ceiling(rem_pred/rem_ind)
    rem_pred <- rem_pred - split_num
    rem_ind <- rem_ind - 1
    split1 <- rep(count, split_num)
    split_list <- c(split_list, split1)
    count <- count + 1
  }
  split_pos <- split(1:length(split_list), split_list) 
  output <- list("index" = split_list, "index_positions" = split_pos)
  return(output)
}