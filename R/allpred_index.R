#' Constructing index coefficient vectors with all predictors in each index
#'
#' Constructs vectors of coefficients for each index including a
#' coefficient for all the predictors that are entering indices. i.e. if a
#' coefficient is not provided for a particular predictor in a particular index,
#' the function will replace the missing coefficient with a zero.
#'
#' @param num_pred Number of predictors.
#' @param num_ind Number of indices.
#' @param ind_pos A list of which the length = `num_ind` that indicates which
#'   predictors belong to which index.
#' @param alpha A vector of index coefficients.

allpred_index <- function(num_pred, num_ind, ind_pos, alpha){
  init_list <- vector(mode = "list", length = num_ind)
  index <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    if(length(ind_pos[[i]]) == 1){
      temp <- i
    }else{
      temp <- unlist(lapply(1:length(ind_pos[[i]]), function(x) paste0(i, x))) 
    }
    if(length(alpha) == (num_ind*num_pred)){
      init_list[[i]] <- alpha[match(temp, names(alpha))]
    }else{
      init_list[[i]] <- numeric(length = num_pred)
      init_list[[i]][ind_pos[[i]]] <- alpha[match(temp, names(alpha))]
    }
    index[[i]] <- rep(i, num_pred)
  }
  alpha_init <- unlist(init_list)
  index <- unlist(index)
  ind_pos <- split(1:length(index), index)
  output <- list("alpha_init_new" = alpha_init,
                 "index" = index,
                 "index_positions" = ind_pos)
  return(output)
}
