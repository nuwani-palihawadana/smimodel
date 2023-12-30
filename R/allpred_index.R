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
#'
#' @export
allpred_index <- function(num_pred, num_ind, ind_pos, alpha){
  init_list <- vector(mode = "list", length = num_ind)
  index <- vector(mode = "list", length = num_ind)
  for(i in 1:num_ind){
    # XQ: Error?
    # if(length(alpha) == (num_ind*num_pred)){
    #   init_list[[i]] <- alpha[ind_pos[[i]]]
    # }else{
    #   init_list[[i]] <- numeric(length = num_pred)
    #   init_list[[i]][ind_pos[[i]]] <- alpha[ind_pos[[i]]]
    # }
    # NP: This will give incorrect "alpha_init" (unnecessary zeros will be added) 
    # when there are more than one index (especially when obtaining 
    # "final_smimodel" in "update_smimodel()").
    # init_list[[i]] <- numeric(length = num_pred)
    # init_list[[i]][ind_pos[[i]]] <- alpha[ind_pos[[i]]]
    if(length(alpha) == (num_ind*num_pred)){
      init_list[[i]] <- alpha[startsWith(names(alpha), paste0(i))]
    }else{
      init_list[[i]] <- numeric(length = num_pred)
      init_list[[i]][ind_pos[[i]]] <- alpha[startsWith(names(alpha), paste0(i))]
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
