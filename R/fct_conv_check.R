#' fct_conv_check
#'
#' @export
#'
fct_conv_check <- function(x, y, k, N, clust_assign, selection, alpha, beta, y_sparse, rank, max_rank, ll){
  new_assign_i <- clust_assign
  new_assign <- clust_assign
  new_ll <- ll
  change <- FALSE
  lik_increase <- rep(-Inf, length(clust_assign))
  for(i in 1:length(clust_assign)){
    if(i %% 10 == 0){print(round(i/length(clust_assign),3))}
    # if(!change){
      assign_i <- clust_assign[i]
      max_ll <- rep(-Inf, k)
      
      for(j in 1:k){
        if(assign_i!=j){
          temp_assign <- clust_assign
          temp_assign[i] <- j
          gamma_model <- fct_gamma(x, y, k, N, temp_assign, selection, alpha, beta, y_sparse, rank, max_rank)
          gamma <- gamma_model$gamma
          weighted_ll <- fct_weighted_ll(gamma)
          if(weighted_ll > ll){
            change <- TRUE
            max_ll[j] <- weighted_ll
            print("CHANGE")
          }
        }
        
      }
      if(change){
        
        new_assign_i[i] <- which.max(max_ll)
        lik_increase[i] <- max(max_ll)
      }
    # }
  }
  new_ll <- max(lik_increase)
  new_assign[which.max(lik_increase)] <- new_assign_i[which.max(lik_increase)]
  return(list(clust_assign = new_assign, new_ll = new_ll))
}