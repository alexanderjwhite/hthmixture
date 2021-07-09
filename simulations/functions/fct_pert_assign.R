#' Perturb Assignments
#'
#' @param clust_assign current assignment
#' @param lik_store likelihood vector
#' @param N observations
#' @param k clusters
#'
#' @return updated assignment
#' @export
#'
#' @import dplyr
#'
#' @examples
fct_pert_assign <- function(clust_assign, lik_store, N, k){
  
  perc_pert <- fct_boltzmann_pert(lik_store)
  num_pert <- (perc_pert*N) %>% floor()
  pert_samples <- sample(1:N, num_pert)
  
  if(num_pert!=0){
    unique_vals <- 1:k
    for(i in 1:length(pert_samples)){
      
      obs <- pert_samples[i]
      orig <- clust_assign[obs]
      vals <- unique_vals[unique_vals != orig]
      
      if(length(vals) == 1){
        
        clust_assign[obs] <- vals
        
      } else if (length(vals) == 0) {
        
        clust_assign[obs] <- orig
        
      } else {
        
        clust_assign[obs] <- sample(vals)
        
      }
    }
  }
  return(clust_assign)
}
