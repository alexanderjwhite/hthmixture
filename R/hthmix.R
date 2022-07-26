#' Main mixrrr function
#'
#' @param x matrix
#' @param y matrix
#' @param k number of clusters
#' @param maxiter maximum number of iterations
#' @param verbose print progress?
#'
#' @return results
#' @export
#' 
#' @import dplyr purrr
#' @importFrom rlang .data
#'
hthmix <- function(x, y, k, maxiter = 1e3, verbose = TRUE){
  
  N <- x %>% nrow()
  clust_assign <- fct_initialize(k, N)

  model <- fct_hthmix_comp(x, y, k, maxiter, clust_assign)
  ll <- model$ll
  clust_assign <- model$assign
  A <- model$A
  sig_vec <- model$sig_vec
  clust_assign_store <- model$assign_store
 
  result <- tibble(llik = ll, assign = list(final_assign = clust_assign, A = A, sig_vec = sig_vec, assign_store = clust_assign_store))
  return(result)
}
