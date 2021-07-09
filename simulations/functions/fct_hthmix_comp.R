#' Compute mixrrr
#'
#' @param x matrix
#' @param y matrix
#' @param k clusters
#' @param maxiter maximum number of iterations
#' @param lam penalization parameter
#' @param rank rank
#' @param clust_assign initial cluster assignment
#' @param val_frac validation fraction
#' @param penal_search penalization search vector
#'
#' @return results
#' @export
#' 
#' @import dplyr
#'
#' @examples
fct_hthmix_comp <- function(x, y, k, maxiter, lam, rank, clust_assign, val_frac, penal_search){
  
  N <- x %>% 
    nrow()
  p <- x %>% 
    ncol()
  m <- y %>% 
    ncol()
  
  iter <- 0
  conv <- Inf
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  ll_store <- tibble(iter = 0, ll = -Inf)
  while(conv > 0 & iter < maxiter){
    iter <- iter + 1
    
    pi_vec <- fct_pi_vec(clust_assign, k, N)
    
    gamma <- fct_gamma(x, y, k, N, p, m, lam, rank, clust_assign, val_frac, penal_search)
    
    weighted_ll <- fct_weighted_ll(gamma)
    ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
    clust_assign_old <- clust_assign
    clust_assign <- fct_update_clust(N, gamma)
    clust_store <- clust_store %>% bind_rows(tibble(iter = rep(iter,N), assign=clust_assign))
    conv <- (clust_assign != clust_assign_old) %>% sum()
    # print(paste("iter", iter,"::: conv", conv,"::: lik",round(weighted_ll,digits = 2)))
    
  }
  return(list(ll = weighted_ll, assign = clust_assign, assign_store = clust_store, ll_store = ll_store, iter = iter))
}
