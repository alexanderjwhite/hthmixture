#' Compute mixrrr
#'
#' @inheritParams hthmix
#' @param clust_assign initial cluster assignment
#'
#' @return results
#' @export
#' 
#' @import dplyr
#'
#' @examples
fct_hthmix_comp <- function(x, y, k, maxiter, clust_assign, selection, alpha, beta, y_sparse, rank){
  
  N <- nrow(x)
  p <- ncol(x)
  m <- ncol(y)
  
  iter <- 0
  conv <- Inf
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  ll_store <- tibble(iter = 0, ll = -Inf)
  while(conv > 0 & iter < maxiter){
    iter <- iter + 1
    # print(paste("Iteration",i,"..."))
    pi_vec <- fct_pi_vec(clust_assign, k, N)
    
    gamma_model <- fct_gamma(x, y, k, N, clust_assign, selection, alpha, beta, y_sparse, rank)
    gamma <- gamma_model$gamma
    # print(gamma)
    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
    clust_assign_old <- clust_assign
    clust_assign <- fct_update_clust(gamma, N)
    clust_store <- clust_store %>% bind_rows(tibble(iter = rep(iter,N), assign=clust_assign))
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    
  }
  return(list(ll = weighted_ll, assign = clust_assign, A = A, sig_vec = sig_vec, assign_store = clust_store, ll_store = ll_store, iter = iter))
}
# NMF::aheatmap(Ahat,Rowv = NA, Colv = NA)
# clust_assign_1 <- clust_assign
# errors_1 <- errors
# grid_1 <- grid
# lam_grid_1 <- lam_grid 
# errors_2 <- errors
# grid_2 <- grid
# lam_grid_2 <- lam_grid
