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
fct_hthmix_comp <- function(x, y, k, maxiter, clust_assign, selection, alpha, beta, y_sparse, rank, max_rank){
  
  N <- nrow(x)
  p <- ncol(x)
  m <- ncol(y)
  
  iter <- 0
  conv <- Inf
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  ll_store <- tibble(iter = 0, ll = -Inf)
  skip_check <- FALSE
  old_ll <- -Inf
  while(conv > 0 & iter < maxiter){
    iter <- iter + 1
    # print(paste("Iteration",i,"..."))
    pi_vec <- fct_pi_vec(clust_assign, k, N)
    
    gamma_model <- fct_gamma(x, y, k, N, clust_assign, selection, alpha, beta, y_sparse, rank, max_rank)
    gamma <- gamma_model$gamma
    # print(gamma)
    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    if(weighted_ll > old_ll){
      ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
      clust_assign_old <- clust_assign
      clust_assign <- fct_update_clust(gamma, N)
      clust_store <- clust_store %>% bind_rows(tibble(iter = rep(iter,N), assign=clust_assign))
    } else {
      clust_assign_old <- clust_assign
      skip_check <- TRUE
    }
    old_ll <- weighted_ll
    # ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
    # clust_assign_old <- clust_assign
    # clust_assign <- fct_update_clust(gamma, N)
    
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    if(conv==0 & !skip_check){
      
      print("checking for final convergence...")
      clust_assign_old <- clust_assign
      new_change <- fct_conv_check(x, y, k, N, clust_assign, selection = "universal", alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL,ll = weighted_ll)
      clust_assign <- new_change$clust_assign
      conv <- (clust_assign != clust_assign_old) %>% sum()
      print(paste("i: ",iter, "| conv: ", conv))
    }
  }
  return(list(ll = weighted_ll, assign = clust_assign, A = A, sig_vec = sig_vec, assign_store = clust_store, ll_store = ll_store, iter = iter))
}
