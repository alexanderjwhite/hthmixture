fct_em <- function(x, y, k, lambda, clust_assign, lik_track, clust_store, em_iter){
  iter <- 0
  N <- nrow(x)
  changed <- Inf
  while(changed > 0 & iter < em_iter){
    iter <- iter + 1

    pi_vec <- fct_pi_vec(clust_assign, k, N)
    
    gamma_model <- fct_gamma(x = x, y = y, k = k, N = N, clust_assign = clust_assign, 
                             lambda = lambda, selection = "universal", alpha = 2*sqrt(3), 
                             beta = 1, y_sparse = TRUE, rank = NULL, max_rank = 3)
    gamma <- gamma_model$gamma
    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    lik_track <- rbind(lik_track,tibble(iter = lik_track$iter[nrow(lik_track)]+1,ll=weighted_ll,type="EM"))

    
    clust_assign_old <- clust_assign
    clust_assign <- fct_update_clust(gamma, N)
    print(clust_assign)
    clust_store <- clust_store %>% bind_rows(tibble(iter = rep((max(clust_store$iter)+1),N), assign=clust_assign))
    old_ll <- weighted_ll
    
    
    changed <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| changed: ", changed, "| ll", weighted_ll))
    lambda <- fct_select_lambda(x, y, k, clust_assign, initial = FALSE)
  }
  return(list(assign = clust_assign, lambda = lambda, lik_track = lik_track, 
              clust_store = clust_store, weighted_ll = weighted_ll))
  
}
