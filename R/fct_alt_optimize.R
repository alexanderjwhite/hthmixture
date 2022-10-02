fct_alt_optimize <- function(x, y, k, clust_assign, lambda, alt_iter, anneal_iter, em_iter, temp, mu, eps, accept_prob, sim_N){
  N <- nrow(x)
  p <- ncol(x)
  m <- ncol(y)

  if(is.null(clust_assign)){
    clust_assign <- fct_initialize(k,N)
  }
  iter <- 0
  lik_track <- tibble(iter=0, ll = -Inf, type = "EM")
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  changed <- Inf
  
  # Select Initial Lambda
  if(is.null(lambda)){
    lambda <- fct_select_lambda(x, y, k, clust_assign = NULL, initial = TRUE)
  }
  
  
  
  
  while(changed > 0 & iter < alt_iter){
    
    # EM 
    clust_assign_old <- clust_assign
    model_em <- fct_em(x = x, y = y, k = k, lambda = lambda, clust_assign = clust_assign, 
                       lik_track = lik_track, clust_store = clust_store, em_iter = em_iter)
    clust_assign <- model_em$assign
    lambda <- model_em$lambda
    lik_track <- model_em$lik_track
    clust_store <- model_em$clust_store
    weighted_ll <- model_em$weighted_ll
    
    
    # Simulated Annealing
    
    model_anneal <- fct_sim_anneal(x = x, y = y, k, init_assign = clust_assign, 
                                   lambda = lambda, t_1 = temp, mu = mu, eps = eps, 
                                   p = accept_prob, N = sim_N, track = lik_track, 
                                   clust_store = clust_store, anneal_iter = anneal_iter)
    clust_assign <- model_anneal$assign
    lambda <- model_anneal$lambda
    lik_track <- model_anneal$lik_track
    clust_store <- model_anneal$clust_store
    weighted_ll <- model_anneal$weighted_ll
    
    changed <- sum(clust_assign != clust_assign_old)
    iter <- iter + 1
    print(paste0("Full Cycle | ", iter, " Changed | ",changed))
    
    
  }
  
  # final fit
  
  final_model <- fct_gamma(x, y, k, N, clust_assign, lambda, selection = "universal", 
                           alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, 
                           max_rank = 3, rank = NULL)
  gamma <- final_model$gamma
  A <- final_model$A
  sig_vec <- final_model$sig_vec
  weighted_ll <- fct_weighted_ll(gamma)
  
  return(list(ll = weighted_ll, assign = clust_assign, A = A, sig_vec = sig_vec, 
              assign_store = clust_store, ll_store = lik_track, iter = iter))

}