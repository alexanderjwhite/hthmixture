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
  lik_track <- tibble(iter=0, ll = -Inf, type = "EM")
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  ll_store <- tibble(iter = 0, ll = -Inf)
  skip_check <- 0
  old_ll <- -Inf
  
  store <- NULL
  for(j in 1:100){
    print(j)
    clust_assign <- fct_initialize(k,N)
    for(i in 1:k){
      x_k <- sim$x[which(clust_assign==i),]
      y_k <- sim$y[which(clust_assign==i),]
      sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
      store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
    }
  }
  lamb <- median(store)
  print(paste("lambda ############# ",lamb))
  
  while(conv > 0 & iter < maxiter){
    iter <- iter + 1
    # print(paste("Iteration",i,"..."))
    pi_vec <- fct_pi_vec(clust_assign, k, N)

    gamma_model <- fct_gamma(x, y, k, N, clust_assign, lamb, selection, alpha, beta, y_sparse, rank, max_rank)
    gamma <- gamma_model$gamma
    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    lik_track <- rbind(lik_track,tibble(iter = lik_track$iter[nrow(lik_track)]+1,ll=weighted_ll,type="EM"))
    ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
    
    clust_assign_old <- clust_assign
    clust_assign <- fct_update_clust(gamma, N)
    clust_store <- clust_store %>% bind_rows(tibble(iter = rep(iter,N), assign=clust_assign))
    old_ll <- weighted_ll


    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    
  }
  
  
  if(conv==0){
    skip_check <- skip_check + 1
    print("checking for final convergence...")
    clust_assign_old <- clust_assign
    new_change <- fct_sim_anneal(sim$x, sim$y, k, clust_assign, lamb, t_1 = 1000, mu = 0.95, eps = 1e-6, p = 0.9, N = 200,track = lik_track , max_iter = maxiter)
    clust_assign <- new_change$assign
    lik_track <- new_change$track
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    old_ll <- -new_change$lik
    store <- NULL
    for(j in 1:100){
      print(j)
      for(i in 1:k){
        x_k <- sim$x[which(clust_assign==i),]
        y_k <- sim$y[which(clust_assign==i),]
        sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
        store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
      }
    }
    lamb <- median(store)
    print("first check")
    print(paste("lambda ############# ",lamb))
  }
  

  
  
  while(conv > 0 & iter < maxiter & skip_check <= 3){
    iter <- iter + 1

    pi_vec <- fct_pi_vec(clust_assign, k, N)
    
    gamma_model <- fct_gamma(x, y, k, N, clust_assign, lamb, selection, alpha, beta, y_sparse, rank, max_rank)
    gamma <- gamma_model$gamma

    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    if(weighted_ll > old_ll){
      ll_store <- ll_store %>% bind_rows(tibble(iter = iter, ll = weighted_ll))
      lik_track <- rbind(lik_track,tibble(iter = lik_track$iter[nrow(lik_track)]+1,ll=weighted_ll,type="EM"))
      clust_assign_old <- clust_assign
      clust_assign <- fct_update_clust(gamma, N)
      clust_store <- clust_store %>% bind_rows(tibble(iter = rep(iter,N), assign=clust_assign))
      old_ll <- weighted_ll
    } else {
      clust_assign_old <- clust_assign
      # skip_check <- TRUE
    }

    
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    if(conv==0 & skip_check <=3){
      skip_check <- skip_check + 1
      print("checking for final convergence...")
      clust_assign_old <- clust_assign
      new_change <- fct_sim_anneal(sim$x, sim$y, k, clust_assign, lamb, t_1 = 1000, mu = 0.95, eps = 1e-6, p = 0.95, N = 200,track = lik_track , max_iter = maxiter)
      clust_assign <- new_change$assign
      lik_track <- new_change$track
      conv <- (clust_assign != clust_assign_old) %>% sum()
      print(paste("i: ",iter, "| conv: ", conv))
      store <- NULL
      for(j in 1:100){
        print(j)
        for(i in 1:k){
          x_k <- sim$x[which(clust_assign==i),]
          y_k <- sim$y[which(clust_assign==i),]
          sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
          store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
        }
      }
      lamb <- median(store)
      print(paste("check ####### ", skip_check))
      print(paste("lambda ############# ",lamb))
    }
    
  }
  
  return(list(ll = weighted_ll, assign = clust_assign, A = A, sig_vec = sig_vec, assign_store = clust_store, ll_store = lik_track, iter = iter))
}
