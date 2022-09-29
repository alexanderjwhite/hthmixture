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
fct_hthmix_comp <- function(x, y, k, maxiter, clust_assign, selection, 
                            alpha, beta, y_sparse, rank, max_rank,
                            temp, p1, p2, sim_N, checks, true){

  N <- nrow(x)
  p <- ncol(x)
  m <- ncol(y)
  safe_rank <- purrr::safely(fct_rank)
  eta_k <- 3
  iter <- 0
  conv <- Inf
  lik_track <- tibble(iter=0, ll = -Inf, type = "EM")
  clust_store <- tibble(iter=rep(iter,N),assign=clust_assign)
  ll_store <- tibble(iter = 0, ll = -Inf)
  skip_check <- 0
  old_ll <- -Inf
  lamb <- rep(0,k)
  
  
  store <- array(0, dim = c(100,k,2))
  for(j in 1:100){
    print(j)
    clust_assign <- fct_initialize(k,N)
    for(i in 1:k){
      x_k <- sim$x[which(clust_assign==i),]
      y_k <- sim$y[which(clust_assign==i),]
      sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
      
      rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta_k)
      rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
      rank_hat <- min(rank_hat, max_rank)
      
      store[j,i,] <- fct_sarrs(y_k, x_k, r = rank_hat, lam = NULL, 
                               alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
                               ptype = "grLasso", y_sparse = TRUE)$lambda_store
      # store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
    }
  }
  store_mat <- rbind(store[,,1],store[,,2])
  # lamb <- median(store)
  lamb <- apply(store_mat, 2, median)
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
    
    store <- array(0, dim = c(50,k,2))
    for(j in 1:50){
      print(j)
      for(i in 1:k){
        x_k <- sim$x[which(clust_assign==i),]
        y_k <- sim$y[which(clust_assign==i),]
        sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
        
        rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta_k)
        rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
        rank_hat <- min(rank_hat, max_rank)
        
        store[j,i,] <- fct_sarrs(y_k, x_k, r = rank_hat, lam = NULL, 
                                 alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
                                 ptype = "grLasso", y_sparse = TRUE)$lambda_store
        # store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
      }
    }
    store_mat <- rbind(store[,,1],store[,,2])
    # lamb <- median(store)
    lamb <- apply(store_mat, 2, median)
    print("first check")
    print(paste("lambda ############# ",lamb))
    
    clust_reorder <- clue::solve_LSAP(table(true, clust_assign), maximum = TRUE)[clust_assign]
    ari <- mclust::adjustedRandIndex(true, clust_reorder)
    print(paste("init ari ###########",ari))
    
    
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
    
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))
    
  }
  
  



  if(conv==0){
    skip_check <- skip_check + 1
    print("checking for final convergence...")
    clust_assign_old <- clust_assign
    new_change <- fct_sim_anneal(sim$x, sim$y, k, clust_assign, lamb, t_1 = temp, mu = 0.95, eps = 1e-6, p = p1, N = sim_N,track = lik_track , max_iter = maxiter)
    clust_assign <- new_change$assign
    lik_track <- new_change$track
    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))

    store <- array(0, dim = c(50,k,2))
    for(j in 1:50){
      print(j)
      for(i in 1:k){
        x_k <- sim$x[which(clust_assign==i),]
        y_k <- sim$y[which(clust_assign==i),]
        sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
        rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta_k)
        rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
        rank_hat <- min(rank_hat, max_rank)
        
        store[j,i,] <- fct_sarrs(y_k, x_k, r = rank_hat, lam = NULL,
                                 alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
                                 ptype = "grLasso", y_sparse = TRUE)$lambda_store
        # store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
      }
    }
    store_mat <- rbind(store[,,1],store[,,2])
    # lamb <- median(store)
    lamb <- apply(store_mat, 2, median)
    print("first check")
    print(paste("lambda ############# ",lamb))
    old_ll <- -fct_j_lik(x,y,k,clust_assign,lamb)
  }
  conv <- 1



  while(conv > 0 & iter < maxiter & skip_check <= checks){
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
    if(conv==0 & skip_check <=checks){
      skip_check <- skip_check + 1
      print("checking for final convergence...")
      clust_assign_old <- clust_assign
      new_change <- fct_sim_anneal(sim$x, sim$y, k, clust_assign, lamb, t_1 = temp, mu = 0.95, eps = 1e-6, p = p2, N = sim_N,track = lik_track , max_iter = maxiter)
      clust_assign <- new_change$assign
      lik_track <- new_change$track
      conv <- (clust_assign != clust_assign_old) %>% sum()
      print(paste("i: ",iter, "| conv: ", conv))
      store <- array(0, dim = c(50,k,2))
      for(j in 1:50){
        print(j)
        for(i in 1:k){
          x_k <- sim$x[which(clust_assign==i),]
          y_k <- sim$y[which(clust_assign==i),]
          sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
          
          rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta_k)
          rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
          rank_hat <- min(rank_hat, max_rank)
          store[j,i,] <- fct_sarrs(y_k, x_k, r = rank_hat, lam = NULL,
                                   alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
                                   ptype = "grLasso", y_sparse = TRUE)$lambda_store
          # store <- c(store,fct_sarrs(y_k, x_k, 1, NULL, sigma = sigma_hat)$lambda_store)
        }
      }
      store_mat <- rbind(store[,,1],store[,,2])
      # lamb <- median(store)
      lamb <- apply(store_mat, 2, median)
      print(paste("check ####### ", skip_check))
      print(paste("lambda ############# ",lamb))
      old_ll <- -fct_j_lik(x,y,k,clust_assign,lamb)
    }
    clust_reorder <- clue::solve_LSAP(table(true, clust_assign), maximum = TRUE)[clust_assign]
    ari <- mclust::adjustedRandIndex(true, clust_reorder)
    print(paste("########### loop ari ###########",ari))


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

    conv <- (clust_assign != clust_assign_old) %>% sum()
    print(paste("i: ",iter, "| conv: ", conv))

  }

  return(list(ll = weighted_ll, assign = clust_assign, A = A, sig_vec = sig_vec, assign_store = clust_store, ll_store = lik_track, iter = iter))
}
