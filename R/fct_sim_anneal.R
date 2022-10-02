fct_j_lik <- function(x, y, k, clust_assign, lambda,  selection = "universal", alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL){
  N <- nrow(x)
  gamma_model <- fct_gamma(x, y, k, N, clust_assign, lambda, selection = "universal", alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL)
  gamma <- gamma_model$gamma
  return(-fct_weighted_ll(gamma))
}

fct_new_assign <- function(assign, k, p){
  n <- length(assign)
  new_assign <- assign
  i <- 1
  flag <- FALSE

  while(!(i == n & flag)){
    if(i > n){i <- 1}
    u <- runif(1)
    if(u > p){
      flag <- TRUE
      cur_assign <- assign[i]
      new_assign[i] <- sample((1:k)[-which(cur_assign==1:k)],1)
    }
    i <- i + 1
  }
  return(new_assign)
}



fct_sim_anneal <- function(x, y, k, init_assign, lambda, t_1, mu, eps, p, N, track, clust_store, anneal_iter = 1e3){
  total_iter <- 0
  count <- 0

  a_t <- a_b <- a_c <- init_assign
  j_b <- j_c <- fct_j_lik(x, y, k, init_assign, lambda)
  t <- t_1

  # step 1 - trial assignment


  while(t >= eps & total_iter <= anneal_iter){
    total_iter <- total_iter + 1
    a_t <- fct_new_assign(a_b, k, p) #step 1
    j_t <- fct_j_lik(x, y, k, a_t, lambda)

    if(j_t <= j_c){ #step 2
      a_c <- a_t
      j_c <- j_t

      if(j_t >= j_b){
        count <- count + 1
      } else {
        changed <- sum(a_b != a_t)
        clust_store <- clust_store %>% bind_rows(tibble(iter = rep(total_iter,nrow(x)), assign=a_b))
        j_b <- j_t
        a_b <- a_t
        count <- 0
        total_iter <- 0
        
        cat(paste("\n i: ",total_iter, "| changed: ", changed, "| ll", -j_b))
      }
    } else {
      u <- runif(1) #step 3
      b <- exp(-(j_t - j_c)/t)
      # print(b)
      if(u <= b){
        j_c <- j_t
        a_c <- a_t
      }
    }
    if(count >= N){ #step 4
      t <- mu*t
    }
    track <- rbind(track,tibble(iter=track$iter[nrow(track)]+1, ll=-j_b, type = "sim"))
    # print(paste(total_iter, "|", t,"|",count, "|", j_b, "|", j_c,"|",j_t))
    if(total_iter %% 100 == 0){cat(".")}
  }
  cat("\n")
  
  lambda <- fct_select_lambda(x, y, k, a_b, initial = FALSE)
  print(paste("sim anneal ll",-fct_j_lik(x, y, k, a_b, lambda)))
  print(a_b)
  return(list(assign = a_b, lik_track = track, weighted_ll = -j_b, lambda = lambda, clust_store = clust_store))
}
