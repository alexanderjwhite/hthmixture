simulate_hthmix <- function(params){
  require(dplyr)
  num_pert <- params %>% purrr::pluck("num_pert")
  chains <- params %>% purrr::pluck("chains")
  lam <- params %>% purrr::pluck("lam")
  maxiter <- params %>% purrr::pluck("maxiter")
  N <- params %>% purrr::pluck("N")
  prob <- params %>% purrr::pluck("prob")
  k <- params %>% purrr::pluck("k")
  nvld <- params %>% purrr::pluck("nvld")
  rho <- params %>% purrr::pluck("rho")
  sigma <- params %>% purrr::pluck("sigma")
  p <- params %>% purrr::pluck("p")
  m <- params %>% purrr::pluck("m")
  s <- params %>% purrr::pluck("s")
  r <- params %>% purrr::pluck("r")
  b <- params %>% purrr::pluck("b")
  
  num_pert <- (num_pert*N) %>% floor()
  int <- prob %>% cumsum()
  rand_assign <- runif(N)
  names <- paste0("c_",1:k)
  if(k*s > p){print("FOR SEPARATION, VERIFY THAT K*S < P")}
  
  # make sure k*s is less than p
  
  clust_assign_true <- (rand_assign) %>% 
    purrr::map_int(.f = function(.x){
      clust <- (.x <= int) %>% 
        which() %>% 
        min()
      return(clust)
    }) %>% 
    sort()
  
  clust_assign_true_key <- clust_assign_true %>% 
    tibble() %>% 
    mutate(order = 1:N) %>% 
    arrange((.)) 
  
  clust_assign_true_vec <- clust_assign_true_key %>% 
    pull(order)
  
  
  n <- clust_assign_true %>% 
    as_tibble() %>% 
    group_by(value) %>% 
    summarize(n = n()) %>% 
    pull(n)
  
  gamma_store <- as_tibble(matrix(rep(0,k),nrow=1)) %>% 
    rename_if(is.numeric,~names) %>% 
    mutate(w_ll = 0, iter = 0)
  
  # chain_store <- tibble(obs = 0, assign = 0, llik = 0, chain = 0)
  clust_iter <<- 1
  clust_min <<- 1
  clust_max <<- s
  data_k <- n %>% 
    list(r,rho,b) %>% 
    purrr::pmap(.f = function(.n,.r,.rho,.b){
      a_rows <- clust_min:clust_max
      clust_iter <<- clust_iter + 1
      clust_min <<- clust_max+1
      clust_max <<- clust_iter*s
      
      sim <- sim_sarrs(.n,nvld,.n,a_rows,p,m,.r,.rho,sigma,.b)
      A <- sim %>% 
        purrr::pluck("A")
      S <- sim %>% 
        purrr::pluck("sig")
      X <- sim %>% 
        purrr::pluck("X")
      Y <- sim %>% 
        purrr::pluck("Y")
      return(list(X=X,Y=Y,A=A,S=S))
    })
  
  stack_X <- data_k %>% 
    purrr::map_dfr(.f = function(.x){
      X <- .x %>% 
        purrr::pluck("X")
      return(as_tibble(X))
    }) %>% 
    slice(clust_assign_true_vec) %>% 
    as.matrix()
  
  stack_Y <- data_k %>% 
    purrr::map_dfr(.f = function(.x){
      Y <- .x %>% 
        purrr::pluck("Y")
      return(as_tibble(Y))
    }) %>% 
    slice(clust_assign_true_vec) %>% 
    as.matrix()
  
  
  init_int <- rep(1/k,k) %>% cumsum()
  init_rand_assign <- runif(N)
  
  clust_assign <<- (init_rand_assign) %>%
    purrr::map_int(.f = function(.x){
      clust <- (.x <= init_int) %>%
        which() %>%
        min()
      return(clust)
    })
  
  change <- 0
  weighted_g_ll_min <<- 0
  clust_assign_g <<- clust_assign
  chain_clust <- 1:chains %>% 
    purrr::map_dfr(
      .f = function(.c){
        
        if(abs(.c - change) > 20){
          init_int <- rep(1/k,k) %>% cumsum()
          init_rand_assign <- runif(N)
          
          clust_assign <<- (init_rand_assign) %>%
            purrr::map_int(.f = function(.x){
              clust <- (.x <= init_int) %>%
                which() %>%
                min()
              return(clust)
            })
          clust_assign_g <<- clust_assign
          
        } else {
          pert_samples <- sample(1:N, num_pert)
          unique_vals <- 1:k
          for(i in 1:length(pert_samples)){
            obs <- pert_samples[i]
            orig <- clust_assign_g[obs]
            vals <- unique_vals[unique_vals != orig]
            if(length(vals) == 1){
              clust_assign_g[obs] <- vals
            } else {
              clust_assign_g[obs] <- sample(vals)
            }
          }
        }
        
        
        
        clust_assign <<- clust_assign_g
        
        conv <- Inf
        main_clust <- Inf
        iter <- 0
        while(conv>0 & iter < maxiter){
          iter <- iter+1
          
          pi_vec <- 1:k %>% 
            purrr::map_dbl(
              .f = function(.x){
                sum(clust_assign==.x)/N
              }
            )
          
          
          gamma <- 1:k %>% 
            list() %>% 
            purrr::pmap_dfc(.f = function(.x){
              
              select <- (clust_assign==.x) %>% 
                which()
              
              nk <- select %>% length()
              
              X_k <- stack_X %>% 
                as_tibble() %>% 
                slice(select) %>% 
                as.matrix()
              
              Y_k <- stack_Y %>% 
                as_tibble() %>% 
                slice(select) %>% 
                as.matrix()
              
              # P=X%*% ginv(t(X)%*%X) %*%t(X)
              # rhat = sum(svd(P%*%Y)$d > sigmahat*(sqrt(2*m)+sqrt(2*min(n,p))))
              
              rhat <- r[.x]
              estimate <- SARRS(Y_k,X_k,rhat, lam, "grLasso")
              
              A_k <- estimate %>% 
                purrr::pluck("Ahat")
              
              sig_vec <- estimate %>% 
                purrr::pluck("sigvec")
              
              mu_mat <- (stack_X %>% 
                           bind_cols(int = rep(1,N)) %>% 
                           as.matrix()) %*% A_k
              
              gam <- array(dim = N)
              for(i in 1:N){
                mu_i <- mu_mat[i,]
                gam[i] <- 1:m %>% 
                  purrr::map_dbl(
                    .f = function(.y){
                      dnorm(stack_Y[i,.y], mean = mu_i[.y], sd = sig_vec[.y])
                    }
                  ) %>%
                  log() %>% 
                  sum() 
              } 
              return(tibble(gam) %>% setNames(names[.x]))
              
            }) 
          
          weighted_ll <- gamma %>% 
            rowwise() %>% 
            mutate(ll = max(c_across(where(is.numeric)))) %>% 
            pull(ll)
          
          
          clust_assign_old <- clust_assign
          clust_assign <<- 1:N %>% 
            purrr::map_int(
              .f = function(.x){
                res <- (gamma %>% as.matrix())[.x,] %>% which.max()
              }
            )
          
          
          
          conv <- (clust_assign!=clust_assign_old) %>% sum()
          
        }
        
        weighted_g_ll <<- weighted_ll %>% sum()
        
        if(weighted_g_ll > weighted_g_ll_min){
          weighted_g_ll_min <<- weighted_g_ll
          clust_assign_g <<- clust_assign
          print(weighted_g_ll)
          change <<- .c
          
        }
        
        shuffle <- clue::solve_LSAP(table(clust_assign_true, clust_assign), maximum = TRUE)
        
        acctbl <- (table(clust_assign_true, clust_assign)[,shuffle])
        acc <- (acctbl %>% diag() %>% sum()) / N
        return(tibble(llik = weighted_g_ll, chain = .c, acc = acc))
      }
    )
  
  best_chain <- chain_clust %>% 
    unique() %>% 
    arrange(desc(llik)) %>% 
    slice(1:3) %>% 
    mutate(top = 1:3) %>% 
    tidyr::pivot_wider(names_from = top, values_from = c("llik", "chain", "acc"))
  
  return(best_chain)
}






