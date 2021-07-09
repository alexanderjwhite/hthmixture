# num_pert <- 0.2
# N <- 50
# PROB <- 0.5
# k <- 2
# sigma <- 1
# dim <- 50
# s <- 5
# R <- 1
# B <- 1
# replication <- 1
# 
# chains <- 500
# lam <- 1
# maxiter <- 100
# nvld <- 1e4
# rho <- rep(0,k)
# b <- rep(B, k)
# r <- rep(R, k)
# prob <- rep(1/k, k)
# p <- m <- dim
# 
# params <- list(
#   chains = chains,
#   lam = lam,
#   maxiter = maxiter,
#   nvld = nvld,
#   num_pert = num_pert,
#   N = N,
#   prob = prob,
#   k = k,
#   sigma = sigma,
#   p = p,
#   rho = rho,
#   m = m,
#   s = s,
#   r = r,
#   b = b)

simulate_kmeans <- function(params){
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
  
  data_kmeans <- stack_Y %>% cbind(stack_X)
  
  k_result <- kmeans(data_kmeans, 2)
  clust_assign <- k_result$cluster
  
  shuffle <- clue::solve_LSAP(table(clust_assign_true, clust_assign), maximum = TRUE)
  
  acctbl <- (table(clust_assign_true, clust_assign)[,shuffle])
  acc <- (acctbl %>% diag() %>% sum()) / N

  return(acc)
}






