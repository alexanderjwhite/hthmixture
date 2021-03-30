paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/Carbonate/r_404_packages/",paths)
.libPaths(paths)
library(dplyr)
future::plan("multiprocess")
source("/geode2/home/u100/whitealj/BigRed3/scripts/20210205_sarrs_alex.R")
source("/geode2/home/u100/whitealj/BigRed3/scripts/20210304_large_sim_alex.R")
set.seed(19921124)

##################### FULL PARAM ##############################

# num_pert <- c(0.05,0.1,0.2)
# N <- c(50, 100, 200, 400)
# PROB <- 0.5
# k <- 2:8
# sigma <- c(0.5,1,2)
# p <- c(50, 100, 1000, 10000)
# m <- c(50, 100, 1000, 10000)
# s <- c(5, 10, 20)
# R <- 1:5 #up to 3
# B <- c(1, 5, 10)
# replication <- 1:2

###############################################################

##################### TRUNCATED PARAM ##############################

num_pert <- c(0.2)
N <- c(50, 100, 200, 400)
PROB <- 0.5
k <- 2
sigma <- c(1)
# p <- c(50, 100, 1000, 10000)
# m <- c(50, 100, 1000, 10000)
dim <- c(50, 100, 1000)
s <- c(5, 10, 20)
R <- 1:2 #up to 3
B <- c(1)
replication <- 1

###############################################################

safe_result <- purrr::safely(simulate_hthmix)

sim_param <- expand.grid(
  num_pert = num_pert,
  N = N,
  PROB = PROB,
  k = k,
  sigma = sigma,
  # p = p,
  # m = m,
  dim = dim,
  s = s,
  R = R,
  B = B, 
  replication = replication
  )

tictoc::tic()
result_final <- sim_param %>% 
  # filter(k==8 & num_pert == 0.2 & N == 50 & sigma == 1 & p == 1000 & m == 1000 & R == 1 & s == 5 & B == 1 & replication == 1) %>% 
  # sample_n(30) %>%
  # slice(c(1)) %>%
  as.list() %>% 
  furrr::future_pmap_dfr(
  # purrr::pmap_dfr(
  .options = furrr::future_options(seed = TRUE),
  .progress = TRUE,
  .f = function( num_pert,
                 N,
                 PROB,
                 k,
                 sigma,
                 dim,
                 s,
                 R,
                 B,
                 replication){
    
    
    chains <- 3
    lam <- 1
    maxiter <- 100
    nvld <- 1e4
    rho <- rep(0,k)
    b <- rep(B, k)
    r <- rep(R, k)
    prob <- rep(1/k, k)
    p <- m <- dim
    
    params <- list(
      chains = chains,
      lam = lam,
      maxiter = maxiter,
      nvld = nvld,
      num_pert = num_pert,
      N = N,
      prob = prob,
      k = k,
      sigma = sigma,
      p = p,
      rho = rho,
      m = m,
      s = s,
      r = r,
      b = b)
    
    print(params)
    
    sim_result <- safe_result(params)
    
    if(is.null(sim_result$result)){
      result <- tibble(
          num_pert,
          N,
          PROB,
          k,
          sigma,
          p,
          m,
          s,
          R,
          B,
          replication,
          error = paste(sim_result$error),
          llik_1 = NA,
          llik_2 = NA,
          llik_3 = NA,
          chain_1 = NA,
          chain_2 = NA,
          chain_3 = NA,
          acc_1 = NA,
          acc_2 = NA,
          acc_3 = NA
        ) 
    } else {
      result <- tibble(
        num_pert,
        N,
        PROB,
        k,
        sigma,
        p,
        m,
        s,
        R,
        B,
        replication,
        error = NULL,
      ) %>%
        bind_cols(sim_result)
    }
    
    return(result)
  }
)
tictoc::toc()

result_final %>% 
  readr::write_rds("20210304_large_sim_results.rds")
