library(dplyr)
library(hthmixture)
N <- c(100,200)
k <- 2:3
sigma <- 1:2
dim <- c(50,200)
s <- c(10,20)
r <- 1:2
rep <- 1:10

future::plan("multiprocess")
grid <- (expand.grid(N=N,k=k,sigma=sigma,dim=dim,s=s,r=r,rep=rep))[1:2,]

result <- list(grid$N, grid$k, grid$sigma, grid$dim, grid$s, grid$r, grid$rep,1:nrow(grid)) %>% 
  furrr::future_pmap(.f = function(.N, .k, .sigma, .dim, .s, .r, .rep, .prog){
    params <- list(
      maxiter = 1e3,
      N = .N,
      k = .k,
      sigma = .sigma,
      p = .dim,
      rho = rep(0,.k),
      m = .dim,
      s = .s,
      r = rep(.r,.k),
      b = (1:.k)*5)
    fct_simulate_run(params)
    fct_push_me(round(.prog/nrow(grid), digits = 2))
  })

result %>% saveRDS("20210709_results.rds")





