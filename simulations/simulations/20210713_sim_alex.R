paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/BigRed3/r_4_0_4_library/",paths)
.libPaths(paths)
library(hthmixture)
N <- c(100,200)
k <- 3
sigma <- 1
dim <- c(50,200)
s <- c(5,15)
r <- 1:2
rep <- 1:3

future::plan("multicore")
grid <- (expand.grid(N=N,k=k,sigma=sigma,dim=dim,s=s,r=r,rep=rep))

safe_run <- purrr::safely(fct_simulate_run)

result <- list(grid$N, grid$k, grid$sigma, grid$dim, grid$s, grid$r, grid$rep,1:nrow(grid)) %>%
  furrr::future_pmap(.options = furrr::future_options(seed = TRUE),
                     .f = function(.N, .k, .sigma, .dim, .s, .r, .rep, .prog){
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
                       if((.prog %% 10) == 0){
                         fct_push_me(round(.prog/nrow(grid), digits = 2))
                       }
                       
                       results <- safe_run(params)
                       
                       if(is.null(results$error)){
                         return(list(N = .N, k = .k, sigma = .sigma, dim = .dim, s = .s, r = .r, rep = .rep, id = .prog, res = results$result))
                       } else {
                         return(list(N = .N, k = .k, sigma = .sigma, dim = .dim, s = .s, r = .r, rep = .rep, id = .prog, res = results$error))
                       }
                       
                       
                     })

result %>% saveRDS("20210722_results.rds")