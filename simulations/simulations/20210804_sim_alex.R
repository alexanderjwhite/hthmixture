paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/BigRed3/r_4_0_4_library/",paths)
.libPaths(paths)
library(hthmixture)
N_lab <- "N"
k_lab <- "k"
sigma_lab <- "sigma"
dim_lab <- "dim"
s_lab <- "s"
r_lab <- "r"
names(N_lab) <- "N"
names(k_lab) <- "k"
names(sigma_lab) <- "sigma"
names(dim_lab) <- "dim"
names(s_lab) <- "s"
names(r_lab) <- "r"

N_param <- dials::new_qual_param(
  type = "character",
  values = c("200","400"),
  label = N_lab
)

k_param <- dials::new_qual_param(
  type = "character",
  values = c("2","3"),
  label = k_lab
)

sigma_param <- dials::new_qual_param(
  type = "character",
  values = c("1","small", "large"),
  label = sigma_lab
)

dim_param <- dials::new_qual_param(
  type = "character",
  values = c("small","large"),
  label = dim_lab
)

s_param <- dials::new_qual_param(
  type = "character",
  values = c("5","10"),
  label = s_lab
)

r_param <- dials::new_qual_param(
  type = "character",
  values = c("1","2"),
  label = r_lab
)

param_dials <- list(N_param, k_param, sigma_param, dim_param, s_param, r_param)
params <- dials::grid_latin_hypercube(param_dials, size = 15) %>% 
  mutate(N = as.numeric(N),
         k = as.numeric(k),
         s = as.numeric(s),
         r = as.numeric(r))

grid <- params %>% 
  slice(rep(1:nrow(params),20)) %>% 
  mutate(rep = rep(1:nrow(params),20))

future::plan("multicore")

safe_run <- purrr::safely(fct_simulate_run)

result <- list(grid$N, grid$k, grid$sigma, grid$dim, grid$s, grid$r, grid$rep,1:nrow(grid)) %>%
  furrr::future_pmap(.options = furrr::future_options(seed = TRUE),
  # purrr::pmap(
                     .f = function(.N, .k, .sigma, .dim, .s, .r, .rep, .prog){
                       
                       if(.dim == "small"){
                         p <- 500
                         m <- 50
                       } else {
                         p <- 1e3
                         m <- 200
                       }
                       
                       params <- list(
                         maxiter = 1e3,
                         N = .N,
                         k = .k,
                         sigma = .sigma,
                         p = p,
                         rho = rep(0,.k),
                         m = m,
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

result %>% saveRDS("20210804_results.rds")