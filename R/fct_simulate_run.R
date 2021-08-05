#' Simulate and run mixrrr
#'
#' @param params list of parameters for simulation
#'
#' @return list of results
#' @export
#' 
#' @import purrr
#' @import dplyr
#' @import stats
#' @import tictoc
#' @import valse
#'
#' @examples 
fct_simulate_run <- function(params){
  
  maxiter <- params %>% purrr::pluck("maxiter")
  N <- params %>% purrr::pluck("N")
  k <- params %>% purrr::pluck("k")
  rho <- params %>% purrr::pluck("rho")
  sigma <- params %>% purrr::pluck("sigma")
  p <- params %>% purrr::pluck("p")
  m <- params %>% purrr::pluck("m")
  s <- params %>% purrr::pluck("s")
  r <- params %>% purrr::pluck("r")
  b <- params %>% purrr::pluck("b")
  
  prob <- rep(1/k,k)
  int <- prob %>% cumsum()
  rand_assign <- stats::runif(N)
  
  clust_assign_true <- (rand_assign) %>% 
    purrr::map_int(.f = function(.x){
      clust <- (.x <= int) %>% 
        which() %>% 
        min()
      return(clust)
    }) %>% 
    sort()
  
  clust_assign_true_key <- clust_assign_true %>% 
    dplyr::tibble() %>% 
    dplyr::mutate(order = 1:N) %>% 
    dplyr::arrange((.)) 
  
  clust_assign_true_vec <- clust_assign_true_key %>% 
    dplyr::pull(order)
  
  
  n <- clust_assign_true %>% 
    dplyr::as_tibble() %>% 
    dplyr::group_by(value) %>% 
    dplyr::summarize(n = dplyr::n()) %>% 
    dplyr::pull(n)
  
  clust_iter <- 1
  clust_min <- 1
  clust_max <- s
  x <- NULL
  y <- NULL
  for(i in 1:k){
    a_rows <- clust_min:clust_max
    clust_iter <- clust_iter + 1
    clust_min <- clust_max+1
    clust_max <- clust_iter*s
    
    sim <- fct_sim_mixrrr(n[i],a_rows,p,m,r[i],rho[i],sigma,b[i])
    x <- x %>% rbind(sim$X)
    y <- y %>% rbind(sim$Y)
  }
  x <- x[clust_assign_true_vec,]
  y <- y[clust_assign_true_vec,]

  
  print("starting hth")
  model_hth <- hthmix(x, y, k, rank = r[1], chains = 1, maxiter = maxiter, penal_search = 1:100/100)
  final_assign_hth <- model_hth %>% purrr::pluck("result","assign","final_assign")
  iter_hth <- model_hth %>% purrr::pluck("result", "iter", 1)
  time_hth <- model_hth %>% purrr::pluck("time")
  print("hth finished")
  
  print("starting kmeans")
  data_kmeans <- x %>% cbind(y)
  tictoc::tic()
  model_kmeans <- kmeans(data_kmeans, k)
  time_kmeans <- tictoc::toc()
  time_kmeans <- time_kmeans$toc-time_kmeans$tic
  print("kmeans finished")
  
  print("starting valse")
  safe_valse <- purrr::safely(valse::runValse)
  tictoc::tic()
  model_valse <- safe_valse(x, y, kmin = k, kmax = k, rank.min = r[1], rank.max = r[1], verbose = TRUE)
  time_valse <- tictoc::toc()
  time_valse <- time_valse$toc-time_valse$tic
  print("valse finished")
  
  return(list(true=clust_assign_true,
              est_hth=final_assign_hth, 
              iter_hth = iter_hth, 
              time_hth = time_hth,
              est_kmeans = model_kmeans,
              time_kmeans = time_kmeans,
              est_valse = model_valse,
              time_valse = time_valse
              ))
}
