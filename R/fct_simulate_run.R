#' Simulate and run mixrrr
#'
#' @param param list of parameters for simulation
#'
#' @return list of results
#' @export
#' 
#' @import purrr
#' @import dplyr
#' @import stats
#'
#' @examples
fct_simulate_run <- function(param){
  
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
      
      sim <- fct_sim_mixrrr(.n,a_rows,p,m,.r,.rho,sigma,.b)
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
  
  x <- data_k %>% 
    purrr::map_dfr(.f = function(.x){
      X <- .x %>% 
        purrr::pluck("X")
      return(as_tibble(X))
    }) %>% 
    dplyr::slice(clust_assign_true_vec) %>% 
    as.matrix()
  
  y <- data_k %>% 
    purrr::map_dfr(.f = function(.x){
      Y <- .x %>% 
        purrr::pluck("Y")
      return(as_tibble(Y))
    }) %>% 
    dplyr::slice(clust_assign_true_vec) %>% 
    as.matrix()
  
  model <- hthmix(x, y, chains = 1, maxiter = maxiter)
  
  final_assign <- model %>% purrr::pluck("result","assign","final_assign")
  
  iter <- model %>% purrr::pluck("result", "iter", 1)
  
  time <- model %>% purrr::pluck("time")
  
  return(list(true=clust_assign_true,est=final_assign, iter = iter, time = time))
}
