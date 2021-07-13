#' Main mixrrr function
#'
#' @param x matrix
#' @param y matrix
#' @param k number of clusters
#' @param lam penalization parameter
#' @param rank rank
#' @param chains number of chains
#' @param maxiter maximum number of iterations
#' @param verbose print progress?
#' @param val_frac validation fraction
#' @param penal_search penalization search parameter
#'
#' @return results
#' @export
#' 
#' @import dplyr
#' @import purrr
#' @import tictoc
#' @importFrom rlang .data
#'
#' @examples
hthmix <- function(x, y, k, lam = NULL, rank = NULL, chains = 50, maxiter = 100, verbose = TRUE, val_frac = 0.2, penal_search = 1:20/20){
  
  tictoc::tic()
  global_opt_ll <- -Inf
  lik_store <- NULL
  
  N <- x %>% nrow()
  clust_assign <- fct_initialize(k, N)
  
  result <- 1:chains %>% 
    purrr::map_dfr(
      .f = function(.c){
        
        clust_assign <- fct_pert_assign(clust_assign, lik_store, N, k)
        model_i <- fct_hthmix_comp(x, y, k, maxiter, lam, rank, clust_assign, val_frac, penal_search)
        ll_i <- model_i %>% purrr::pluck("ll")
        ll_store_i <- model_i %>% purrr::pluck("ll_store")
        clust_assign <- model_i %>% purrr::pluck("assign")
        A <- model_i %>% purrr::pluck("A")
        clust_assign_store <- model_i %>% purrr::pluck("assign_store")
        iter_i <- model_i %>% purrr::pluck("iter")
        lik_store <- c(lik_store, ll_i)
        global_opt_ll <- fct_global_opt(ll_i, global_opt_ll)

        return(tibble(llik = ll_i, chain = .c, iter = iter_i, assign = list(final_assign = clust_assign, A = A, assign_store = clust_assign_store, lik_store = ll_store_i)))
      }
    )
  time <- tictoc::toc()
  return(list(result = result, time = as.numeric(time$toc - time$tic)))
}
