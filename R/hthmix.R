#' Main mixrrr function
#'
#' @param x matrix
#' @param y matrix
#' @param k_min minimum k to search
#' @param k_max maximum k to search
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
hthmix <- function(x, y, k_min = 2, k_max = 8, lam = NULL, rank = NULL, chains = 50, maxiter = 100, verbose = TRUE, val_frac = 0.2, penal_search = 1:20/20){
  
  tictoc::tic()
  global_opt_ll <- -Inf
  lik_store <- NULL
  
  # Substitute
  k <- k_min

  clust_assign <- fct_initialize(k, N)
  
  result <- 1:chains %>% 
    purrr::map_dfr(
      .f = function(.c){
        
        clust_assign <- fct_pert_assign(clust_assign, lik_store, N, k)
        model_i <- fct_hthmix_comp(x, y, k, maxiter, lam, rank, clust_assign, val_frac, penal_search)
        ll_i <- model_i %>% purrr::pluck("ll")
        ll_store_i <- model_i %>% purrr::pluck("ll_store")
        clust_assign <- model_i %>% purrr::pluck("assign")
        clust_assign_store <- model_i %>% purrr::pluck("assign_store")
        iter_i <- model_i %>% purrr::pluck("iter")
        lik_store <- c(lik_store, ll_i)
        global_opt_ll <- fct_global_opt(ll_i, global_opt_ll)
        # if(verbose){
        #   
        #   # print(paste("Chain:",.c,"ll:",scales::number(ll_i, accuracy = 0.01, big.mark = ",")))
        # }
        return(tibble(llik = ll_i, chain = .c, iter = iter_i, assign = list(final_assign = clust_assign, assign_store = clust_assign_store, lik_store = ll_store_i)))
      }
    )
  time <- tictoc::toc()
  return(list(result = result, time = as.numeric(time$toc - time$tic)))
}
