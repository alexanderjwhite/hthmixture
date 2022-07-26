#' Compute Model Error
#'
#' @param x matrix
#' @param y matrix
#' @param k number of clusters
#' @param model model object from hthmix
#'
#' @return error
#' @import purrr
#' @import dplyr
#' @export
#'
fct_full_model_error <- function(x, y, k, model){
  
  N <- x %>% nrow()
  
  m <- y %>% ncol()
  
  A <- model %>% 
    purrr::pluck("result", "assign", "A")
  
  sig_vec <- model %>% 
    purrr::pluck("result", "assign", "sig_vec")
  
  est_assign <- 1:k %>% 
    list() %>% 
    purrr::pmap_dfc(.f = function(.x){
      
      A_k <- A %>% purrr::pluck(.x)
      sig_k <- sig_vec %>% purrr::pluck(.x)
      
      if(is.null(A_k)|is.null(sig_k)){return(dplyr::tibble(-Inf))}
      
      mu_mat <- (x %>% 
                   dplyr::bind_cols(int = rep(1,N)) %>% 
                   as.matrix()) %*% A_k
      
      gam <- fct_log_lik(mu_mat, sig_k, y, N, m)
      return(dplyr::tibble(gam))
      
    }) %>% 
    fct_update_clust(N)
  
  1:k %>% 
    as.list() %>% 
    purrr::map_dbl(.f = function(.x){
      
      rows_k <- (est_assign==.x) %>% 
        which()
      
      if(length(rows_k) < 1){
        return(Inf)
      }
      
      A_k <- A %>% 
        purrr::pluck(.x)
      
      x_k <- x %>% 
        dplyr::as_tibble(.name_repair = "universal") %>% 
        dplyr::slice(rows_k) %>% 
        as.matrix()
      
      y_k <- y %>% 
        dplyr::as_tibble(.name_repair = "universal") %>% 
        dplyr::slice(rows_k) %>% 
        as.matrix()
      
      mean((y_k-(cbind(x_k,1) %*% A_k))^2)
      
    }) %>% 
    sum()
}