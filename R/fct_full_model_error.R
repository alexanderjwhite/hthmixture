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
#' @examples
fct_full_model_error <- function(x, y, k, model){
  est_assign <- model %>% 
    purrr::pluck("result", "assign", "final_assign")
  
  A <- model %>% 
    purrr::pluck("result", "assign", "A")
  
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