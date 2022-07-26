#' Initial Assignment
#'
#' @param k number of clusters
#' @param N number of observations
#'
#' @return assignment vector
#' @export
#' 
#' @import dplyr
#' @import purrr
#' @import stats
#'
fct_initialize <- function(k, N){
  
  init_int <- cumsum(rep(1/k,k))
  init_rand_assign <- stats::runif(N)
  
 init_rand_assign %>% 
    purrr::map_int(.f = function(.x){
      clust <- (.x <= init_int) %>% 
        which() %>% 
        min()
      return(clust)
    })
  
}
