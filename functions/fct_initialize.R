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
#'
#' @examples
fct_initialize <- function(k, N){
  
  init_int <- rep(1/k,k) %>% cumsum()
  init_rand_assign <- runif(N)
  
 init_rand_assign %>% 
    purrr::map_int(.f = function(.x){
      clust <- (.x <= init_int) %>% 
        which() %>% 
        min()
      return(clust)
    })
  
}
