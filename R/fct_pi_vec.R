#' Compute pi vector
#'
#' @param clust_assign current cluster assignment
#' @param k clusters
#' @param N observations
#'
#' @return pi vector
#' @export
#' 
#' @import purrr
#' @import dplyr
#'
#' @examples
fct_pi_vec <- function(clust_assign, k, N){
  
  1:k %>% 
    purrr::map_dbl(
      .f = function(.x){
        sum(clust_assign == .x)/N
      }
    )
      
}
