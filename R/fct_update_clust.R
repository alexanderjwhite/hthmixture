#' Update Cluster
#'
#' @param N observations
#' @param gamma posterior
#'
#' @return updated cluster assignments
#' @export
#' 
#' @import dplyr
#' @import purrr
#'
#' @examples
fct_update_clust <- function(gamma, N){
  1:N %>% 
    purrr::map_int(
      .f = function(.x){
        res <- (gamma %>% as.matrix())[.x,] %>% which.max()
      }
    )
}
