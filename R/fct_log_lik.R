#' Compute Log Likelihood
#'
#' @param mu_mat matrix
#' @param sig_vec vector
#' @param y value
#' @param N observations
#' @param m dimension
#'
#' @return array
#' @export
#' 
#' @import purrr dplyr stats
#'
#' @examples
fct_log_lik <- function(mu_mat, sig_vec, y, N, m){
  
  gam <- array(dim = N)
  for(i in 1:N){
    mu_i <- mu_mat[i,]
    gam[i] <- 1:m %>% 
      purrr::map_dbl(
        .f = function(.y){
          stats::dnorm(y[i,.y], mean = mu_i[.y], sd = sig_vec[.y])
        }
      ) %>%
      log() %>% 
      sum() 
  }
  return(gam)
}
