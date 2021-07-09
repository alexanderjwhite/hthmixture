#' Compute weighted log likelihood
#'
#' @param gamma posterior
#'
#' @return weighted log likelihood
#' @export
#' 
#' @import dplyr
#'
#' @examples
fct_weighted_ll <- function(gamma){
  
  llik <- gamma %>% 
    rowwise() %>% 
    mutate(ll = max(c_across(where(is.numeric)))) 
  
  llik <- llik %>% 
    mutate(ll = ifelse(is.infinite(ll), min((llik %>% pull(ll))[is.finite((llik %>% pull(ll)))]), ll)) %>% 
    pull(ll) %>% 
    sum()
  
  return(llik)
}
