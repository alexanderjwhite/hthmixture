#' Compute weighted log likelihood
#'
#' @param gamma posterior
#'
#' @return weighted log likelihood
#' @export
#' 
#' @import dplyr
#' @import tidyselect
#'
#' @examples
fct_weighted_ll <- function(gamma){
  
  llik <- gamma %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(ll = max(c_across(tidyselect:::where(is.numeric)))) 
  
  llik <- llik %>% 
    dplyr::mutate(ll = ifelse(is.infinite(ll), min((llik %>% pull(ll))[is.finite((llik %>% pull(ll)))]), ll)) %>% 
    dplyr::pull(ll) %>% 
    sum()
  
  return(llik)
}
