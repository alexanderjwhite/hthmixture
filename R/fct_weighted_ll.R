#' Compute weighted log likelihood
#'
#' @param gamma posterior
#'
#' @return weighted log likelihood
#' @export
#'
fct_weighted_ll <- function(gamma){
  
  llik <- apply(gamma, 1, max)
  llik[is.infinite(llik)] <- min(llik[is.finite(llik)])
  
  return(sum(llik))
}
