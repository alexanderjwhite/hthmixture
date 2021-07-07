#' Boltzmann Perturbation
#'
#' @param lik vector of likelihoods
#'
#' @return perturbation fraction
#' @export
#' 
#' @import dplyr
#'
#' @examples
fct_boltzmann_pert <- function(lik){
  if(is.null(lik)|length(lik)==1){
    return(1)
  }
  rel <- ((lik %>% lead())/lik) %>% .[!is.na(.)]
  kt <- 1
  pert <- 1-exp(rel[length(rel)]/kt)/sum(exp(rel/kt))
  return(pert)
  
}
