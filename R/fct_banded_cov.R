#' Create Banded Covariance Matrix
#'
#' @param dim integer dimension of matrix
#' @param type large or small
#' @param perc_large if large, percentage of large entries
#'
#' @return matrix
#' @export 
#'
#' @examples
fct_banded_cov <- function(dim, type = c("large", "small"), perc_large = 0.01){
  sigma <- diag(nrow = dim, ncol = dim)
  if(type=="small"){
    coef <- function(n){runif(n, 0.01, 0.1)} 
  } else {
    coef <- function(n){
      if(rbinom(1, size=1, prob=perc_large)==1){
        runif(n, 0.9, 0.99)
      } else {
        runif(n, 0.01, 0.1)
      }
    }
  }
  for(i in 1:dim){
    for(j in c((i-1),(i+1))){
      if(i <= dim & j <= dim){
        sigma[i,j] <- coef(1)
      }
      
    }
  }
  return(sigma)
}