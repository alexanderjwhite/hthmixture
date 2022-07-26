#' Simulate SARRS setting
#'
#' @param n sample size
#' @param a_rows non-zero coefficient rows
#' @param p number of x features
#' @param m number of y features
#' @param r rank
#' @param rho correlation coefficient
#' @param sigma error sd
#' @param b coefficient
#'
#' @return
#' @export
#' 
#' @import MASS
#'
#' @examples
fct_sim_mixrrr <- function(n,a_rows, a_cols, p, m, rank, rho,sigma,b){
  
  A <- matrix(0,p,m)
  nrow <- length(a_rows)
  ncol <- length(a_cols)
  B0 <- matrix(rnorm(nrow*rank),nrow,rank)
  B1 <- matrix(rnorm(rank*ncol),rank,ncol)
  A[a_rows,a_cols] <- b*B0%*%B1
  
  if(sigma == "small"){
    sigma <- 1
    Sigma <- fct_banded_cov(p, type="small")
  } else if (sigma == "large"){
    sigma <- 1
    Sigma <- fct_banded_cov(p, type="large")
  } else {
    sigma <- as.numeric(sigma)
    Sigma <- matrix(1,p,p)
    for(j in 1:p) for (k in 1:p) Sigma[j,k] <- rho^abs(j-k)
  }
  
  X <- MASS::mvrnorm(n,rep(0,p),Sigma)	
  E <- matrix(rnorm(n*m, sd = sigma),n,m)
  Y <- X%*%A + E
  
  return(list(Y=Y,X=X,A=A,rank=rank,sig=Sigma,E=E))
}