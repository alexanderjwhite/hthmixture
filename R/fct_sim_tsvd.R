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
fct_sim_tsvd <- function(n = 100, p = 200, m = 250, b = 1, d = 20, rank = 1, h = 0.2, case = "independent"){
  # u <- rep(0,p)
  # v <- rep(0,m)
  # u_index <- sample(1:p,25, replace = FALSE)
  # v_index <- sample(1:m,25, replace = FALSE)
  
  
  u1 <- c(sample(c(1,-1), 5, replace = TRUE)*runif(5,h,1),rep(0,20))
  u2 <- c(rep(0,5),sample(c(1,-1), 5, replace = TRUE)*runif(5,h,1),rep(0,15))
  u3 <- c(rep(0,10),sample(c(1,-1), 5, replace = TRUE)*runif(5,h,1),rep(0,10))
  u <- matrix(cbind(u1,u2,u3)[,1:rank],ncol = rank)
  # u[u_index] <- u1
  
  v1 <- c(sample(c(1,-1), 10, replace = TRUE),rep(0,15))
  v2 <- c(rep(0,12),sample(c(1,-1), 10, replace = TRUE),rep(0,3))
  v3 <- c(sample(c(1,-1), 6, replace = TRUE),v1[7:8],-v1[9:10],sample(c(1,-1), 2, replace = TRUE),
          -v1[13:14],v1[15:16],rep(0,9))
  v <- matrix(cbind(v1,v2,v3)[,1:rank],ncol = rank)
  # v[v_index] <- v1
  # u <- matrix(u)
  # v <- matrix(v)
  if(p > 25){u <- rbind(u,matrix(0,nrow = p-25,ncol=rank))}
  if(m > 25){v <- rbind(v,matrix(0,nrow = m-25,ncol=rank))}
  # 
  # if(rank > 1){
  #   d <- diag(c(20,10,15)[1:rank])
  # } else {
  #   d <- 20
  # }
  
  c <- u%*%d%*%t(v)
  
  
  if(case=="independent"){rho = 0} else {rho = sample(seq(0.125,0.625,0.125),1)}
  
  Sigma <- matrix(1,p,p)
  for(j in 1:p) for (k in 1:p) Sigma[j,k] <- rho^abs(j-k)
  
  
  
  x <- MASS::mvrnorm(n,rep(0,p),Sigma)	
  e <- matrix(rnorm(n*m, sd = sum(diag(t(c)%*%Sigma%*%c))/(n*m*b)),n,m)
  y <- x%*%c + e
  
  return(list(Y=y,X=x,A=c,rank=rank,sig=Sigma,E=e))
}
