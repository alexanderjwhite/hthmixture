#' Title
#'
#' @param X_k matrix
#' @param Y_k matrix
#' @param n_k observations
#' @param p dimension
#' @param m dimension
#'
#' @return estimated noise and rank
#' @export
#'
#' @examples
fct_rank_var <- function(X_k, Y_k, n_k, p, m){
  
  # proj_mat <- X_k %*% MASS::ginv(t(X_k) %*% X_k) %*% t(X_k)
  # sigmahat <- sqrt(sum((Y_k-proj_mat%*%Y_k)^2)/(n_k*m-min(n_k,p)*m))
  # rhat <- sum(svd(proj_mat %*% Y_k)$d > sigmahat*(sqrt(2*m) + sqrt(2*min(n_k, p))))
  
  return(list(sigmahat = 1, rhat = 1))
}
