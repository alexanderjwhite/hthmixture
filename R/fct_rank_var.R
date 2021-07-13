#' Estimate noise variance
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
  
  sig <- svd(Y_k)$d
  sig <- sig[sig > 0]
  sigmahat <- median(sig)/sqrt(max(m,n_k))
  
  return(list(sigmahat = sigmahat, rhat = 1))
}
