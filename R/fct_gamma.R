#' Compute Posterior
#'
#' @inheritParams fct_hthmix_comp
#' @param N doc
#' @param clust_assign doc
#' 
#' @import dplyr stats purrr
#'
#' @return doc
#' @export
#'
fct_gamma <- function(x, y, k, N, clust_assign){
  
  alpha <- 2*sqrt(3)
  beta <- 1
  p <- dim(x)[2]
  m <- dim(y)[2]
 
  gamma <- NULL
  A <- NULL
  sig_vec <- NULL
  for (i in 1:k){
    
    cluster_rows <- which((clust_assign==i))
    n_k <- length(cluster_rows)
    x_k <- x[cluster_rows,] 
    y_k <- y[cluster_rows, ] 
    eta_k <- sqrt(2*m) + sqrt(2*min(n_k,p))

    if (n_k > 1){
      
      sigma_hat <- fct_sigma(y_k, n_k, m)
      rank_hat <- fct_rank(x, y, sigma_hat, eta_k)
      lam_univ <- fct_lambda(sigma_hat, p, n_k)
      model_k <- fct_sarrs(y_k, x_k, rank_hat, lam_univ, alpha, beta, sigma_hat, "grLasso")
      
      A_k <- model_k$Ahat 
      sigvec <- model_k$sigvec
      mu_mat <- cbind(x,1) %*% A_k
      gam <- fct_log_lik(mu_mat, sigvec, y, N, m)
      
      gamma <- cbind(gamma, gam)
      A <- c(A,list(A_k))
      sig_vec <- c(sig_vec,list(sigvec))
      
    } else {
      # return(rep(-Inf,k))
      # return(list(gamma = matrix(rep(-Inf,k),nrow=1), A = NULL, sig_vec = NULL))
      gamma <- cbind(gamma, rep(-Inf,N))
      A <- c(A,list(NULL))
      sig_vec <- c(sig_vec,list(NULL))
      
    }
    
   
    
    
  }
  return(list(gamma = gamma, A = A, sig_vec = sig_vec))
}