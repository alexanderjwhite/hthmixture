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
  val_frac <- 0.2
  grid <- seq(-8,8,1)
 
  gamma <- NULL
  A <- NULL
  sig_vec <- NULL
  for (i in 1:k){
    
    cluster_rows <- which((clust_assign==i))
    n_k <- length(cluster_rows)
    x_k <- x[cluster_rows,] 
    y_k <- y[cluster_rows, ] 
    eta_k <- sqrt(2*m) + sqrt(2*min(n_k,p))
    
    
    if (n_k > 3){
      val_size <- ifelse((val_frac*n_k) < 1, (n_k*0.5), (val_frac*n_k))
      val_rows <- sample(1:length(cluster_rows), size = val_size)
      train_rows <- which(!((1:length(cluster_rows)) %in% val_rows))

      x_train <- x_k[train_rows,]
      x_test <- x_k[val_rows,]
      y_train <- y_k[train_rows,]
      y_test <- y_k[val_rows]

      sigma_hat <- fct_sigma(y_k, n_k, m)
      rank_hat <- fct_rank(x_k, y_k, sigma_hat, eta_k)
      lam_univ <- fct_lambda(sigma_hat, p, n_k)
      lam_grid <- (2^(grid/2))*lam_univ
      
      models <- NULL
      errors <- rep(0,length(lam_grid))
      for (j in 1:length(lam_grid)){
        print(j/length(lam_grid))
        model_j <- fct_sarrs(y_train, x_train, rank_hat, lam_grid[j], alpha, beta, sigma_hat, "grLasso")
        errors[j] <- mean((y_test-(cbind(x_test,1) %*% model_j$Ahat))^2)
        models <- c(models,list(model_j))
      }

    }
    else if (n_k > 1){
      # if (n_k > 1){
      
      sigma_hat <- fct_sigma(y_k, n_k, m)
      rank_hat <- fct_rank(x_k, y_k, sigma_hat, eta_k)
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