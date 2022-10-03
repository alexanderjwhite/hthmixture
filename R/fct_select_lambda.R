fct_select_lambda <- function(x, y, k, clust_assign = NULL, initial = FALSE, type = "single"){
  max_rank <- 3
  safe_rank <- purrr::safely(fct_rank)
  if(initial){
    M <- 100
    print("Initializing Lambda...")
    clust_assign <- fct_initialize(k,nrow(x))
  } else{
    M <- 50
    print("Selecting Lambda...")
  }
  store <- array(0, dim = c(M,k,2))
  
  for(i in 1:M){
    if(i%%5 == 0){
      print(paste0(round(i*100/M,2),"% Complete"))
    }
    if(initial){
      clust_assign <- fct_initialize(k,nrow(x))
    }
    for(j in 1:k){
      x_k <- sim$x[which(clust_assign==j),]
      y_k <- sim$y[which(clust_assign==j),]
      sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))
      
      rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta = 3)
      rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
      rank_hat <- min(rank_hat, max_rank)
      
      store[i,j,] <- fct_sarrs(y_k, x_k, r = rank_hat, lam = NULL, 
                               alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
                               ptype = "grLasso", y_sparse = TRUE)$lambda_store
    }
  }
  store_mat <- rbind(store[,,1],store[,,2])
  if(type == "single"){
    lambda <- median(store_mat)
  } else {
    lambda <- apply(store_mat, 2, median)
  }
  
  print(paste("Selected Lambda ",lambda))
  return(lambda)
}