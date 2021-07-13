#' Compute Posterior
#'
#' @param x matrix
#' @param y matrix
#' @param k clusters
#' @param N observations
#' @param p dimension
#' @param m dimension
#' @param lam penalty param
#' @param rank rank
#' @param clust_assign current assignment
#' @param val_frac fraction validation
#' @param penal_search search vector
#'
#' @return tibble of posterior
#' @export
#' 
#' @import dplyr
#' @import purrr
#'
#' @examples
fct_gamma <- function(x, y, k, N, p, m, lam, rank, clust_assign, val_frac, penal_search){
  
  names <- paste0("c_",1:k)
  1:k %>% 
    list() %>% 
    purrr::pmap_dfc(.f = function(.x){
      
      cluster_rows <- (clust_assign==.x) %>% 
        which()
      
      n_k <- cluster_rows %>% 
        length()
      
      X_k <- x %>% 
        as_tibble() %>% 
        slice(cluster_rows) %>% 
        as.matrix()
      
      Y_k <- y %>% 
        as_tibble() %>% 
        slice(cluster_rows) %>% 
        as.matrix()
      
      
      
      if(is.null(lam)|is.null(rank)){
        
        # Cross Validate to select penalization
        
        val_rows <- cluster_rows %>% 
          tibble() %>% 
          slice_sample(prop = val_frac) %>% 
          pull()
        
        train_rows <- cluster_rows %>% 
          tibble() %>% 
          filter(!((.) %in% val_rows)) %>% 
          pull()

        split_data <- fct_data_split(X_k, Y_k, val_frac)
        
        x_train <- split_data %>% 
          purrr::pluck("x_train")
        
        x_test <- split_data %>% 
          purrr::pluck("x_test")
        
        y_train <- split_data %>% 
          purrr::pluck("y_train")
        
        y_test <- split_data %>% 
          purrr::pluck("y_test")
        
        n_train <- split_data %>% 
          purrr::pluck("n_train")
        
        n_test <- split_data %>% 
          purrr::pluck("n_test")
        
        rank_var_test <- fct_rank_var(x_train, y_train, n_train, p, m)
        rank_search <- 1

        sigmahat_test <- rank_var_test %>% 
          purrr::pluck("sigmahat")
        
        
        grid_search <- expand.grid(lam = penal_search, r = rank_search)
        model_k <- list(grid_search$lam, grid_search$r) %>% 
          purrr::pmap_dfr(.f = function(.l, .r){
            # lam_0 <- 2*sigmahat_test*max(sqrt(colSums(x_train^2)))/n_train/.r*(sqrt(.r)+2*sqrt(log(p)))
            lam_0 <- fct_lam_coef(x_train, sigmahat_test, m, p)
            lam <- .l*lam_0
            model <- fct_sarrs(y_train,x_train,.r, lam, "grLasso")
            error <- mean((y_test-(cbind(x_test,1) %*% model$Ahat))^2)
            
            tibble(lam = lam, rank = .r, error = error, model = list(model))
          }) %>% 
          arrange(error) %>% 
          slice(1) %>% 
          pull(model) %>% 
          purrr::pluck(1)

      } else {
        
        model_k <- fct_sarrs(Y_k, X_k, rank, lam, "grLasso")
      }
      
      
      A_k <- model_k %>% 
        purrr::pluck("Ahat")
      
      sig_vec <- model_k %>% 
        purrr::pluck("sigvec")
      
      mu_mat <- (x %>% 
                   bind_cols(int = rep(1,N)) %>% 
                   as.matrix()) %*% A_k
      
      
      
      gam <- fct_log_lik(mu_mat, sig_vec, y, N, m)
      return(tibble(gam) %>% setNames(names[.x]))
      
    })
  
}
