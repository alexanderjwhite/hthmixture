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
#' @import stats
#'
#' @examples
fct_gamma <- function(x, y, k, N, p, m, lam, rank, clust_assign, val_frac, penal_search){
  
  names <- paste0("c_",1:k)
  gamma_calc <- 1:k %>% 
    list() %>% 
    purrr::pmap(.f = function(.x){
      
      cluster_rows <- (clust_assign==.x) %>% 
        which()
      
      n_k <- cluster_rows %>% 
        length()
      
      X_k <- x %>% 
        dplyr::as_tibble(.name_repair = "universal") %>% 
        dplyr::slice(cluster_rows) %>% 
        as.matrix()
      
      Y_k <- y %>% 
        dplyr::as_tibble(.name_repair = "universal") %>% 
        dplyr::slice(cluster_rows) %>% 
        as.matrix()
      
      
      
      if((is.null(lam)|is.null(rank)) & nrow(X_k) > 3){
        
        # Cross Validate to select penalization
        
        val_rows <- cluster_rows %>% 
          dplyr::tibble() %>% 
          dplyr::slice_sample(prop = ifelse((val_frac*nrow(X_k)) < 1, 0.5, 1)) %>% 
          dplyr::pull()
        
        train_rows <- cluster_rows %>% 
          dplyr::tibble() %>% 
          dplyr::filter(!((.) %in% val_rows)) %>% 
          dplyr::pull()
        
        split_data <- fct_data_split(X_k, Y_k, ifelse((val_frac*nrow(X_k)) < 1, 0.5, val_frac))
        
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
            lam_0 <- fct_lam_coef(x_train, sigmahat_test, m, p)
            # lam_0 <- 2*sigmahat_test*max(sqrt(colSums(x_train^2)))/n_train/.r*(sqrt(.r)+2*sqrt(log(p)))
            # lam_0 <- 1
            lam <- .l*lam_0
            model <- fct_sarrs(y_train,x_train,.r, lam, "grLasso")
            error <- mean((y_test-(cbind(x_test,1) %*% model$Ahat))^2)
            dplyr::tibble(lam = lam, rank = .r, error = error, model = list(model))
          }) %>% 
          dplyr::arrange(error) %>% 
          dplyr::slice(1) %>% 
          dplyr::pull(model) %>% 
          purrr::pluck(1)
        
      } else if (nrow(X_k) > 1){
        rank_var <- fct_rank_var(x, y, n_k, p, m)
        sigmahat <- rank_var %>% 
          purrr::pluck("sigmahat")
        lam_0 <- fct_lam_coef(x, sigmahat, m, p)
        model_k <- fct_sarrs(Y_k, X_k, rank, lam_0, "grLasso")
        
        
      } else if (nrow(X_k) <= 1){
        
        # model_k <- fct_sarrs(Y_k, X_k, 1, 1, "grLasso")
        return(dplyr::tibble(-Inf) %>% stats::setNames(names[.x]))
        
      } else {
        
        model_k <- fct_sarrs(Y_k, X_k, rank, lam, "grLasso")
      }
      
      
      A_k <- model_k %>% 
        purrr::pluck("Ahat")
      
      sig_vec <- model_k %>% 
        purrr::pluck("sigvec")
      
      mu_mat <- (x %>% 
                   dplyr::bind_cols(int = rep(1,N)) %>% 
                   as.matrix()) %*% A_k
      
      
      
      gam <- fct_log_lik(mu_mat, sig_vec, y, N, m)
      return(list(gamma = dplyr::tibble(gam) %>% stats::setNames(names[.x]), A_k = A_k))
      
    })
  gamma <- gamma_calc %>% 
    purrr::map_dfc(.f = function(.x){.x %>% purrr::pluck("gamma")})
  A <- gamma_calc %>% 
    purrr::map(.f = function(.x){.x %>% purrr::pluck("A_k")})
  return(list(gamma = gamma, A = A))
}