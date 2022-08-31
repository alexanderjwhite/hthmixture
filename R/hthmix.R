#' Main mixrrr function
#'
#' @param x matrix
#' @param y matrix
#' @param k number of clusters
#' @param nstart number of times to initialize
#' @param init_assign initial assignment
#' @param selection selection method
#' @param alpha sparsity tuning parameter
#' @param beta sparsity tuning parameter
#' @param y_sparse treat y as sparse?
#' @param maxiter maximum number of iterations
#' @param verbose print progress?
#'
#' @return results
#' @export
#' 
#' @import dplyr purrr
#' @importFrom rlang .data
#'
hthmix <- function(x, y, k, nstart = 1, init_assign = NULL, selection = c("universal", "cv"), alpha = 2*sqrt(3), beta = 1, y_sparse = FALSE, maxiter = 1e3, verbose = TRUE){
  
  N <- x %>% nrow()
  
  assignments <- NULL
  likelihood <- NULL
  A <- NULL
  for (i in 1:nstart){
    print(paste("start: ",i))
    if(is.null(init_assign)){
      clust_assign <- fct_initialize(k, N)
    } else {
        clust_assign <- init_assign
    }
    # print(clust_assign)
    model <- fct_hthmix_comp(x, y, k, maxiter, clust_assign, selection, alpha, beta, y_sparse)
    likelihood <- c(likelihood,model$ll)
    assignments <- c(assignments,list(model$assign))
    A <- c(A,list(model$A))
    assign_store <- model$assign_store
    
  }
  best <- which.max(likelihood)
 
  result <- list(llik = likelihood[best], assign = assignments[[best]], A = A[[best]], assign_store = assign_store)
  return(result)
}
