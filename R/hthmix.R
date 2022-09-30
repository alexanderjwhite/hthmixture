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
#' @param rank force to be known rank
#' @param max_rank limit rank
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
# hthmix <- function(x, y, k, nstart = 1, init_assign = NULL, selection = "universal", "cv"), 
#                    alpha = 2*sqrt(3), beta = 1, y_sparse = FALSE, rank = NULL, max_rank = 3, maxiter = 1e3, verbose = TRUE,
#                    temp = 1000, p1 = 0.9, p2 = 0.95, sim_N = 500, checks = 5, true = NULL){
hthmix <- function(x, y, k, nstart, init_assign = NULL, init_lambda, alt_iter, anneal_iter, 
                   em_iter, temp, mu, eps, accept_prob, sim_N){
                     
  
  N <- x %>% nrow()
  
  assignments <- NULL
  likelihood <- NULL
  ll_store <- NULL
  A <- NULL
  for (i in 1:nstart){
    print(paste("start: ",i))
    # if(is.null(init_assign)){
    #   clust_assign <- fct_initialize(k, N)
    # } else {
    #     clust_assign <- init_assign
    # }
    # print(clust_assign)
    # model <- fct_hthmix_comp(x, y, k, maxiter, clust_assign, selection, alpha, beta, y_sparse, rank, max_rank, temp, p1, p2, sim_N, checks, true)
    model <- fct_alt_optimize(x, y, k, init_assign, init_lambda, alt_iter, anneal_iter, em_iter, temp, mu, eps, accept_prob, sim_N)
    likelihood <- c(likelihood,model$ll)
    assignments <- c(assignments,list(model$assign))
    ll_store <- c(ll_store,list(model$ll_store))
    A <- c(A,list(model$A))
    assign_store <- model$assign_store
    
  }
  best <- which.max(likelihood)
 
  result <- list(llik = likelihood[best], assign = assignments[[best]], A = A[[best]], assign_store = assign_store, ll_store = ll_store[[best]])
  return(result)
}
