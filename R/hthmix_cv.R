#' Cross Validation for Rank/K Selection in HTHMix
#'
#' @param x matrix
#' @param y matrix
#' @param k_min minimum k to search
#' @param k_max maximum k to search
#' @param r_min minimum rank to search
#' @param r_max maximum rank to search
#' @param lam penalization parameter
#' @param chains number of chains
#' @param maxiter maximum iterations
#' @param verbose TRUE or FALSE
#' @param val_frac fraction for validation
#' @param penal_search penalty search vector
#'
#' @return list of results
#' @import purrr
#' @export
#'
#' @examples
hthmix_cv <- function(x, y, k_min = 2, k_max = 8, r_min = 1, r_max = 3, lam = NULL, chains = 50, maxiter = 100, verbose = TRUE, val_frac = 0.2, penal_search = 1:20/20){
  search_grid <- expand.grid(k = k_min:k_max, r = r_min:r_max)
  
  list(search_grid$k, search_grid$r) %>% 
    purrr::pmap(.f = function(.k, .r){
      model <- hthmix(x, y, .k, lam = NULL, rank = .r, chains = chains, maxiter = maxiter, verbose = verbose, val_frac = val_frac, penal_search = penal_search)
      error <- fct_full_model_error(x, y, .k, model)
      return(list(k = .k, r = .r, model = model, error = error))
    })
}