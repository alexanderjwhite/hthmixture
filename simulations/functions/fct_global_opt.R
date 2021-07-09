#' Global Likelihood Update
#'
#' @param ll_i current likelihood
#' @param global_opt_ll global likelihood
#'
#' @return updated global likelihood
#' @export
#'
#' @examples
fct_global_opt <- function(ll_i, global_opt_ll){

  if(ll_i > global_opt_ll){
    return(ll_i)
  } else {
    return(global_opt_ll)
  }
}
