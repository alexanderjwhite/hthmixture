#' Compute penalization coefficient
#'
#' @param x matrix
#' @param sigma value
#' @param m integer
#' @param p integer
#'
#' @return penalization coefficient
#' 
#' @export 
#'
#' @examples fct_lam_coef(matrix(1:9, nrow = 3), sigma = 1, m = 3, p = 3)
fct_lam_coef <- function(x, sigma, m, p){
 xmax <- apply(x, 2, function(.x){norm(.x, type = "2")}) %>% 
   max()
 return(4*sigma*xmax*(sqrt(m) + sqrt(4*log(p))))
}