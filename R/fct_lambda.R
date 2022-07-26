#' lam
#'
#' @param sigma doc
#' @param p doc
#' @param n doc
#'
#' @return doc
#' @export
#'
fct_lambda <- function(sigma, p, n){
  sigma*sqrt(2*log(p)/n)
}