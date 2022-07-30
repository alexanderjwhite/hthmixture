#' y
#'
#' @inheritParams fct_gamma
#' @param n doc
#' @param m doc
#'
#' @return doc
#' @export
#'
fct_sigma <- function(y, n, m){
  sv <- svd(y)$d
  sigma_hat <- median(sv[sv!=0])/sqrt(max(n,m))
}