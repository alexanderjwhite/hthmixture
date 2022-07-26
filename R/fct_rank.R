#' rank
#'
#' @inheritParams fct_gamma
#' @param sigma doc
#' @param eta doc
#' 
#' @import MASS
#'
#' @return doc
#' @export
#'
fct_rank <- function(x, y, sigma, eta){
  m <- t(x)%*%x
  m_inv <- MASS::ginv(m)
  p <- x%*%m_inv%*%t(x)
  sv <- svd(p%*%y)$d 
  max(which(sv >= sigma*eta))
}