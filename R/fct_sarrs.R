#' Compute SARRS algorithm
#'
#' @param Y response matrix
#' @param X design matrix
#' @param r rank
#' @param lam penalization param
#' @param ptype penalty type
#' @param V0 v0
#'
#' @return estimate of coefficients
#' @export
#' 
#' @import grpreg
#' @import stats
#'
#' @examples
fct_sarrs <- function(Y,X,r,lam,ptype,V0=NULL){
  n= dim(X)[1]
  p= dim(X)[2]
  m= dim(Y)[2]
  group = rep(1:(p+1),r)
  
  if(is.null(V0)){
    V0 = svd(Y,nu=r,nv=r)$v
  }
  
  XX = kronecker(diag(rep(1,r)),cbind(X,1))
  YY = Y %*% V0
  YY = as.vector(YY)
  fit1 = grpreg::grpreg(XX,YY,group,lambda=lam, penalty= ptype, family="gaussian")
  B1= matrix(fit1$beta[-1],nrow=p+1,ncol=r)
  B1[p+1,]=B1[p+1,]+fit1$beta[1]
  XB=cbind(X,1) %*% matrix(B1,nrow=p+1,ncol=r)
  
  U1= svd(XB,nu=r,nv=r)$u
  tmp= U1%*%t(U1)%*%Y
  V1= svd(tmp,nu=r,nv=r)$v
  
  YY = Y %*% V1
  YY = as.vector(YY)
  fit2 = grpreg::grpreg(XX,YY,group,lambda=lam, penalty= ptype, family="gaussian")
  B2= matrix(fit2$beta[-1],nrow=p+1,ncol=r)
  B2[p+1,]=B2[p+1,]+fit2$beta[1]
  
  Ahat= B2 %*% t(V1)
  sigvec = apply(Y- cbind(X,rep(1,nrow(X)))%*% Ahat, 2, stats::sd)
  return(list(Ahat=Ahat,sigvec=sigvec))
}