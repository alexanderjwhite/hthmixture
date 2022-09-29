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
#' @import grpreg stats
#'
fct_sarrs <- function(Y,X,r,lam=NULL, alpha = 2*sqrt(3), beta = 1, sigma, ptype = "grLasso", y_sparse = TRUE){
  n= dim(X)[1]
  p= dim(X)[2]
  m= dim(Y)[2]
  group = rep(1:(p+1),r)
  Y_thresh <- matrix(0, ncol = m, nrow = n)
  lambda_store <- c(0,0)
  
  
  thresh_1 <- sigma^2*(n+alpha*sqrt(n*log(max(p,m))))
  j0 <- which(apply(Y, 2, function(x){sum(x^2)}) >= thresh_1)
  Y0 <- Y_thresh
  if(y_sparse){
    Y0[,j0] <- Y[,j0]
  } else {
    Y0 <- Y
  }
  
  V0 = svd(Y0,nu=r,nv=r)$v
  
  XX = kronecker(diag(rep(1,r)),cbind(X,1))
  YY = Y %*% V0
  YY = as.vector(YY)
  if(is.null(lam)){
    fit1_cv <- grpreg::cv.grpreg(XX,YY,group, penalty= ptype, nfolds = 10, family="gaussian")
    fit1 = grpreg::grpreg(XX,YY,group,lambda=fit1_cv$lambda.min, penalty= ptype, family="gaussian")
    # print(fit1_cv$lambda.min)
    lambda_store[1] <- fit1_cv$lambda.min
  } else {
    fit1 = grpreg::grpreg(XX,YY,group,lambda=lam, penalty= ptype, family="gaussian")
    lambda_store[1] <- lam
  }
  
  B1= matrix(fit1$beta[-1],nrow=p+1,ncol=r)
  B1[p+1,]=B1[p+1,]+fit1$beta[1]
  XB=cbind(X,1) %*% matrix(B1,nrow=p+1,ncol=r)
  U1= svd(XB,nu=r,nv=r)$u
  
  thresh_2 <- beta*sigma^2*(r + 2*sqrt(3*r*log(max(p,m))) + 6*log(max(p,m)))
  j1_tmp <- which(apply(Y, 2, function(x){sum((t(U1)%*%matrix(x))^2)}) > thresh_2)
  j1 <- sort(unique(c(j0, j1_tmp)))
  Y1 <- Y_thresh
  
  if(y_sparse){
    Y1[,j1] <- Y[,j1]
  } else {
    Y1 <- Y
  }
  
  
  tmp= U1%*%t(U1)%*%Y1
  V1= svd(tmp,nu=r,nv=r)$v
  
  YY = Y %*% V1
  YY = as.vector(YY)
  
  if(is.null(lam)){
    fit2_cv <- grpreg::cv.grpreg(XX,YY,group, penalty= ptype, family="gaussian")
    fit2 = grpreg::grpreg(XX,YY,group,lambda=fit2_cv$lambda.min, penalty= ptype, nfolds = 10, family="gaussian")
    # print(fit2_cv$lambda.min)
    lambda_store[2] <- fit2_cv$lambda.min
  } else {
    fit2 = grpreg::grpreg(XX,YY,group,lambda=lam, penalty= ptype, nfolds = 10, family="gaussian")
    lambda_store[2] <- lam
  }
  B2= matrix(fit2$beta[-1],nrow=p+1,ncol=r)
  B2[p+1,]=B2[p+1,]+fit2$beta[1]
  
  Ahat= B2 %*% t(V1)
  sigvec = apply(Y- cbind(X,rep(1,nrow(X)))%*% Ahat, 2, stats::sd)
  return(list(Ahat=Ahat,sigvec=sigvec,lambda_store=lambda_store))
}