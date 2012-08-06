#----------------------------------------------------------------------------------#
# Package: bigmatrx                                                                #
# bigmatrix.slasso(): Scaled Lasso method for sparse precision matrix              #
#             estimation                                                           #
# Authors: Xingguo Li                                                              #
# Emails: <xingguo.leo@gmail.com>                                                  #
# Date: Aug 6th, 2012                                                              #
# Version: 0.9.2                                                                   #
#----------------------------------------------------------------------------------#


tiger.slasso <- function(data, lambda, standardize=FALSE){
  nlambda = length(lambda)
  X=data
  d=ncol(X)
  n=nrow(X)
  X = X - matrix(rep(colMeans(X),n), nrow=n, byrow=TRUE)
  Gamma=diag(diag(t(X)%*%X/n))
  Omega=array(0,dim=c(d,d,nlambda))
  Z=X%*%(diag(1/sqrt(diag(Gamma))))
  for(j in 1:d){
    tau0 = rep(1,nlambda)
    gap = 1
    iter = 1
    while(gap > 1e-4){
      out  = glmnet(Z[,-j],Z[,j],lambda=tau0*lambda,standardize=FALSE)
      beta = as.matrix(out$beta)
      tau1 = sqrt(colSums((matrix(rep(Z[,j],nlambda),ncol=nlambda)-Z[,-j]%*%beta)^2)/n)
      gap  = max(abs(tau1-tau0)/tau0)
      tau0 = tau1
      if(iter>300){
        cat(sprintf("MaxIter achieved!--\n"))
        break
      }  		
      iter=iter+1
    }
    beta  = as.matrix(out$beta)
    Omega[ j,j,]=1/(tau0^2)
    Omega[-j,j,]=-beta/matrix(rep(tau0^2,d-1),byrow=TRUE,nrow=d-1)
  }  
  Q = diag(1/sqrt(diag(Gamma)))   
  for(j in 1:nlambda)
    Omega[,,j] = Q%*%Omega[,,j]%*%Q
  return(Omega)
}