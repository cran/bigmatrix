#----------------------------------------------------------------------------------#
# Package: bigmatrx                                                                #
# bigmatrix.slasso(): Scaled Lasso method for sparse precision matrix              #
#             estimation                                                           #
# Authors: Xingguo Li                                                              #
# Emails: <xingguo.leo@gmail.com>                                                  #
# Date: July 27th 2012                                                             #
# Version: 1.0.7                                                                   #
#----------------------------------------------------------------------------------#


tiger.slasso <- function(x, lambda){
  d=ncol(x)
  x=scale(x)
  Omega=matrix(0,d,d)
  for(j in 1:d){
    out   = SLasso(x[,j], x[,-j], lambda)
    beta  = out$beta
    sigma = out$sigma
    Omega[j,j]=1/(sigma^2)
    Omega[-j,j]=-as.vector(Omega[j,j]*beta)
    #if(j%%100==0)
    #cat("*")
    #Omega[abs(Omega)<1e-4]=0
  }
  return(Omega)
}

SLasso <- function(y, x, lambda){
  d =ncol(x)
  n = nrow(x)
  sigma.gap = 1
  sigmahat0 = 1
  iter = 1
  while(sigma.gap>1e-4){
    fit=glmnet(x,y,lambda=2*sigmahat0*lambda*sqrt(2*log(d)/n))
    sigmahat=sqrt(fit$nulldev*(1-fit$dev.ratio[1])/n)
    sigma.gap = abs(sigmahat-sigmahat0)/abs(sigmahat)
    sigmahat0 = sigmahat
    iter=iter+1
    if(iter>100){
      cat("maxIter achieved!")
      break
    }			
  }
  return(list(beta=fit$beta,sigma=sigmahat0))
}
