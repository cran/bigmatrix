#----------------------------------------------------------------------------------#
# Package: bigmatrx                                                                #
# bigmatrix.slasso(): Scaled Lasso method for sparse precision matrix              #
#             estimation                                                           #
# Authors: Xingguo Li                                                              #
# Emails: <xingguo.leo@gmail.com>                                                  #
# Date: Oct 15th, 2012                                                             #
# Version: 0.9.6                                                                   #
#----------------------------------------------------------------------------------#


tiger.slasso.cdadm <- function(data, n, d, maxdf, lambda, shrink, prec, max.ite){
  
  nlambda = length(lambda)
  lambda = lambda-shrink*prec
  X1=data
  d_sq=d^2
  X1 = X1 - matrix(rep(colMeans(X1),n), nrow=n, byrow=TRUE)
  Gamma=diag(diag(t(X1)%*%X1/n))
  Q = diag(1/sqrt(diag(Gamma)))
  Omega=array(0,dim=c(d,d,nlambda))
  Z=X1%*%(diag(1/sqrt(diag(Gamma))))
  
  icov = array(0,dim=c(d,d,nlambda))
  ite_ext = rep(0,d*nlambda)
  ite_int = rep(0,d*nlambda)
  x = array(0,dim=c(d,maxdf,nlambda))
  col_cnz = rep(0,d+1)
  row_idx = rep(0,d*maxdf*nlambda)
  icov_list = vector("list", nlambda)
  icov_list1 = vector("list", nlambda)
  str=.C("tiger_slasso_cdadm", as.double(Z), as.double(icov), as.double(x), as.integer(d), as.integer(n), 
         as.integer(ite_ext), as.integer(ite_int), as.double(lambda), as.integer(nlambda), as.integer(max.ite), 
         as.integer(col_cnz),as.integer(row_idx), as.double(prec), PACKAGE="bigmatrix")
  for(i in 1:nlambda){
    icov_i = Q%*%(matrix(unlist(str[2])[((i-1)*d_sq+1):(i*d_sq)], byrow = FALSE, ncol = d))%*%Q
    icov_list1[[i]] = icov_i
    icov_list[[i]] = icov_i*(abs(icov_i)<=abs(t(icov_i)))+t(icov_i)*(abs(t(icov_i))<abs(icov_i))
  }
  ite_ext = matrix(unlist(str[6]), byrow = FALSE, ncol = nlambda)
  ite_int = matrix(unlist(str[7]), byrow = FALSE, ncol = nlambda)
  x = unlist(str[3])
  col_cnz = unlist(str[11])
  row_idx = unlist(str[12])
  ite = vector("list", 2)
  ite[[1]] = ite_ext
  ite[[2]] = ite_int
  return(list(icov=icov_list, icov1=icov_list1, ite=ite, x=x, col_cnz=col_cnz, row_idx=row_idx))
}