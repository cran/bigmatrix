#----------------------------------------------------------------------------------#
# Package: bigmatrix                                                               #
# tiger.adm(): Coordinate descent method for sparse precision matrix               #
#             estimation                                                           #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 17th, 2012                                                             #
# Version: 0.9.4                                                                   #
#----------------------------------------------------------------------------------#


tiger.clime.cdadm <- function(Sigma, d, maxdf, lambda, rho, shrink, prec, max.ite){
  
  gamma = 1/rho
  d_sq = d^2
  lambda = lambda-shrink*prec
  nlambda = length(lambda)
  icov = array(0,dim=c(d,d,nlambda))
  ite_ext = rep(0,d*nlambda)
  ite_int1 = rep(0,d*nlambda)
  ite_int2 = rep(0,d*nlambda)
  x = array(0,dim=c(d,maxdf,nlambda))
  col_cnz = rep(0,d+1)
  row_idx = rep(0,d*maxdf*nlambda)
  icov_list = vector("list", nlambda)
  icov_list1 = vector("list", nlambda)
  str=.C("tiger_clime_cdadm", as.double(Sigma), as.double(icov), as.double(x), as.integer(d), 
         as.integer(ite_ext), as.integer(ite_int1), as.integer(ite_int2),as.double(lambda), 
         as.integer(nlambda), as.double(gamma), as.integer(max.ite), as.integer(col_cnz), 
         as.integer(row_idx), as.double(prec), PACKAGE="bigmatrix")
  for(i in 1:nlambda){
    icov_i = matrix(unlist(str[2])[((i-1)*d_sq+1):(i*d_sq)], byrow = FALSE, ncol = d)
    icov_list1[[i]] = icov_i
    icov_list[[i]] = icov_i*(abs(icov_i)<=abs(t(icov_i)))+t(icov_i)*(abs(t(icov_i))<abs(icov_i))
  }
  ite_ext = matrix(unlist(str[5]), byrow = FALSE, ncol = nlambda)
  ite_int1 = matrix(unlist(str[6]), byrow = FALSE, ncol = nlambda)
  ite_int2 = matrix(unlist(str[7]), byrow = FALSE, ncol = nlambda)
  x = unlist(str[3])
  col_cnz = unlist(str[12])
  row_idx = unlist(str[13])
  ite = vector("list", 3)
  ite[[1]] = ite_ext
  ite[[2]] = ite_int1
  ite[[3]] = ite_int2
  return(list(icov=icov_list, icov1=icov_list1, ite=ite, x=x, col_cnz=col_cnz, row_idx=row_idx))
}
