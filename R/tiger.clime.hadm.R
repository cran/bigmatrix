#----------------------------------------------------------------------------------#
# Package: bigmatrix                                                               #
# tiger.adm(): Coordinate descent method for sparse precision matrix               #
#             estimation                                                           #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Sep 4th, 2012                                                              #
# Version: 0.9.5                                                                   #
#----------------------------------------------------------------------------------#


tiger.clime.hadm <- function(Sigma, d, maxdf, lambda, rho, shrink, prec, max.ite){
  
  gamma = 1/rho
  d_sq = d^2
  lambda = lambda-shrink*prec
  nlambda = length(lambda)
  icov = array(0,dim=c(d,d,nlambda))
  ite_ext = rep(0,d*nlambda)
  ite_ext2 = rep(0,d*nlambda)
  ite_int1 = rep(0,d*nlambda)
  ite_int2 = rep(0,d*nlambda)
  x = array(0,dim=c(d,maxdf,nlambda))
  col_cnz = rep(0,d+1)
  row_idx = rep(0,d*maxdf*nlambda)
  icov_list = vector("list", nlambda)
  icov_list1 = vector("list", nlambda)
  SS = Sigma%*%Sigma
  gamma2 = eigen(SS)$values[1]
  str=.C("tiger_clime_hadm", as.double(Sigma), as.double(SS), as.double(icov), as.double(x), as.integer(d), 
         as.integer(ite_ext), as.integer(ite_ext2), as.integer(ite_int1), as.integer(ite_int2),as.double(lambda), 
         as.integer(nlambda), as.double(gamma), as.double(gamma2), as.integer(max.ite), as.double(rho), as.integer(col_cnz), 
         as.integer(row_idx), as.double(prec), PACKAGE="bigmatrix")
  for(i in 1:nlambda){
    icov_i = matrix(unlist(str[3])[((i-1)*d_sq+1):(i*d_sq)], byrow = FALSE, ncol = d)
    icov_list1[[i]] = icov_i
    icov_list[[i]] = icov_i*(abs(icov_i)<=abs(t(icov_i)))+t(icov_i)*(abs(t(icov_i))<abs(icov_i))
  }
  ite_ext = matrix(unlist(str[6]), byrow = FALSE, ncol = nlambda)
  ite_ext2 = matrix(unlist(str[7]), byrow = FALSE, ncol = nlambda)
  ite_int1 = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite_int2 = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  x = unlist(str[4])
  col_cnz = unlist(str[16])
  row_idx = unlist(str[17])
  ite = vector("list", 4)
  ite[[1]] = ite_ext
  ite[[2]] = ite_ext2
  ite[[3]] = ite_int1
  ite[[4]] = ite_int2
  return(list(icov=icov_list, icov1=icov_list1,ite=ite, x=x, col_cnz=col_cnz, row_idx=row_idx))
}
