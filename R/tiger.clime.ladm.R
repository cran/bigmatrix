#----------------------------------------------------------------------------------#
# Package: bigmatrix                                                               #
# tiger.adm(): Coordinate descent method for sparse precision matrix               #
#             estimation                                                           #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 12th, 2012                                                             #
# Version: 0.9.4                                                                   #
#----------------------------------------------------------------------------------#


tiger.clime.ladm <- function(Sigma, d, maxdf, lambda, rho, shrink, prec, max.ite, verbose){
  
  SS = Sigma%*%Sigma
  gamma = eigen(SS)$values[1]
  d_sq = d^2
  lambda = lambda-shrink*prec
  nlambda = length(lambda)
  icov = array(0,dim=c(d,d,nlambda))
  ite = rep(0,d*nlambda)
  x = array(0,dim=c(d,maxdf,nlambda))
  col_cnz = rep(0,d+1)
  row_idx = rep(0,d*maxdf*nlambda)
  icov_list = vector("list", nlambda)
  icov_list1 = vector("list", nlambda)
  if (verbose) verbose=1
  else verbose=0
  str=.C("tiger_clime_ladm", as.double(Sigma), as.double(SS), as.double(icov), as.double(x), 
         as.integer(d), as.integer(ite), as.double(lambda), as.integer(nlambda), 
         as.double(gamma), as.integer(max.ite), as.double(rho), as.integer(col_cnz), 
         as.integer(row_idx), as.double(prec), PACKAGE="bigmatrix")
  for(i in 1:nlambda){
    icov_i = matrix(unlist(str[3])[((i-1)*d_sq+1):(i*d_sq)], byrow = FALSE, ncol = d)
    icov_list1[[i]] = icov_i
    icov_list[[i]] = icov_i*(abs(icov_i)<=abs(t(icov_i)))+t(icov_i)*(abs(t(icov_i))<abs(icov_i))
  }
  ite = matrix(unlist(str[6]), byrow = FALSE, ncol = nlambda)
  x = unlist(str[4])
  col_cnz = unlist(str[12])
  row_idx = unlist(str[13])
  return(list(icov=icov_list, icov1=icov_list1,ite=ite, x=x, col_cnz=col_cnz, row_idx=row_idx))
}
