#----------------------------------------------------------------------------------#
# Package: bigmatrix                                                               #
# tiger(): The user interface for tiger()                                          #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 17th, 2012                                                             #
# Version: 0.9.4                                                                   #
#----------------------------------------------------------------------------------#

tiger <- function(data,
                  lambda = NULL,
                  nlambda = NULL,
                  lambda.min.ratio = NULL,
                  rho = NULL,
                  method = "clime",
                  alg = NULL,
                  sym = "or",
                  shrink = NULL,
                  prec = 1e-3,
                  max.ite = NULL,
                  wmat = NULL,
                  standardize = FALSE,
                  correlation = FALSE,
                  perturb = TRUE,
                  verbose = TRUE)
{
  if(method!="clime" && method!="slasso" && method!="adp") {
    cat("\"method\" must be either \"clime\", \"slasso\" or \"adp\" \n")
    return(NULL)
  }
  
  if(method=="clime") {
    if(is.null(alg)) alg="hybrid"
    if(alg!="cdadm" && alg!="ladm" && alg!="hybrid"){
      cat("\"alg\" must be either \"cdadm\", \"ladm\" or \"hybrid\" when method is \"clime\".\n")
      return(NULL)
    }
  } else {
    if(alg != "cdadm" || is.null(alg))
      cat(sprintf("\"alg\" is \"cdadm\" or leave it empty when method is \"%s\"", method))
  }
  
  n = nrow(data)
  d = ncol(data)
  maxdf = max(n,d)
  est = list()
  est$cov.input = isSymmetric(data)
  if(est$cov.input)
  {
    if(verbose) {
      cat("The input is identified as the covriance matrix.\n")
    }
    if(method=="slasso") {
      cat("The input for \"slasso\" cannot be covriance matrix.\n")
      return(NULL)
    }
    if(correlation)
      S = cov2cor(data)
    else
      S = data
  }
  if(!est$cov.input)
  {
    if(standardize)
      data = scale(data)
    
    if(correlation)
      S = cor(data)
    else
      S = cov(data)
    
    if(!correlation)
      S = S*(1-1/n)
  }
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda))
  {
    if(method == "slasso") {
      if(is.null(nlambda) || nlambda==1){
        nlambda = 1
        lambda=pi*sqrt(log(d)/(2*n))
      }
      else{
        lambda0=pi*sqrt(log(d)/(2*n))
        lambda = seq(lambda0,0.5*lambda0,length=nlambda)
      }
    }
    else {
      if(is.null(nlambda))
        nlambda = 10
      if(is.null(lambda.min.ratio))
        lambda.min.ratio = 0.4
      lambda.max = min(max(S-diag(diag(S))),-min(S-diag(diag(S))))
      lambda.min = lambda.min.ratio*lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      rm(lambda.max,lambda.min,lambda.min.ratio)
      gc()
    }
  }
  
  est$lambda = lambda
  est$nlambda = nlambda
  
  if(method == "clime" || method == "adp"){
    if(method == "clime"){
      if(is.null(max.ite)){
        if(alg=="cdadm") max.ite=1e2
        else max.ite=1e3
      }
    }
    if(is.null(rho))
      rho = sqrt(d)
    if (is.logical(perturb)) {
      if (perturb) { 
        eigvals = eigen(S, only.values=T)$values
        #perturb = max(max(max(eigvals) - d*min(eigvals), 0)/(d-1), 1/n)
        perturb = 1/sqrt(n)
      } else {
        perturb = 0
      }
    }
    S = S + diag(d)*perturb
    if(method == "clime"){
      if(is.null(shrink)) shrink=1.5
      if(alg=="cdadm")
        re_tiger = tiger.clime.cdadm(S, d, maxdf, lambda, rho, shrink, prec, max.ite)
      if(alg=="hybrid")
        re_tiger = tiger.clime.hadm(S, d, maxdf, lambda, rho, shrink, prec, max.ite)
      if(alg=="ladm")
        re_tiger = tiger.clime.ladm(S, d, maxdf, lambda, rho, shrink, prec, max.ite)
    }
    if(method == "adp"){
      if(is.null(shrink)) shrink=0
      if(is.null(max.ite)) max.ite=3e2
      if(is.null(wmat)){
        #         if(rankMatrix(S)<d)
        #           S = S + diag(d)*1/n
        wmat = 1/(abs(solve(S))+1/n)
      }
      re_tiger = tiger.adp.cdadm(S, wmat, n, d, maxdf, lambda, rho, shrink, prec, max.ite)
      est$wmat = wmat
    }
  }
  
  if(method == "slasso"){
    if(is.null(shrink)) shrink=0
    if(is.null(max.ite)) max.ite=1e2
    re_tiger = tiger.slasso.cdadm(data, n, d, maxdf, lambda, shrink, prec, max.ite)
  }
  est$ite = re_tiger$ite
  
  for(j in 1:d) {
    if(re_tiger$col_cnz[j+1]>re_tiger$col_cnz[j])
    {
      idx.tmp = (re_tiger$col_cnz[j]+1):re_tiger$col_cnz[j+1]
      ord = order(re_tiger$row_idx[idx.tmp])
      re_tiger$row_idx[idx.tmp] = re_tiger$row_idx[ord + re_tiger$col_cnz[j]]
      re_tiger$x[idx.tmp] = re_tiger$x[ord + re_tiger$col_cnz[j]]
    }
  }
  G = new("dgCMatrix", Dim = as.integer(c(d*nlambda,d)), x = as.vector(re_tiger$x[1:re_tiger$col_cnz[d+1]]),
          p = as.integer(re_tiger$col_cnz), i = as.integer(re_tiger$row_idx[1:re_tiger$col_cnz[d+1]]))
  
  est$beta = list()
  est$path = list()
  est$df = matrix(0,d,nlambda)
  est$rss = matrix(0,d,nlambda)  
  est$sparsity = rep(0,nlambda)  
  for(i in 1:nlambda) {
    est$beta[[i]] = G[((i-1)*d+1):(i*d),]
    est$path[[i]] = abs(est$beta[[i]])
    est$df[,i] = apply(sign(est$path[[i]]),2,sum)
    
    if(sym == "or")
      est$path[[i]] = sign(est$path[[i]] + t(est$path[[i]]))
    if(sym == "and")
      est$path[[i]] = sign(est$path[[i]] * t(est$path[[i]]))
    est$sparsity[i] = sum(est$path[[i]])/d/(d-1)
  }
  rm(G)
  est$icov = re_tiger$icov
  est$icov1 = re_tiger$icov1
  est$sigma = S
  est$data = data
  est$method = method
  est$sym = sym
  est$verbose = verbose
  est$standardize = standardize
  est$correlation = correlation
  est$perturb = perturb
  class(est) = "tiger"
  return(est)
}

print.tiger <- function(x, ...)
{  
  cat("\n tiger options summary: \n")
  cat(x$nlambda, " lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Method=", x$method, "\n")
  cat("Path length:",x$nlambda,"\n")
  cat("Graph dimension:",ncol(x$data),"\n")
  cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}

plot.tiger = function(x, align = FALSE, ...){
  gcinfo(FALSE)
  
  if(x$nlambda == 1)	par(mfrow = c(1, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  if(x$nlambda == 2)	par(mfrow = c(1, 3), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  if(x$nlambda >= 3)	par(mfrow = c(1, 4), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  
  if(x$nlambda <= 3)	z.final = 1:x$nlambda
  
  if(x$nlambda >=4){
    z.max = max(x$sparsity)
    z.min = min(x$sparsity)
    z = z.max - z.min
    z.unique = unique(c(which(x$sparsity>=(z.min + 0.03*z))[1],which(x$sparsity>=(z.min + 0.07*z))[1],which(x$sparsity>=(z.min + 0.15*z))[1]))
    
    
    if(length(z.unique) == 1){
      if(z.unique<(x$nlambda-1))	z.final = c(z.unique,z.unique+1,z.unique+2)
      if(z.unique==(x$nlambda-1)) z.final = c(z.unique-1,z.unique,z.unique+1)
      if(z.unique==x$nlambda) 	z.final = c(z.unique-2,z.unique-1,z.unique)
    }
    
    if(length(z.unique) == 2){
      if(diff(z.unique)==1){
        if(z.unique[2]<x$nlambda) z.final = c(z.unique,z.unique[2]+1) 
        if(z.unique[2]==x$nlambda) z.final = c(z.unique[1]-1,z.unique)
      }
      if(diff(z.unique)>1) z.final = c(z.unique[1],z.unique[1]+1,z.unique[2])
    }
    
    if(length(z.unique) == 3) z.final = z.unique
    
    rm(z.max,z.min,z,z.unique)
    gc()
    
  }
  plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
  
  lines(x$lambda[z.final],x$sparsity[z.final],type = "p")
  
  if(align){
    layout.grid = layout.fruchterman.reingold(graph.adjacency(as.matrix(x$path[[z.final[length(z.final)]]]), mode="undirected", diag=FALSE))
    for(i in z.final){
      g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
      plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
      rm(g)
      gc()
    }
    rm(layout.grid)
  }
  if(!align){
    for(i in z.final){
      g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
      layout.grid = layout.fruchterman.reingold(g)
      plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
      rm(g,layout.grid)
      gc()
    }
  }
  if(align) cat("Three plotted graphs are aligned according to the third graph\n")
}
