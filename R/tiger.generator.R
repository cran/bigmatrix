#-----------------------------------------------------------------------------------#
# Package: bigmatrix                                                                #
# tiger.generator(): Data generator                                                 #
# Author: Xingguo Li                                                                #
# Email: <xingguo.leo@gmail.com>                                                    #
# Date: July 29th 2012                                                              #
# Version: 0.9.2                                                                    #
#-----------------------------------------------------------------------------------#

## Main function
tiger.generator <- function(n = 200, d = 50, graph = "random", v = NULL, u = NULL, g = NULL, prob = NULL, vis = FALSE, verbose = TRUE){	
  gcinfo(FALSE)
  if(verbose) cat("Generating data from the multivariate normal distribution with the", graph,"graph structure....")
  if(is.null(g)){
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(d > 40)	g = ceiling(d/20)
      if(d <= 40) g = 2
    }
    if(graph == "block"){
      if(d > 100)  g = ceiling(d/20)
      if(d <= 100) g = 5
    }
  }
  
  if(graph == "random"){
    if(is.null(prob))	prob = min(1, 3/d)
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }
  
  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)	prob = 0.3
      if(d/g <= 30)	prob = min(1,6*g/d)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }  
  
  
  # parition variables into groups
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()
  
  # build the graph structure
  theta = matrix(0,d,d);
  if(graph == "band"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      diag(theta[1:(d-i),(1+i):d]) = 1
      diag(theta[(1+i):d,1:(d-1)]) = 1
    }	
  }
  if(graph == "cluster"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      tmp = which(g.ind==i)
      tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
      tmp2 = tmp2 + t(tmp2)		 	
      theta[tmp,tmp][tmp2<prob] = 1
      rm(tmp,tmp2)
      gc()
    }
  }
  if(graph == "hub"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      tmp = which(g.ind==i)
      theta[tmp[1],tmp] = 1
      theta[tmp,tmp[1]] = 1
      rm(tmp)
      gc()
    }
  }
  if(graph == "random"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    
    tmp = matrix(runif(d^2,0,0.5),d,d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    #theta[tmp >= tprob] = 0
    rm(tmp)
    gc()
  }
  
  if(graph == "scale-free"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    out = .C("SFGen",dd0=as.integer(2),dd=as.integer(d),G=as.integer(theta),package="tiger")
    theta = matrix(as.numeric(out$G),d,d)
  }
  if(graph=="band"||graph=="cluster"||graph=="hub"||graph=="random"||graph=="scale-free") {
    diag(theta) = 0
    omega = theta*v
    
    # make omega positive definite and standardized
    diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
    sigma = cov2cor(solve(omega))
    omega = solve(sigma)
  }
  # decay omega
  if (graph == "decay") {
    if(is.null(prob)) {
      prob = 0.6
    } else {
      if(prob<=0 || prob>=1) {
        cat("For \"decay\" model, the prob need to be between 0 and 1 as a base number\n")
        break
      }
    }
    omega = matrix(0,nrow=d,ncol=d)
    for (i in 1:d){
      for (j in 1:d){
        omega[i,j] = prob^abs(i-j)
      }
    }
    #omega[omega<1e-4] = 0
  }
  
  # dense omega
  if (graph == "dense") {
    omega = matrix(.5, ncol = d, nrow = d)
    omega[row(omega)==col(omega)] = 1
  }
  
  # random sparse omega
  if (graph == "sparse") {
    if(is.null(prob)) prob = 0.02
    while(TRUE){
      ndlt = 300
      rand_num = runif(d^2)
      rand_vec = rand_num
      rand_vec[rand_num>prob] = 0
      rand_vec[rand_num<=prob] = .5
      omega = matrix(rand_vec, ncol = d, nrow = d)
      ind = lower.tri(omega)
      omega[ind] = t(omega)[ind]
      dlt1 = seq(0.5001,sqrt(d),length.out=ndlt)
      for(j in 1:5){
        dlt=dlt1
        dif_cond = matrix(0,ncol=ndlt, nrow=1)
        for(i in 1:ndlt) {
          omega[row(omega)==col(omega)] = dlt[i]
          dif_cond[i] = abs(kappa(omega)-d)
        }
        idx_min = which(dif_cond == min(dif_cond))
        
        idx1=ifelse(idx_min[1]-1<1,idx_min[1],idx_min[1]-1)
        idx2=ifelse(idx_min[1]+1>ndlt,idx_min[1],idx_min[1]+1)
        dlt1 = seq(dlt[idx1],dlt[idx2],length.out=ndlt)
      }
      omega[row(omega)==col(omega)] = dlt[idx_min[1]]
      omega = omega/dlt[idx_min[1]]
      if(isPosDef(omega)) 
        break
    }
  }
  
  # block omega
  if (graph == "block") {
    omega = matrix(0,nrow=d,ncol=d)
    blocksize = g
    nblock = floor(d/blocksize)
    res = d%%blocksize
    for(i in 1:nblock) {
      omega[((i-1)*blocksize+1):(i*blocksize),((i-1)*blocksize+1):(i*blocksize)] = .5
    }
    if(res != 0) omega[(nblock*blocksize+1):d,(nblock*blocksize+1):d] = .5
    omega[row(omega)==col(omega)] = 1
    while(TRUE) {
      per = sample(1:dim(omega)[1])
      omega_tmp = matrix(data=omega[per,per],nrow=dim(omega)[1],ncol=dim(omega)[2])
      if(isPosDef(omega_tmp)) break
    }
    omega = omega_tmp
  }
  
  if(graph=="decay"||graph=="dense"||graph=="sparse"||graph=="block"){
    theta[which(omega>1e-3)] = 1
    theta[row(theta)==col(theta)] = 0
    sigma = solve(omega)
  }
  
  # generate multivariate normal data
  x = mvrnorm(n,rep(0,d),sigma)
  sigmahat = cor(x)
  
  # graph and covariance visulization
  if(vis == TRUE){
    fullfig = par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    fullfig[1] = image(theta, col = gray.colors(256),  main = "Adjacency Matrix")
    
    fullfig[2] = image(sigma, col = gray.colors(256), main = "Covariance Matrix")
    g = graph.adjacency(theta, mode="undirected", diag=FALSE)
    layout.grid = layout.fruchterman.reingold(g)
    
    fullfig[3] = plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")
    
    fullfig[4] = image(sigmahat, col = gray.colors(256), main = "Empirical Matrix")
    rm(fullfig,g,layout.grid)
    gc()
  }
  if(verbose) cat("done.\n")
  rm(vis,verbose)
  gc()
  
  sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega, theta = Matrix(theta,sparse = TRUE), sparsity= sum(theta)/(d*(d-1)), graph.type=graph)
  class(sim) = "sim" 
  return(sim)
}


print.sim = function(x, ...){
  cat("Simulated data generated by tiger.generator()\n")
  cat("Sample size: n =", nrow(x$data), "\n")
  cat("Dimension: d =", ncol(x$data), "\n")
  cat("Graph type = ", x$graph.type, "\n")
  cat("Sparsity level:", sum(x$theta)/ncol(x$data)/(ncol(x$data)-1),"\n")
}

plot.sim = function(x, ...){
  gcinfo(FALSE)	
  par = par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  image(as.matrix(x$theta), col = gray.colors(256),  main = "Adjacency Matrix")
  image(x$sigma, col = gray.colors(256), main = "Covariance Matrix")
  g = graph.adjacency(x$theta, mode="undirected", diag=FALSE)
  layout.grid = layout.fruchterman.reingold(g)
  
  plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")
  rm(g, layout.grid)
  gc()
  image(x$sigmahat, col = gray.colors(256), main = "Empirical Covariance Matrix")
}

isPosDef <- function(M) { 
  if (all(M == t(M) ) ) {  # first test  symmetricity
    if (all(eigen(M)$values > 0) ) {TRUE}
    else {FALSE} 
  } # not symmetric
  else {FALSE}
}
