\name{tiger.generator}
\alias{tiger.generator}

\title{
Data generator
}

\description{
Implements the data generation from multivariate normal distributions with different graph structures, including \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"band"}, \code{"scale-free"}, \code{"decay"}, \code{"sparse"} and \code{"block"}.
}

\usage{
tiger.generator(n = 200, d = 50, graph = "random", v = NULL, u = NULL, 
g = NULL, prob = NULL, vis = FALSE, verbose = TRUE)
}

\arguments{
  \item{n}{
The number of observations (sample size). The default value is \code{200}.
}
  \item{d}{
The number of variables (dimension). The default value is \code{50}.
}
  \item{graph}{
The graph structure with 8 options: \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"band"}, \code{"scale-free"}, \code{"decay"}, \code{"sparse"} and \code{"block"}.
}
  \item{v}{
The off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with \code{u}. The default value is \code{0.3}.
}
  \item{u}{
A positive number being added to the diagonal elements of the precision matrix, to control the magnitude of partial correlations. The default value is \code{0.1}.
}
  \item{g}{
For \code{"cluster"} or \code{"hub"} graph, \code{g} is the number of hubs or clusters in the graph. The default value is about \code{d/20} if \code{d >= 40} and \code{2} if \code{d < 40}. For \code{"band"} graph, \code{g} is the bandwidth and the default value is \code{1}. For \code{"block"} graph, \code{g} is the size of the block. NOT applicable to \code{"random"} graph.
}
  \item{prob}{
For \code{"random"} graph, it is the probability that a pair of nodes has an edge. The default value is \code{3/d}. For \code{"cluster"} graph, it is the probability that a pair of nodes has an edge in each cluster. The default value is \code{6*g/d} if \code{d/g <= 30} and \code{0.3} if \code{d/g > 30}. NOT applicable to \code{"hub"} or \code{"band"} graphs.
}
  \item{vis}{
Visualize the adjacency matrix of the true graph structure, the graph pattern, the covariance matrix and the empirical covariance matrix. The default value is \code{FALSE}
}
  \item{verbose}{
If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
}
}

\details{
Given the adjacency matrix \code{theta}, the graph patterns are generated as below:\cr\cr
(I) \code{"random"}: Each pair of off-diagonal elements are randomly set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with probability \code{prob}, and \code{0} other wise. It results in about \code{d*(d-1)*prob/2} edges in the graph.\cr\cr
(II)\code{"hub"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each group is associated with a "center" row \code{i} in that group. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} if \code{j} also belongs to the same group as \code{i} and \code{0} otherwise. It results in \code{d - g} edges in the graph.\cr\cr
(III)\code{"cluster"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with the probability \code{prob}if both \code{i} and \code{j} belong to the same group, and \code{0} other wise. It results in about \code{g*(d/g)*(d/g-1)*prob/2} edges in the graph.\cr\cr
(IV)\code{"band"}: The off-diagonal elements are set to be \code{theta[i,j]=1} if \code{1<=|i-j|<=g} and \code{0} other wise. It results in \code{(2d-1-g)*g/2} edges in the graph.\cr\cr
(V) \code{"scale-free"}: The graph is generated using B-A algorithm. The initial graph has two connected nodes and each new node is connected to only one node in the existing graph with the probability proportional to the degree of the each node in the existing graph. It results in \code{d} edges in the graph.\cr\cr
(VI) \code{"decay"}: \eqn{\omega^0_{ij} = 0.6^{| i-j |}}. \cr\cr
(VII) \code{"sparse"}: Let \eqn{\Omega_0 = O + \delta I}, where each off-diagonal entry in \eqn{O} is generated independently and equals to 0.5 with probability \code{prob} or 0 with probability 1-\code{prob}. \eqn{\delta} is chosen such that the conditional number (the ratio of maximal and minimal singular values of a matrix) is equal to d. Finally, the matrix is standardized to have unit diagonals. The default value of \code{prob} = 0.1. \cr\cr
(VIII) \code{"block"}: A block diagonal matrix with block size \code{g} where each block has off-diagonal entries equal to 0.5 and diagonal 1. The resulting matrix is then randomly permuted. Teh default value block size is 5. \cr\cr

The adjacency matrix \code{theta} has all diagonal elements equal to \code{0}. To obtain a positive definite precision matrix, the smallest eigenvalue of \code{theta*v} (denoted by \code{e}) is computed. Then we set the precision matrix equal to \code{theta*v+(|e|+0.1+t)I}. The covariance matrix is then computed to generate multivariate normal data.
}

\value{
An object with S3 class "sim" is returned:
\item{data}{
The \code{n} by \code{d} matrix for the generated data
}
\item{sigma}{
The covariance matrix for the generated data
}
  \item{omega}{
The precision matrix for the generated data
}
  \item{sigmahat}{
The empirical covariance matrix for the generated data
}
  \item{theta}{
The adjacency matrix of true graph structure (in sparse matrix representation) for the generated data
}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang and Han Liu \cr
Maintainers: Xingguo Li<xingguo.leo@gmail.com>
}

%\references{
%}
\seealso{
\code{\link{tiger}} and \code{\link{bigmatrix-package}}
}
\examples{
## band graph with bandwidth 3
L = tiger.generator(graph = "band", g = 3)
plot(L)

## random sparse graph
L = tiger.generator(vis = TRUE)

## random dense graph
L = tiger.generator(prob = 0.5, vis = TRUE)

## hub graph with 6 hubs
L = tiger.generator(graph = "hub", g = 6, vis = TRUE)

## hub graph with 8 clusters
L = tiger.generator(graph = "cluster", g = 8, vis = TRUE)

## scale-free graphs
L = tiger.generator(graph="scale-free", vis = TRUE)

## decay graphs
L = tiger.generator(graph="decay", vis = TRUE)

## block graphs
L = tiger.generator(graph="block", g=6, vis = TRUE)

## sparse graphs
L = tiger.generator(graph="sparse", prob=0.05, vis = TRUE)
}
