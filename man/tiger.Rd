\name{tiger}
\alias{tiger}

\title{
Tuning Insensitive Graph Estimation and Regression
}

\description{
The package "bigmatrix" estimates high dimensional precision matrices for Gaussian graphical models. We adopt three estimation procedures based on column by column regression: (1) Tuning-Insensitive Graph Estimation and Regression (TIGER); (2) The Constrained L1 Minimization for Sparse Precision Matrix Estimation (CLIME); (3) Adaptive CLIME for Sparse Precision Matrix Estimation (ADP). The optimization algorithms for all three methods are implemented based on Alternative Direction Method combined with coordinate descent. The computation is memory-optimized using the sparse matrix output.
}

\usage{
tiger(data, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
     rho = NULL, method = "clime", alg = NULL, sym = "or", 
     shrink=NULL, prec = 1e-3, max.ite = NULL, wmat = NULL, 
     standardize = FALSE, correlation = FALSE, perturb = TRUE, 
     verbose = TRUE)
}

\arguments{
  \item{data}{There are 2 options for \code{"clime"}: (1) \code{data} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension). For \code{"slasso"}, covariance input is not supported. }
  \item{lambda}{A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is \eqn{ (log(d)/(2n))^{1/2}} for \code{"slasso"}.}
  \item{nlambda}{The number of values used in \code{lambda}. Default value is 10 for \code{"clime"} and 1 for \code{"slasso"}.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda} for \code{clime}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.4}.}
  \item{rho}{Penalty parameter used in the optimization algorithm for \code{clime}. The default value is \eqn{\sqrt{d}}.}
  \item{method}{CLIME is applied if \code{method="clime"}, TIGER is applied if \code{method = "slasso"} and Adaptive Lasso is applied if \code{method = "adp"}. Default value is \code{"clime"}.}
  \item{alg}{Optimization algorithms to solve subproblems if and only if \code{method = "clime"}. The coordinate descent is applied if \code{"alg = cdadm"}, and the linearization is applied if \code{"alg = ladm"}. The combination of coordinate descent and linearization is applied if \code{"alg = hybrid"}. Default value is \code{"hybrid"} for \code{method = "clime"}.}
  \item{sym}{Symmetrization of output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}.}
  \item{shrink}{Shrinkage of regularization parameter based on precision of estimation. The default value is 1.5 if \code{method = "clime"} and the default value is 0 if \code{method="slasso"} or \code{method = "adp"}.}
  \item{prec}{Stopping criterion. The default value is 1e-3.}
  \item{max.ite}{The iteration limit. The default value is 1e3 if \code{method = "clime"} and \code{alg = "ladm"}, and is 1e2 if \code{method = "clime"} and \code{alg = "cdadm"}, or \code{method = "slasso"}.}
  \item{wmat}{Weight matrix for adaptive Lasso. It is only used when \code{method = "adp"}.}
  \item{standardize}{Variables are standardized to have mean zero and unit standard deviation if \code{standardize = TRUE}. The default value is \code{FALSE}.}
  \item{correlation}{Correlation matrix is used as the input of \code{Sigma} for \code{method = "clime"} if \code{correlation = TRUE}. The default value is \code{FALSE}.}
  \item{perturb}{The diagonal of \code{Sigma} is added by a positive value to guarantee that \code{Sigma} is positive definite if \code{perturb = TRUE}. User can specify a numeric value for \code{perturbe}. The default value is TRUE.}
  \item{verbose}{Tracing information printing is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{ 
CLIME solves the following objected function\cr
\deqn{
    \min || \Omega ||_1 \quad \textrm{s.t. } || S \Omega - I ||_\infty \le \lambda
  }\cr
where \eqn{||\cdot||_1} and \eqn{||\cdot||_\infty} are element-wise 1-norm and \eqn{\infty}-norm respectively. Our implementation is based on the alternative direction method, and is further accelertated by active set and warm start tricks. 
}
\value{
An object with S3 class \code{"tiger"} is returned:  
  \item{data}{
  The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input.
  }
  \item{cov.input}{
  An indicator of the sample covariance. 
  }
  \item{lambda}{
  The sequence of regularization parameters \code{lambda} used in the program.
  }
  \item{nlambda}{
  The number of values used in \code{lambda}.
  }
  \item{icov}{
  A list of \code{d} by \code{d} precision matrices corresponding to regularization parameters.
  }
  \item{sym}{
  The \code{sym} from the input.
  }
  \item{method}{
  The \code{method} from the input. 
  }
  \item{path}{
  A list of \code{d} by \code{d} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
  }
  \item{sparsity}{
  The sparsity levels of the graph path.
  }
  \item{ite}{
  If \code{method = "clime"}, it is a list of two matrices where ite[[1]] is the number of external iteration and ite[[2]] is the number of internal iteration with the entry of (i,j) as the number of iteration of i-th column and j-th lambda. If \code{method="slasso"}, it is a matrix of iteration with the entry of (i,j) as the number of iteration of i-th column and j-th lambda. 
  }
  \item{df}{
  It is a \code{d} by \code{nlambda} matrix. Each row contains the number of nonzero coefficients along the lasso solution path. 
  }
  \item{standardize}{The \code{standardize} from the input.}
  \item{correlation}{The \code{correlation} from the input.}
  \item{perturb}{The \code{perturb} from the input.}
  \item{verbose}{The \code{verbose} from the input.}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011. 
}

\seealso{
\code{\link{bigmatrix-package}}, \code{\link{tiger.generator}}, \code{\link{tiger.select}}, \code{\link{tiger.plot}}, \code{\link{tiger.roc}}, \code{\link{plot.tiger}}, \code{\link{plot.select}}, \code{\link{plot.roc}}, \code{\link{plot.sim}}, \code{\link{print.tiger}}, \code{\link{print.select}}, \code{\link{print.roc}} and \code{\link{print.sim}}.
}


\examples{

## generating data
n = 100
d = 100
D = tiger.generator(n=n,d=d,graph="hub",g=10)
plot(D)

## sparse precision matrix estimation with method "clime"
out1 = tiger(D$data, method = "clime")
plot(out1)
tiger.plot(out1$path[[7]])

## sparse precision matrix estimation with method "slasso"
out2 = tiger(D$data, method = "slasso", nlambda = 5)
plot(out2)
tiger.plot(out2$path[[3]])

## sparse precision matrix estimation with method "adp"
out3 = tiger(D$data, lambda=seq(0.25, 0.1, length.out=10), 
             method = "adp", wmat=1/(abs(D$omega)+1/n))
plot(out3)
tiger.plot(out3$path[[10]])
}
