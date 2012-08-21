\name{tiger}
\alias{tiger}

\title{
Tuning Insensitive Graph Estimation and Regression
}

\description{
Tiger provides two estimation procedures for high dimensional precision matrix estimation based column by column scheme: (1) Tuning Insensitive Graph Estimation and Regression (TIGER). (2) The Constrained L1 minimization for sparse precision Matrix Estimation (CLIME). Both TIGER and CLIME have better performance than graphical lasso in parameter estimation under the spectral norm. Moreover, TIGER is less sensitive to tuning parameter than graphical lasso. The algorithm is implemented based on Alternative Direction Method combined with coordinate descent, which scales up to large problems. The computation is memory-optimized using the sparse matrix output.
}

\usage{
tiger(data, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
     rho = NULL, method = "clime", alg = "cdadm", sym = "or", 
     prec = 1e-3, max.ite = NULL, standardize = FALSE, 
     correlation = FALSE, perturb = TRUE, verbose = TRUE)
}

\arguments{
  \item{data}{There are 2 options for \code{"clime"}: (1) \code{data} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension). For \code{"slasso"}, covariance input is not supported. }
  \item{lambda}{A sequence of decresing positive numbers to control the regularization. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. Default value is \eqn{ (log(d)/(2n))^{1/2}} for \code{"slasso"}.}
  \item{nlambda}{Number of the constraint parameter \code{lambda}. Default value is 10 for \code{"clime"} and 1 for \code{"slasso"}.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda} for \code{clime}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.4}.}
  \item{rho}{Penalization parameter used in alternative direction method for \code{clime}. The default value is \eqn{\sqrt(d)}. }
  \item{method}{CLIME is applied if \code{method="clime"} and TIGER is applied if \code{method="slasso"}. Default value is \code{"clime"}.}
  \item{alg}{Algorithms of \code{"clime"} with two options: \code{"cdadm"} and \code{"ladm"}. Default value is \code{"cdadm"}.}
  \item{sym}{Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}.}
  \item{prec}{Stopping criterion. The default value is 1e-3.}
  \item{max.ite}{Maximal value of iterations. The default value is 1e3 if \code{method="clime"} and \code{alg="ladm"}, and is 1e2 if \code{method="clime"} and \code{alg="cdadm"}, or \code{method="slasso"}.}
  \item{standardize}{Variables are standardized to have mean zero and unit standard deviation if \code{standardize=TRUE}. The default value is \code{FALSE}.}
  \item{correlation}{Correlation matrix is used as the input of \code{Sigma} for \code{mathod="clime"} if \code{correlation=TRUE}. The default value is \code{FALSE}.}
  \item{perturb}{The diagonal of \code{Sigma} is added by a positive value to guarantee that \code{Sigma} is positive definite if \code{perturb=TRUE}. The default value is FALSE.}
  \item{verbose}{Tracing information printing is disabled if \code{verbose = FALSE}. The default value is \code{TRUE}.}
  }
\details{ 
All methods proceed column-wise estimation of the precision matrix, which allows for parallel computation of large scale of data. For \code{"clime"} based algorithms \code{"cdadm"} and \code{"ladm"}, they solve the following objected function\cr
\deqn{
    \min || \Omega ||_1 \quad \textrm{s.t. } || S \Omega - I ||_\infty \le \lambda
  }\cr
where \eqn{||\cdot||_1} is element-wise 1-norm of a matrix and \eqn{||\cdot||_\infty} is element-wise \eqn{\infty}-norm of a matrix. For \code{"cdadm"}, it solves \code{"clime"} via coordinate descent implementation of alternative direction method of multipliers, and for \code{"ladm"} the solver is implemented via linearized alternative direction method. We also implemented active set of features and warm start for all methods for further acceleration. 
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
  Sequence of regularization parameters \code{lambda} used in the program.
  }
  \item{nlambda}{
  Number of values for program generated \code{lambda}.
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
  If \code{method="clime"}, it is a list of two matrices where ite[[1]] is the number of external iteration and ite[[2]] is the number of internal iteration with the entry of (i,j) as the number of iteration of i-th column and j-th lambda. If \code{method="slasso"}, it is a matrix of iteration with the entry of (i,j) as the number of iteration of i-th column and j-th lambda. 
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
1. J. Friedman, T. Hastie and R. Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2010. \cr
2. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011. 
3. H. Liu and L. Wang. TIGER: A tuning-insenstiive approach for optimally estimating large undirected graphs. \emph{Technical Report}, 2012.
}

\seealso{
\code{\link{bigmatrix-package}}, \code{\link{tiger.generator}}, \code{\link{tiger.select}}, \code{\link{tiger.plot}}, \code{\link{tiger.roc}}, \code{\link{plot.tiger}}, \code{\link{plot.select}}, \code{\link{plot.roc}}, \code{\link{plot.sim}}, \code{\link{print.tiger}}, \code{\link{print.select}}, \code{\link{print.roc}} and \code{\link{print.sim}}.
}

%\note{
%}

\examples{
## generate data
n = 100
d = 50
D = tiger.generator(n=n,d=d,graph="random", prob=0.03)
plot(D)

## sparse precision matrix estimation with method "clime" and algorithm "ladm"
out1 = tiger(D$data, method = "clime", alg="ladm")
plot(out1)
out1.select = tiger.select(out1, criterion="stars", stars.thresh = 0.05)
plot(out1.select)
tiger.plot(out1$path[[7]])

## sparse precision matrix estimation with method "slasso"
out2 = tiger(D$data, method = "slasso", nlambda = 5)
plot(out2)
tiger.plot(out2$path[[3]])
}
