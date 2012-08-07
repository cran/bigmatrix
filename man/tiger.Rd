\name{tiger}
\alias{tiger}

\title{
Tuning Insensitive Graph Estimation and Regression
}

\description{
Solve for the precision matrix estimations with two different methods. The first is the tuning insensitive graph estimation and regression implemented with scaled Lasso. Another method is by estimating at a grid of values for the constrained lambda and the augmented Lagrangian multiplier rho based on constrained \eqn{\ell_1} minimization for sparse precision matrix estimation, which guarantees positive definite output by projected augmented Lagrangian method and further accelerates the algorithm via coordinate descent and active set methods.
}

\usage{
tiger(data, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
     rho = NULL, method = "clime", sym = "or", standardize = FALSE, 
     correlation = FALSE, biased = TRUE, perturb = TRUE, 
     verbose = TRUE, prec = 1e-3)
}

\arguments{
  \item{data}{There are 2 options for \code{"clime"}: (1) \code{data} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension). For \code{"slasso"}, covariance input is not supported. }
  \item{lambda}{Grid of non-negative values for the constraint parameter \code{lambda}. If missing, \code{nlambda}  values from \code{lambda.min} to \code{lambda.max} will be generated. Default value is \eqn{ (log(p)/(2n))^{1/2}} for \code{"slasso"}.}
  \item{nlambda}{Number of the constraint parameter \code{lambda}. Default value is 10 for \code{"clime"} and 1 for \code{"slasso"}.}
  \item{lambda.min.ratio}{The smallest value for \code{lambda} for \code{clime}, as a fraction of the uppperbound (\code{MAX}) of the regularization parameter. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. The default value is \code{0.2}.}
  \item{rho}{Augmented Lagrangian multiplier for \code{clime}. The default value is the square root of the dimension of data. }
  \item{method}{Precision matrix estimation with two options: \code{"clime"} and \code{"slasso"}. Default value is \code{"clime"}.}
  \item{sym}{Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}.}
  \item{standardize}{Whether the variables will be standardized to have mean zero and unit standard deviation. The default value is \code{FALSE}.}
  \item{correlation}{Whether to use the correlation matrix rather than the covariance matrix as the \code{Sigma} for \code{clime}.}
  \item{biased}{Whether the sample covariance matrix is biased estimator of covariance matrix or not for \code{clime}. It is only applicable when \code{correlation=FALSE}. The default value is \code{TRUE}.}
  \item{perturb}{Whether a perturbed \code{Sigma} should be used or the positive perturbation added if it is numerical for \code{clime}. The default value is TRUE.}
  \item{verbose}{If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.}
  \item{prec}{Stopping criterion. The default value is 1e-3.}
  }
\details{ 
A constrained \eqn{\ell_1} minimization approach for sparse precision matrix estimation solved by projected augmented Lagrangian method. It solves the following objected function\cr
\deqn{
    \min | \Omega |_1 \quad \textrm{s.t. } || S
  \Omega - I ||_\infty \le \lambda, \tau_{\ell} \le \Lambda_{min}(\Omega) \le \Lambda_{max}(\Omega) \le \tau_{u}
  }\cr
}
\value{
An object with S3 class \code{"tiger"} is returned:  
  \item{data}{
  The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input
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
  Sequence of estimated omega, corresponding to regularization parameters.
  }
  \item{sym}{
  The \code{sym} from the input.
  }
  \item{method}{
  Whether \code{"clime"} or \code{"slasso"} is used as the solver. 
  }
  \item{path}{
  A list of \code{d} by \code{d} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
  }
  \item{sparsity}{
  The sparsity levels of the graph path.
  }
  \item{ite}{
  For \code{"clime"}, it is a list of two matrices where ite[[1]] is the number of external iteration and ite[[2]] is the number of internal iteration with the entry of (i,j) as the number of iteration of i-th column of j-th lambda. No output for \code{"slasso"}. 
  }
  \item{df}{
  It is a \code{d} by \code{nlambda} matrix. Each row contains the number of nonzero coefficients along the lasso solution path. 
  }
  \item{standardize}{Flag of whether the variables are standardized to have mean zero and unit standard deviation.}
  \item{correlation}{Flag of whether the correlation matrix rather than the covariance matrix is used as the \code{Sigma}.}
  \item{biased}{Flag of whether the sample covariance matrix is biased estimator of covariance matrix or not. }
  \item{perturb}{Flag of whether a perturbed \code{Sigma} should be used or the positive perturbation added if it is numerical for \code{clime}.}
  \item{verbose}{
  Flag of whether tracing information is printed.
  }
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang and Han Liu \cr
Maintainers: Xingguo Li <xingguo.leo@gmail.com>
}

\references{
1. Jerome Friedman, Trevor Hastie and Rob Tibshirani. \emph{Regularization Paths for Generalized Linear Models via Coordinate Descent.} Journal of Statistical Software, 2010.\cr
2. Tony Cai, Weidong Liu and Xi Luo. \emph{A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation.} Journal of the American Statistical Association, 2011.
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
D = tiger.generator(n=n,d=d,graph="random", prob=0.02)
plot(D)

## sparse precision matrix estimation with method "clime"
out1 = tiger(D$data, method = "clime")
out1.select = tiger.select(out1, criterion="cv")
plot(out1)
plot(out1.select)
tiger.plot(out1$path[[10]])

## sparse precision matrix estimation with method "slasso"
out2 = tiger(D$data, method = "slasso", nlambda = 5)
plot(out2)
tiger.plot(out2$path[[2]])
}
