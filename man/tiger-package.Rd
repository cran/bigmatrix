\name{bigmatrix-package}
\alias{bigmatrix-package}
\docType{package}
\title{
Tuning Insensitive Graph Estimation and Regression
}
\description{
Solve for the precision matrix estimations with two different methods. The first is the tuning insensitive graph estimation and regression implemented with scaled Lasso. Another method is by estimating at a grid of values for the constrained lambda and the augmented Lagrangian multiplier rho based on constrained L1 minimization for sparse precision matrix estimation, which guarantees positive definite output by projected augmented Lagrangian method and further accelerates the algorithm via coordinate descent and active set methods.
}
\details{
\tabular{ll}{
Package: \tab tiger\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2012-07-28\cr
License: \tab GPL-2\cr
}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang and Han Liu \cr
Maintainers: Xingguo Li<xingguo.leo@gmail.com>
}
\references{
1. Jerome Friedman, Trevor Hastie and Rob Tibshirani. \emph{Regularization Paths for Generalized Linear Models via Coordinate Descent.} Journal of Statistical Software, 2010. \cr
2. Tony Cai, Weidong Liu and Xi Luo. \emph{A constrained L1 minimization approach to sparse precision matrix estimation.} Journal of the American Statistical Association, 2011. 
}

\keyword{ package }
\seealso{
\code{\link{tiger}}, \code{\link{tiger.generator}}, \code{\link{tiger.select}}, \code{\link{tiger.plot}}, \code{\link{tiger.roc}},\code{\link{plot.tiger}}, \code{\link{plot.select}}, \code{\link{plot.roc}},\code{\link{plot.sim}}, \code{\link{print.tiger}}, \code{\link{print.select}}, \code{\link{print.roc}} and \code{\link{print.sim}}.
}
