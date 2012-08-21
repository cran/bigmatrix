\name{bigmatrix-package}
\alias{bigmatrix-package}
\docType{package}
\title{
Tuning Insensitive Graph Estimation and Regression
}
\description{
The package "bigmatrix" targets at high dimensional precision matrix estimation problem. We adopt two estimation procedures based column by column scheme: (1) Tuning Insensitive Graph Estimation and Regression (TIGER). (2) The Constrained L1 minimization for sparse precision Matrix Estimation (CLIME). Both TIGER and CLIME have better performance than graphical lasso in parameter estimation under the spectral norm. Moreover, TIGER is less sensitive to tuning parameter than graphical lasso. The algorithm is implemented based on Alternative Direction Method combined with coordinate descent, which scales up to large problems. The computation is memory-optimized using the sparse matrix output.
}
\details{
\tabular{ll}{
Package: \tab bigmatrix\cr
Type: \tab Package\cr
Version: \tab 0.9.4\cr
Date: \tab 2012-08-21\cr
License: \tab GPL-2\cr
}
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

\keyword{ package }
\seealso{
\code{\link{tiger}}, \code{\link{tiger.generator}}, \code{\link{tiger.select}}, \code{\link{tiger.plot}}, \code{\link{tiger.roc}},\code{\link{plot.tiger}}, \code{\link{plot.select}}, \code{\link{plot.roc}},\code{\link{plot.sim}}, \code{\link{print.tiger}}, \code{\link{print.select}}, \code{\link{print.roc}} and \code{\link{print.sim}}.
}
