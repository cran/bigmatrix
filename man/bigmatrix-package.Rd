\name{bigmatrix-package}
\alias{bigmatrix-package}
\alias{bigmatrix}
\docType{package}
\title{
TIGER: Tuning-Insensitive Graph Estimation and Regression
}
\description{
The package "bigmatrix" estimates high dimensional precision matrices for Gaussian graphical models. We adopt three estimation procedures based on column by column regression: (1) Tuning-Insensitive Graph Estimation and Regression (TIGER); (2) The Constrained L1 Minimization for Sparse Precision Matrix Estimation (CLIME); (3) Adaptive CLIME for Sparse Precision Matrix Estimation (ADP). The optimization algorithms for all three methods are implemented based on Alternative Direction Method combined with coordinate descent. The computation is memory-optimized using the sparse matrix output.
}
\details{
\tabular{ll}{
Package: \tab bigmatrix\cr
Type: \tab Package\cr
Version: \tab 0.9.6\cr
Date: \tab 2012-10-16\cr
License: \tab GPL-2\cr
}
}

\author{
Xingguo Li, Tuo Zhao, Lie Wang and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}
\references{
1. T. Cai, W. Liu and X. Luo. A constrained \eqn{\ell_1} minimization approach to sparse precision matrix estimation. \emph{Journal of the American Statistical Association}, 2011. 
}

\seealso{
\code{\link{tiger}}, \code{\link{tiger.generator}}, \code{\link{tiger.select}}, \code{\link{tiger.plot}}, \code{\link{tiger.roc}}.
}
