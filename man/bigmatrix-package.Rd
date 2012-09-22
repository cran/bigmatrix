\name{bigmatrix-package}
\alias{bigmatrix-package}
\alias{bigmatrix}
\docType{package}
\title{
TIGER: Tuning-Insensitive Graph Estimation and Regression
}
\description{
The package "bigmatrix" estimates high dimensional precision matrices for Gaussian graphical models. We adopt two estimation procedures based column by column regression: (1) Tuning-Insensitive Graph Estimation and Regression (TIGER). (2) The Constrained L1 Minimization for Sparse Precision Matrix Estimation (CLIME). The optimization algorithm for CLIME is based on alternative direction method combined with coordinate decent. The computation is memory-optimized using the sparse matrix output.
}
\details{
\tabular{ll}{
Package: \tab bigmatrix\cr
Type: \tab Package\cr
Version: \tab 0.9.5\cr
Date: \tab 2012-08-23\cr
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
