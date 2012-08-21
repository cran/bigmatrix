\name{tiger-internal}
\alias{tiger.likelihood}
\alias{tiger.tracel2}
\alias{tiger.cv}
\alias{part.cv}
\alias{tiger.clime.cdadm}
\alias{tiger.clime.ladm}
\alias{tiger.slasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal tiger functions}
\description{
  Internal tiger functions
}
\usage{
tiger.likelihood(Sigma, Omega)
tiger.tracel2(Sigma, Omega)
tiger.cv(obj, loss=c("likelihood", "tracel2"), fold=5)
part.cv(n, fold)
tiger.clime.cdadm(Sigma, d, maxdf, lambda, rho, prec, max.ite, verbose)
tiger.clime.ladm(Sigma, d, maxdf, lambda, rho, prec, max.ite, verbose)
tiger.slasso(data, lambda, max.ite)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{Covariance matrix.}
  \item{Omega}{Inverse covariance matrix.}
  \item{obj}{An object with S3 class returned from \code{"tiger"}.}
  \item{loss}{Type of loss function for cross validation.}
  \item{fold}{The number of fold for cross validatio.}
  \item{n}{The number of data.}
  \item{d}{Dimension of data.}
  \item{maxdf}{Maximal degree of freedom.}
  \item{lambda}{Grid of non-negative values for the constraint parameter lambda.}
  \item{rho}{Value of augmented Lagrangian multipiler.}
  \item{prec}{Stopping criterion.}
  \item{max.ite}{Maximal value of iterations.}
  \item{verbose}{Tracing information printing.}
  \item{data}{\code{n} by \code{d} data matrix. }
}

\details{
  These are not intended for use by users.   
}

%\references{
%}

\author{
Xingguo Li, Tuo Zhao, Lie Wang and Han Liu \cr
Maintainer: Xingguo Li <xingguo.leo@gmail.com>
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tiger}} and \code{\link{bigmatrix-package}}.
}
