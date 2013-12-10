\name{LogGammaDist}
\alias{LogGammaDist}
\alias{dloggamma}
\alias{ploggamma}
\alias{qloggamma}
\alias{rloggamma}

\title{
The Extended Log Gamma Distribution
}
\description{
Density, distribution function, quantile function and random
  generation for the Extended Log Gamma distribution with parameters
  \code{mu} (location), \code{sigma} (scale) and \code{lambda} (shape).
}
\usage{
dloggamma(x, mu=0, sigma=1, lambda, log=FALSE)
ploggamma(q, mu=0, sigma=1, lambda)
qloggamma(p, mu=0, sigma=1, lambda)
rloggamma(n, mu=0, sigma=1, lambda)
}

\arguments{
  \item{x,q}{numeric. Vector of quantiles.}
  \item{p}{numeric. Vector of probabilities.}
  \item{n}{numeric. Number of observations.}
  \item{mu}{numeric. Location parameter.}
  \item{sigma}{numeric. Scale parameter.}
  \item{lambda}{numeric. Shape parameter.}
  \item{log}{logical. If \code{TRUE}, the log of the density is reported.}
}

\value{
  \code{dloggamma} gives the density, \code{ploggamma} gives the distribution
  function, \code{qloggamma} gives the quantile function, and \code{rloggamma}
  generates random deviates.
}

\seealso{
  \code{\link{GammaDist}}
}

\examples{
  set.seed(1234)
  x <- rloggamma(10, lambda=1)
  x
}

\keyword{distribution}
