\name{LogGammaDist}
\alias{LogGammaDist}
\alias{dloggamma}
\alias{ploggamma}
\alias{qloggamma}
\alias{rloggamma}

\title{
The Extended Log Gamma dDistribution
}
\description{
Density, distribution function, quantile function and random
  generation for the Extended Log Gamma distribution with parameters
  \code{mu} (location), \code{sigma} (scale) and \code{lambda} (shape).
}
\usage{
dloggamma(x, mu=0, sigma=1, lambda)
ploggamma(q, mu=0, sigma=1, lambda)
qloggamma(p, mu=0, sigma=1, lambda)
rloggamma(n, mu=0, sigma=1, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{mu}{location parameter.}
  \item{sigma}{scale parameter.}
  \item{lambda}{shape parameter.}
}
% \details{
% TO BE DONE!
% }

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
