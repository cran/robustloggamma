\name{summary.loggammarob}
\alias{summary.loggammarob}
\alias{print.summary.loggammarob}

\title{
Summary Method for \code{"loggammarob"} Objects
}
\description{
  Summary method for \R object of class \code{"loggammarob"} and
  \code{\link{print}} method for the summary object.
}
\usage{
\method{summary}{loggammarob}(object, p = NULL, conf.level = 0.95,
  prob = 1e-05, ...)
\method{print}{summary.loggammarob}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{
    an \R object of class \code{loggammarob}, typically created by
    \code{\link{loggammarob}}.
  }  
  \item{p}{
    numeric. Values in the interval [0,1]. Quantile orders for
    which point estimation and confidence interval are require.   
  }
  \item{conf.level}{
    numeric. A scalar or a vector of length \code{4+length(p)}. This is
    the confidence level used to construct confidence intervals for the
    four parameters \code{mu}, \code{sigma}, \code{lambda} and
    \code{eta} and for the quantiles.
  }
  \item{prob}{
    numeric. Value in the interval [0,1]. This is used to determine the
    interval of numerical integration in the evaluation of the
    asymptotic variance and covariance matrix. See details below.
  }
  \item{x}{
    an \R object of class \code{summary.loggammarob}, typically
    resulting from \code{summary(\link{loggammarob}(..),..)}.
  }
  \item{digits}{
    number of digits for printing, see \code{digits} in
    \code{\link{options}}.
  }
  \item{\dots}{
    potentially more arguments passed to methods.
  }
}
\details{
  The \code{prob} argument determines the interval of the numerical
  integration in the evaluation of the asymptotic variance and
  covariance matrix with the following code
  \code{qloggamma(p=prob/2,lambda=lambda)} for the lower limit and
  \code{qloggamma(p=1-prob/2,lambda=lambda)} for the upper limit. Here
  \code{lambda} is the estimate.
}
\value{
 An object of class \code{summary.loggammarob}. A list that includes the
  following components:
  \item{muse}{standard error for the \code{mu} estimate.}
  \item{sigmase}{standard error for the \code{sigma} estimate.}
  \item{lambdase}{standard error for the \code{lambda} estimate.}
  \item{etase}{standard error for the \code{eta} estimate.}
  \item{muconf.int}{confidence iterval for the \code{mu} parameter.}
  \item{sigmaconf.int}{confidence iterval for the \code{sigma} parameter.}
  \item{lambdaconf.int}{confidence iterval for the \code{lambda} parameter.}
  \item{etaconf.int}{confidence iterval for the \code{eta} parameter.}

  If \code{p} is not \code{NULL} then

  \item{q}{quantiles estimate.}
  \item{qse}{standard error for the require quantiles estimate.}
  \item{qconf.int}{confidence iterval for the require quantiles parameter.}  
}
\references{
  C. Agostinelli, A. Marazzi and V.J. Yohai (2013) Robust estimates of
  the generalized loggamma distribution, Technometrics, Accepted.
}
\author{
 C. Agostinelli, A. Marazzi, V.J. Yohai and A. Randriamiharisoa
}

\seealso{
\code{\link{loggammarob}}
}
\examples{
  set.seed(1234)
  x <- rloggamma(n=50, lambda=1)
  res <- loggammarob(x, control=loggammarob.control(lower=0, upper=2, n=30))
  summary(res)
}

\keyword{robust}
\keyword{models}
