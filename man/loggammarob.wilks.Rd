\name{loggammarob.wilks}
\alias{loggammarob.wilks}
\title{
  Robust Wilks test for generalized log gamma model.
}
\description{
  Performs a robust Wilks test to check equality between scale and
  shape parameters of a generalized log gamma model.
}
\usage{
loggammarob.wilks(x, thetainit = NULL, method = "L-BFGS-B",
  lower = c(-Inf, 1e-04), upper = c(Inf, Inf), ...)
}

\arguments{
  \item{x}{
    list. An object of class 'loggammarob', typically created by
    \code{\link{loggammarob}}.
  }
  \item{thetainit}{
    \code{NULL} or numeric. A vector of length 2 to be used as initial
    values, first element corresponds to mu, second element to
    sigma/lambda under the Null Hypothesis that sigma is equal to lambda.
  }
  \item{method}{
    parameter passed to function \code{\link{optim}}.
}
\item{lower}{
      parameter passed to function \code{\link{optim}}.
}
  \item{upper}{
    parameter passed to function \code{\link{optim}}.
}
  \item{\dots}{
    further parameters passed to function \code{\link{optim}}.
  }
}

\value{
   an object of class 'htest'.
}

\references{
  C. Agostinelli, A. Marazzi and V.J. Yohai (2015).
  Robust estimates of the generalized loggamma distribution.
  Technometrics, Volume 56, Issue 1, 2014.
  doi:10.1080/00401706.2013.818578
}

\author{
  C. Agostinelli, A. Marazzi, V.J. Yohai and A. Randriamiharisoa
}

\seealso{
\code{\link{loggammarob.test}} for robust (weighted) t-test on one or
more parameters of the generalized loggamma model.
}

\examples{
  set.seed(1234)
  x <- sort(log(rgamma(n=30, shape=2, scale=2)))
  res <- loggammarob(x)
  loggammarob.wilks(res)
}

\keyword{robust}
\keyword{models}
