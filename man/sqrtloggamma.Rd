\name{sqrtloggamma}
\docType{data}
\alias{sqrtloggamma}
\title{
A reparametrization for the sigma parameter in loggamma model. 
}
\description{
A reparametrization for the sigma parameter in loggamma model. This is
used in the function \code{loggammarob.control}.
}
\usage{sqrtloggamma}

\format{
  The format is a list of 3 elements
  
  $ gam   :function (sigma)
  
  $ gaminv:function (gam)
  
  $ delta :function (sigma)  
}

\details{
The object must be a length of dimension 3. Each component is a
  function. The first component \code{gam} performs the transformation
  for the parameter \code{sigma}; \code{gaminv} is the inverse function
  and \code{delta} is the derivative of the \code{gam} function wrt
  \code{sigma} to the -1 power.  
}

\examples{
  str(sqrtloggamma)
  set.seed(1234)
  x <- sort(rloggamma(n=50, lambda=1))
  res <- loggammarob(x, control=loggammarob.control(lower=0, upper=2,
  n=30, reparam=sqrtloggamma))
  print(res)
}

\keyword{datasets}
