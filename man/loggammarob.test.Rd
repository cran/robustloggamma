\name{loggammarob.test}
\alias{loggammarob.test}

\title{
Robust inference for the generalized loggamma model. 
}
\description{
Robust Tests and confidence intervals for the parameters of the generalized
loggamma model.
}
\usage{
loggammarob.test(x, mu = NULL, sigma = NULL, lambda = NULL,
eta = NULL, type = "Wald", conf.level = 0.95, prob = 0.00001)
}
\arguments{
  \item{x}{list. An object of class 'loggammarob'}
  \item{mu}{numeric. Null value for the location parameter}
  \item{sigma}{numeric. Null value for the scale parameter}
  \item{lambda}{numeric. Null value for the shape parameter}
  \item{eta}{numeric. Null value for E(exp(X)) parameter}
  \item{type}{character. Type of inference, for now only Wald test is performed}
  \item{conf.level}{numeric. Level of the confidence interval.}
  \item{prob}{quantile order of the loggamma to be used in the
    numerical calculation of the expected Fisher Information.}
}
\details{
If no null values are provided, the function report the test for
\code{mu=0}, \code{sigma=1} and \code{lambda=0}.
}

\value{
An object of class 'htest' containing the following components:
\item{statistic:}{the value of the t-statistic.}
\item{parameter:}{the degrees of freedom for the statistic.}
\item{p.value:}{the p-value for the test.}
\item{conf.int:}{ a confidence interval for the parameter appropriate to
  the specified alternative hypothesis.}
\item{estimate:}{the estimate of the parameter(s).}
\item{null.value:}{the specified hypothesized value of the parameter.}
\item{alternative:}{a character string describing the alternative hypothesis.}
\item{method:}{a character string indicating what type of test was performed.}
\item{data.name:}{a character string giving the name(s) of the data.}
}

\references{
  C. Agostinelli, A. Marazzi and V.J. Yohai (2015) Robust estimates of
  the generalized loggamma distribution, Technometrics, Volume 56, Issue
  1, 2014. DOI: 10.1080/00401706.2013.818578

  Agostinelli C., Marazzi A., Yohai V.J., Randriamiharisoa A. (2016).
  Robust Estimation of the Generalized Loggamma Model:
  The R Package robustloggamma.
  Journal of Statistical Software, 70(7), 1-21.
  doi:10.18637/jss.v070.i07
}

\author{
  A. Marazzi, C. Agostinelli, V.J. Yohai and A. Randriamiharisoa
}

\seealso{
  \code{\link{loggammarob}}
}
\examples{
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  res <- loggammarob(x, control=loggammarob.control(lower=0, upper=2, n=30))
  loggammarob.test(res, mu=0) #only location
  loggammarob.test(res, mu=0, sigma=1) #location and scale
  loggammarob.test(res, eta=1) #E(exp(X))
}

\keyword{robust}
\keyword{models}
