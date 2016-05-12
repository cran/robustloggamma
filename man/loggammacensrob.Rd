\name{loggammacensrob}
\alias{loggammacensrob}

\title{
  Robust estimation of the three parameters extended Log Gamma model in
  presence of censored observations.
}
\description{
  Five different type of robust procedures are provided for the
  estimation of the parameters in the three parameters extended Log
  Gamma model in presence of censored observations.
}
\usage{
loggammacensrob(x, delta, start=NULL, weights=rep(1, length(x)),
   method=c("oneTML", "oneWL", "TWQTau", "TQTau", "ML"), control, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric. A vector with the dataset.}
  \item{delta}{numeric, 0 or 1. 0 indicated censored observations.}  
  \item{start}{\code{NULL} or numeric. A vector of length 3 to be used
  when \code{method} is \code{WL}, \code{oneWL} and \code{ML}, otherwise
  starting values are obtained from \code{WQTau} in the first two
  methods and \code{QTau} in the last one.}
  \item{weights}{numeric. A vector of weights used in the method \code{TQTau}.}
  \item{method}{character. The method used, see Details below. Default is
  \code{oneTML} a one step Trunceted Maximum Likelihood estimates
  starting 
  from \code{TQTau}.}
  \item{control}{list. An object from function
    \code{\link{loggammarob.control}}.}
  \item{\dots}{further arguments can be passed directly to the function
    instead of using the \code{control} argument.}
}
% \details{
% TO BE DONE!
% }
\value{
  An object of class 'loggammacensrob'. A list that includes the
  following components:
\item{mu}{location parameter.}
\item{sigma}{scale parameter.}
\item{lambda}{shape parameter.}
\item{eta}{estimate of E(exp(x)) parameter.}
\item{weights}{the final weights.}
\item{iterations}{number of iterations.}
\item{error}{if not \code{NULL} or 0, an error occured.}
\item{data}{the original dataset.}
\item{delta}{the original dataset.}
\item{method}{the method used.}
}

\references{
  C. Agostinelli, I. Locatelli, A. Marazzi and V.J. Yohai (2015)
  Robust estimators of accelerated failure time
  regression with generalized log-gamma errors. Submitted.
}

\author{
  C. Agostinelli, A. Marazzi and V.J. Yohai
}

\seealso{
  \code{\link{LogGammaDist}} and \code{\link{loggammarob}} for the case
  without censored observations.
}
\examples{
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  cens <- rloggamma(n=80, mu=1.5, lambda=1)
  delta <- as.numeric(x <= cens)
  x[delta==0] <- cens[delta==0]  
  res <- loggammacensrob(x, delta, method="oneTML",
    control=loggammarob.control(lower=0, upper=2, n=30))
  print(res)
  summary(res)
}
\keyword{robust}
\keyword{models}
