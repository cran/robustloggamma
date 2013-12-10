\name{loggammarob}
\alias{loggammarob}
\alias{print.loggammarob}

\title{
  Robust estimation of the three parameters extended Log Gamma model.
}
\description{
Four different type of robust procedures are provided for the estimation
  of the parameters in the three parameters extended Log Gamma model  
}
\usage{
loggammarob(x, start=NULL, weights = rep(1, length(x)),
  method=c("oneWL", "WQTau", "WL", "QTau", "ML"), control, ...)
\method{print}{loggammarob}(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric. A vector with the dataset.}
  \item{start}{\code{NULL} or numeric. A vector of length 3 to be used
  when \code{method} is \code{WL}, \code{oneWL} and \code{ML}, otherwise
  starting values are obtained from \code{WQTau} in the first two
  methods and \code{QTau} in the last one.}
  \item{weights}{numeric. A vector of weights used in the method \code{QTau}.}
  \item{method}{character. The method used, see Details below. Default is
  \code{oneWL} a one step weighted likelihood estimates starting from \code{WQTau}.}
  \item{control}{list. An object from function
    \code{\link{loggammarob.control}}.}
  \item{digits}{minimal number of _significant_ digits, see \code{\link{print.default}}.}
  \item{\dots}{further arguments can be passed directly to the function
    instead of using the \code{control} argument.}
}
% \details{
% TO BE DONE!
% }
\value{
An object of class 'loggammarob'. A list that includes the following components:
\item{mu}{location parameter.}
\item{sigma}{scale parameter.}
\item{lambda}{shape parameter.}
\item{eta}{estimate of E(exp(x)) parameter.}
\item{weights}{the final weights.}
\item{iterations}{number of iterations.}
\item{error}{used only in method \code{oneWL}. If \code{1} then the
  Jacobian matrix is not invertible.}
\item{data}{the original dataset.}
\item{method}{the method used.}
}

\references{
  C. Agostinelli, A. Marazzi and V.J. Yohai (2013) Robust estimates of
  the generalized loggamma distribution, Technometrics, Accepted.
}

\author{
  C. Agostinelli, A. Marazzi, V.J. Yohai and A. Randriamiharisoa
}

\seealso{
  \code{\link{LogGammaDist}}.
}
\examples{
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  res <- loggammarob(x, control=loggammarob.control(lower=0, upper=2, n=30))
  print(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}
\keyword{models}
