\name{alcoa}
\alias{alcoa}
\docType{data}
\title{
Quality assurance measurement on aluminium.
}
\description{
Three samples (with labels A,B and C) from measurement quality assurance
(QA) data base of ALCOA aluminium refineries in Western Australia.
}
\usage{data(alcoa)}
\format{
  A data frame with 668 observations on the following 3 variables.
  \describe{
    \item{\code{dist}}{a numeric vector}
    \item{\code{ratio}}{a numeric vector}
    \item{\code{label}}{a factor with levels \code{A} \code{B} \code{C}}
  }
}

\details{
Under ALCOA's QA program, several thousand bauxite ore samples are
routinely submitted to Fourier transform infrared spectroscopy
(Eyer and Riley, 1999). Part of the quality assurance is the need to
automatically highlight unusual spectra and this is obtained with the
help of special statistical diagnostics - called representation
indicators - derived from the Fourier transform.
}

\source{
  Alcoa World Alumina

  Clarke B.R., McKinnon P.L., Riley G. (2012). A fast robust method
  for fitting gamma distributions. Statistical Papers, 53, 4, 1001-1014.
}

\references{
  Agostinelli C., Marazzi A, Yohai V.J. (2012). Robust estimators
  of the generalized loggamma distribution. Submitted.

  Clarke B.R., McKinnon P.L., Riley G. (2012). A fast robust method
  for fitting gamma distributions. Statistical Papers, 53, 4, 1001-1014.

  Eyer S., Riley G. (1999). Measurement quality assurance in a
  production system for bauxite analysis by FTIR. North American Chapter
  of the International Chemometrics Society, Newsl No. 19.
}

\examples{
data(alcoa)
par(mfcol=c(1,2))
boxplot(I(log(alcoa$ratio))~alcoa$label)
boxplot(I(log(alcoa$dist))~alcoa$label)
}

\keyword{datasets}
