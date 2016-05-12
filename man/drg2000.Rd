\name{drg2000}
\alias{drg2000}
\docType{data}
\title{
Length of stay, cost and AP-DRG for several hospital in Switzerland in 2000.
}
\description{
 The data refer to 70323 stays that were observed in year 2000 in a group of Swiss hospitals within a pilot study aimed at the implementation of a diagnosis-related grouping (DRG) system. DRG systems are used in modern hospital management to classify each individual stay into a group according to the patient characteristics.
}
\usage{data(drg2000)}
\format{
  A data frame with 70323 observations on the following 4 variables.
  \describe{
    \item{\code{LOS}}{numeric. Lenght of the hospitalization (LOS).}
    \item{\code{Cost}}{numeric. Cost of the hospitalization in
      Switzerland franc.}
    \item{\code{MDC}}{numeric. Major Diagnostic Category.}
    \item{\code{APDRG}}{numeric. All Patient Diagnosis Related Group.}
  }
}
%% \source{

%% }
\references{
  Agostinelli C., Marazzi A., Yohai V.J., Randriamiharisoa A. (2016).
  Robust Estimation of the Generalized Loggamma Model:
  The R Package robustloggamma.
  Journal of Statistical Software, 70(7), 1-21.
  doi:10.18637/jss.v070.i07
}

\examples{
data(drg2000)
str(drg2000)
}

\keyword{datasets}
