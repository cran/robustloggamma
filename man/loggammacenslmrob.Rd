\name{loggammacenslmrob}
\alias{loggammacenslmrob}
\title{
  Robust estimation of an Accelerated Failure Time model with
  extended Log Gamma errors in presence of censored observations.
}
\description{
  Three different type of robust procedures are provided for the
  estimation of the parameters in an Accelerated Failure Time
  model with extended Log Gamma errors in presence of censored
  observations. Maximum Likelihood is also provided.
}
\usage{
loggammacenslmrob(formula, delta, data, subset, weights,
  na.action, method = c("oneTML", "oneWL", "TWQTau", "TQTau",
  "ML"), model = TRUE, x = !control$compute.rd, y = FALSE,
  singular.ok = TRUE, contrasts = NULL, offset = NULL,
  control = NULL, init = NULL, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.  See
    \code{\link{lm}} and \code{\link{formula}} for more details.}
  \item{delta}{numeric, 0 or 1. 0 indicated censored observations.}  
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing the
    variables in the model. If not found in \code{data}, the variables
    are taken from \code{environment(formula)}, typically the
    environment from which \code{loggammacenslmrob} is called.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{weights}{an optional vector of weights to be used in the
    fitting process (in addition to the robustness weights computed
    in the fitting process).}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  The default is set by the
    \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}. Another possible value is
    \code{NULL}, no action. Value \code{\link{na.exclude}} can be
    useful.}
  \item{method}{string specifying the estimator-chain. Default is TML
    one step.}
  \item{model, x, y}{logicals. If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response) are returned.}
  \item{singular.ok}{logical. If \code{FALSE} (the default in S but not
    in R) a singular fit is an error.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    An \code{\link{offset}} term can be included in the formula instead
    or as well, and if both are specified their sum is used.}
  \item{control}{a \code{\link{list}} specifying control parameters; use the
    function \code{\link{loggammarob.control}} and see its help page.}
  \item{init}{an optional argument to specify or supply the initial
    estimate. See \emph{Details}.}
  \item{\dots}{additional arguments can be used to specify control
    parameters directly instead of (but not in addition to!) via
    \code{control}.}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
 An object of class \code{lmrob}; a list including the following
  components:
  \item{coefficients}{The estimate of the coefficient vector for the
    regression part. First element is the intercept and it would be
    equals to the parameter \code{mu}.}
  \item{mu}{The estimate of the intercept parameter.}
  \item{sigma}{The estimate of the scale parameter.}
  \item{lambda}{The estimate of the shape parameter.}
  \item{fitted.values}{Fitted values associated with the estimator.}  
  \item{residuals}{Residuals associated with the estimator.}
  \item{cut.lower, cut.upper}{Cut points for the method based on TML
    (missing for the other methods).}
  \item{iter}{number of iterations.}
  \item{weights}{the specified weights (missing if none were used).}
  \item{errors}{errors messages.}  
  \item{n.ret}{number of non zero robust weights for the method
    based on TML.}
  \item{control}{\code{control} argument.}  
  \item{converged}{\code{TRUE} if the procedure converged.}
  \item{method}{method used during the fit.}
  \item{rank}{the numeric rank of the fitted linear model.}
  \item{rweights}{the \dQuote{robustness weights}.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{degree.freedom}{the same as df.residual}
  \item{delta}{as in input.}
  \item{df}{a vector with 3 components: (number of linearly independent
    regressors, df.residual, number of regressors).}
  \item{xlevels}{(only where relevant) a record of the levels of the
    factors used in fitting.}
  \item{call}{the matched call.}
  \item{terms}{the \code{terms} object used.}
  %assign
  \item{model}{if requested (the default), the model frame used.}
  \item{x}{if requested, the model matrix used.}
  \item{y}{if requested, the response used.}
  \item{scale}{square root of the \code{sigma} parameter.}
  \item{na.action}{(where relevant) information returned by
     \code{\link{model.frame}} on the special handling of \code{NA}s.}
  \item{offset}{the offset used (missing if none were used).}
  \item{contrasts}{(only where relevant) the contrasts used.}
  
  In addition, non-null fits will have \code{qr} relating to the linear
  fit, for use by extractor functions such as \code{summary}.
}
\references{
  C. Agostinelli, I. Locatelli, A. Marazzi and V.J. Yohai (2015)
  Robust estimators of accelerated failure time
  regression with generalized log-gamma errors. Submitted.
}
\author{
  C. Agostinelli, A. Marazzi and V.J. Yohai
}
%% \note{
%% %%  ~~further notes~~
%% }
\seealso{
  \code{\link{loggammacensrob}} for the case of censored observations
  without covariates. 
}
\examples{
  \dontrun{
    n <- 50
    p <- 2
    set.seed(1234)
    X <- matrix(rnorm(p*n, sd=2), ncol=p)
    mu <- 2
    beta <- rep(2,p)
    sigma <- 2
    lambda <- 1
    linear <- mu + X%*%beta
    y <- rloggamma(n=n, mu=linear, sigma=sigma, lambda=lambda)
    cens <- rloggamma(n=n, mu=linear+3, sigma=sigma, lambda=lambda)
    delta <- as.numeric(y <= cens)
    y[delta==0] <- cens[delta==0]
    x <- data.frame(y=as.vector(y), x1=X[,1], x2=X[,2])
    res <- loggammacenslmrob(y~x1+x2, delta=delta, data=x)
    summary(res)
  }
}
\keyword{robust}
\keyword{models}

