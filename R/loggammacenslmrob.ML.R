loggammacenslmrob.ML <- function(x, y, delta, init=NULL, w, control, ...) {
  ## par[1] location
  ## par[p+1] sigma
  ## par[p+2] lambda
  ## resto beta
  x <- as.matrix(x)
  p <- NCOL(x)
  if (is.null(init)) {
    keep <- w >= control$minw
    MM <- loggammacenslmrob.MM(x=cbind(1,x[keep,]),y=y[keep],delta=delta[keep],control=loggammacenslmrob.MM.control())
    init <- MM$coefficients
    init2 <- loggammacensrob.TQTau(x=MM$residuals, delta=delta, w=w, control=control)
    init <- c(init2$mu, init[-1], init2$sigma, init2$lambda)
  }
  lower <- c(rep(-Inf,p+1), 0.001, control$lower)
  upper <- c(rep(Inf,p+2), control$upper)
  pos <- w >= control$minw
  y <- y[pos]
  delta <- delta[pos]
  x <- x[pos,,drop=FALSE]
  w <- w[pos]
  res <- optim(par=init, fn=loggamma.loglikreg, y=y, delta=delta, X=x, w=w, lower=lower, upper=upper, method="L-BFGS-B", ...)
  ans <- list()
  ans$mu <- res$par[1]
  ans$coefficients <- res$par[1:(p+1)]
  ans$sigma <- res$par[p+2]
  ans$lambda <- res$par[p+3]
  ans$fitted.values <- drop(x%*%ans$coefficients[-1] + ans$mu)
  ans$residuals <- y - ans$fitted.values
  ans$iter <- res$counts[1]
  ans$weights <- w
  ans$control <- control
  ans$converged <- res$convergence == 0
  return(ans)
}

loggamma.loglikreg <- function(par, y, delta, X, w, xmax=1e100) {
##  cat(par, "\n")  
  ## par[1] location
  ## par[p+2] sigma
  ## par[p+3] lambda
  ## resto beta
  p <- NCOL(X) ## number of regressors without intercept
  mu <- par[1]+X%*%par[2:(p+1)]
  sigma <- par[p+2]
  lambda <- par[p+3]
  uno <- delta==1
  zero <- !uno
  fi <- w[uno]*dloggamma(y[uno], mu=mu[uno], sigma=sigma, lambda=lambda, log=TRUE)
  Si <- w[zero]*ploggamma(y[zero], mu=mu[zero], sigma=sigma, lambda=lambda, log.p=TRUE, lower.tail=FALSE)
  res <- -(sum(fi)+sum(Si))
  if (res > xmax | is.na(res))
    res <- xmax
    ##res <- .Machine$double.xmax/2
##  cat(res, "\n")
  return(res)
}

VCOV.WMLCensloggammareg <- function(X, y, delta, mu, sigma, lambda, beta, w, xmax=1e100, minw=0.05) {
  par <- c(mu,beta,sigma,lambda)
  pos <- w >= minw
  y <- y[pos]
  delta <- delta[pos]
  X <- X[pos,,drop=FALSE]
  w <- w[pos]
  p <- length(par)
  n <- length(y)
  Var <- matrix(NA, nrow=p, ncol=n)
  H <- matrix(0, nrow=p, ncol=p)
  for (i in 1:n) {
    Var[,i] <- grad(func=function(x) loggamma.loglikreg(par=x, y=y[i], delta=delta[i], X=X[i,,drop=FALSE], w=w[i], xmax=xmax), x=par)
    H <- H + hessian(func=function(x) loggamma.loglikreg(par=x, y=y[i], delta=delta[i], X=X[i,,drop=FALSE], w=w[i], xmax=xmax),  x=par)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
  return(vcov)
}

