loggamma.loglik <- function(par, y, delta, w, xmax=1e100) {
  ## par[1] location
  ## par[2] sigma
  ## par[3] lambda
  mu <- par[1]
  sigma <- par[2]
  lambda <- par[3]
  uno <- delta==1
  zero <- !uno
  fi <- w[uno]*dloggamma(y[uno], mu=mu, sigma=sigma, lambda=lambda, log=TRUE)
  Si <- w[zero]*ploggamma(y[zero], mu=mu, sigma=sigma, lambda=lambda, log.p=TRUE, lower.tail=FALSE)
  res <- -(sum(fi)+sum(Si))
  if (res > xmax | is.na(res))
    res <- xmax
  return(res)
}

loggammacensrob.ML <- function(x, delta, start=NULL, w=rep(1, length(x)), control, ...) {
  if (is.null(start)) {
    TQTau <- loggammacensrob.TQTau(x=x, delta=delta, w=w, control=control)
    start <- c(TQTau$mu, TQTau$sigma, TQTau$lambda)
  } else {
    if (!is.numeric(start))
      stop("'start' must be a numeric vector")
    if (length(start)!=3)
      stop("'start' must be a vector of length 3: mu, sigma2, lambda")    
    TQTau <- NULL
  }
  ## start[1] location
  ## start[2] sigma
  ## start[3] lambda
  w[w > 1-control$minw] <- 1  
  pos <- w >= control$minw
  y <- y[pos]
  delta <- delta[pos]
  w <- w[pos]
  control$lower <- c(-Inf, 0.001, control$lower)
  control$upper <- c(Inf, Inf, control$upper)
  res <- optim(par=start, fn=loggamma.loglik, y=x, delta=delta, w=w, lower=control$lower, upper=control$upper, method="L-BFGS-B", ...)
  result <- list(mu=res$par[1], sigma=res$par[2], lambda=res$par[3], weights=w, iterations=res$counts[1], error=res$convergence)
  return(result)
}


WMLCensloggammaWeights.control <- function(subdivisions=1000,tol=1e-4,bw=0.3,raf="NED",tau=0.5) {
  control <- list(subdivisions=subdivisions,tol=tol,bw=bw,raf=raf,tau=tau)
  return(control)
}

WMLCensloggammaWeights <- function(y, delta, mu0, sig0, lam0, control) {
  qemp <- SemiQuant(y=y, delta=delta, mu0, sig0, lam0, nmod=control$subdivisions, tol=control$tol) # 1000 quantili della semi-empirica basata sul WQtau
  w <- weightsCensloggamma(y=y, delta=delta, qemp=qemp, mu=mu0, sigma=sig0, lambda=lam0, bw=control$bw, raf=control$raf, tau=control$tau, nmod=control$subdivisions)
}

VCOV.WMLCensloggamma <- function(par, y, delta, w, xmax=1e100, minw=0.005) {
  ## par[1] location
  ## par[2] sigma
  ## par[3] lambda
  w[w > 1-minw] <- 1    
  pos <- w >= minw
  y <- y[pos]
  delta <- delta[pos]
  w <- w[pos] 
  p <- length(par)
  n <- length(y)
  Var <- matrix(NA, nrow=p, ncol=n)
  H <- matrix(0, nrow=p, ncol=p)
  for (i in 1:n) {
    Var[,i] <- grad(func=function(x) loggamma.loglik(par=x, y=y[i], delta=delta[i], w=w[i], xmax=xmax), x=par)
    H <- H + hessian(func=function(x) loggamma.loglik(par=x, y=y[i], delta=delta[i], w=w[i], xmax=xmax),  x=par)
  }
  Var <- cov(t(Var))
  H <- solve(H)
  vcov <- H%*%Var%*%t(H)*n
  return(vcov)
}
