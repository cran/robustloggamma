# TML_cens_loggamma_functions

# Distributions, Scores, and derivatives
# ======================================

tdloggamma <- function(x,lambda,tl,tu,zero) {
  if (missing(zero))
    zero <- 1e-4
  cond <- (tl < x & x < tu)*1
  dens <- dloggamma(x=x,lambda=lambda,zero=zero)
  prob <- ploggamma(q=tu,lambda=lambda,zero=zero)-ploggamma(q=tl,lambda=lambda,zero=zero)
  dens <- dens/prob
  ans <- dens*cond
  return(ans)
}

tploggamma <- function(x,lambda,tl,tu,zero) {
  if (missing(zero))
    zero <- 1e-4
  if (x  <= tl) cdf <- 0
  if (tu <= x)  cdf <- 1
  if (tl < x & x < tu) {
    cdf  <- ploggamma(q=x,lambda=lambda,zero=zero)-ploggamma(q=tl,lambda=lambda,zero=zero)
    prob <- ploggamma(q=tu,lambda=lambda,zero=zero)-ploggamma(q=tl,lambda=lambda,zero=zero)
    cdf  <- cdf/prob
  } 
  return(cdf)
}

tqloggamma <- function(p,lambda,tl,tu,zero) {
  if (missing(zero))
    zero <- 1e-4
  pu <- ploggamma(q=tu,lambda=lambda,zero=zero)
  pl <- ploggamma(q=tl,lambda=lambda,zero=zero)
  qn <- qloggamma(p=pl+(pu-pl)*p,lambda=lambda,zero=zero)
  return(qn)
}

l.score.loggamma <- function(u,lambda,zero) {
  if (missing(zero))
    zero <- 1e-4   
# negative score for lambda in loggamma model
  if (abs(lambda) > zero) {
    lamu   <- u*lambda
    elamu  <- exp(lamu)
    lamm1  <- lambda^(-1)
    lamm2  <- lambda^(-2)
    llamm2 <- log(lamm2)
    lamm3  <- lambda^(-3)
    arg1d  <- lamm1-2*lamm3*llamm2-2*lamm3-lamm2*u+elamu*(2*lamm3-lamm2*u)
    digl   <- digamma(lamm2)
    res    <- -(arg1d+2*lamm3*digl)
  } else {
    res <- (u^3)/6
  }
  return(res)
}

####Exp.response <- function(mu, sigma, lambda, eps=0.0001, npoints=100000) Questa e' l'api nel pacchetto!
## Exp.response <- function(lambda,mu,sigma,zero=zro){
## nm <- 100000
## if (is.na(lambda)) return(NA)
## if (lambda > zero) {
##   alpha <- 1/lambda^2
##   ca    <- lambda/sigma
##   lMu   <- mu - log(alpha)/ca + lgamma(alpha+1/ca)-lgamma(alpha)
##   Mu    <- exp(lMu) }
## if (abs(lambda) <= zero) Mu <- exp(mu + sigma^2/2)
## if (lambda < -zero) {
##   ql  <- qloggamma(ppoints(nm),lambda)
##   Mu  <- mean(exp(mu+sigma*ql)) }
## Mu}

fp.over.f.loggamma <- function(u,lambda,zero) {
# score csi_lambda
  if (missing(zero))
    zero <- 1e-4
  if (abs(lambda) > zero) {
    res <- (1-exp(lambda*u))/lambda
  } else {
    res <- -u
  }
  return(res)
}

dersig.fp.over.f.loggamma <- function(u,lambda,sigma,zero) {
  if (missing(zero))
    zero <- 1e-4
# score derivative
  if (abs(lambda) > zero) {
    res <- u*exp(lambda*u)/sigma
  } else {
    res <- u/sigma
  }
  return(res)
}

ld.score.loggamma <- function(u,lambda,zero) {
  if (missing(zero))
    zero <- 1e-4
# derivative with respect of lambda of the score for lambda
  lamu  <- u*lambda
  elamu <- exp(lamu)
  if (abs(lambda) > zero) {
    lamm2  <-  lambda^(-2)
    llamm2 <-  log(lamm2)
    lamm3  <-  lambda^(-3)
    lamm4  <-  lambda^(-4)
    lamm6  <-  lambda^(-6)
    arg1dd <- -lamm2+6*lamm4*llamm2+10*lamm4+2*lamm3*u-6*lamm4*elamu+4*lamm3*elamu*u-lamm2*elamu*u^2
    trigl  <-  trigamma(lamm2)
    digl   <-  digamma(lamm2)
    res    <- -arg1dd+6*lamm4*digl+4*lamm6*trigl
  } else {
    res <-(u^4)/12
  }
  return(res)
}

# Derivatives of cdf F_lambda(u) and density f_lambda(u)
# ======================================================

phiLG <- function(u,lambda){
  csiLG(u,lambda)*u
}

phiLG.dot <- function(u,lambda){
  csiLG.dot(u,lambda)*u
}

phiLG.prime <- function(u,lambda) {
  csiLG.prime(u,lambda)*u +csiLG(u,lambda)
}
                                   
# Weight functions
# ================

weightsHD <- function(r,tl,tu){ # rectangular weight function
nr   <- length(r); tmp  <- rep(0,nr)
ind  <- (1:nr)[r > tl & r < tu] 
tmp[ind] <- 1; tmp}

rho.weightsHD <- function(r,rhotu){ # hard weights on data scale
nr   <- length(r); tmp  <- rep(0,nr)
ind  <- (1:nr)[r < rhotu]
tmp[ind] <- 1; tmp}

rhsHD <- function(lambda,zl,zu,zero,npoints=100000) {
  if (missing(zero))
    zero <- 1e-4
  if (abs(lambda) > zero)
    ui <- qloggamma(ppoints(npoints),lambda=lambda)
  else
    ui <- qnorm(ppoints(npoints),mean=0,sd=1)
  wi <- weightsHD(ui,zl,zu)
  hh <- mean(wi*fp.over.f.loggamma(ui,lambda) )
  kk <- mean(wi*fp.over.f.loggamma(ui,lambda)*ui)
  gg <- mean(wi*l.score.loggamma(ui,lambda))
  res <- list(hh=hh,kk=kk,gg=gg)
  return(res)
}

