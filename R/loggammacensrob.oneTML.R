loggammacensrob.oneTML <- function(x, delta, start=NULL, w=rep(1, length(x)), control) {
  if (is.null(start)) {
    TQTau <- loggammacensrob.TQTau(x=x, delta=delta, w=w, control=control)
    mu0 <- TQTau$mu
    sig0 <- TQTau$sigma
    lam0 <- TQTau$lambda
  } else {
    if (!is.numeric(start))
      stop("'start' must be a numeric vector")
    if (length(start)!=3)
      stop("'start' must be a vector of length 3: mu, sigma2, lambda")    
    TQTau <- NULL
    mu0 <- start[1]
    sig0 <- start[2]
    lam0 <- start[3]
  }
  rs0 <- (x-mu0)/sig0
  max.iter <- min(control$max.iter,control$iter)
  dif <- 1+control$refine.tol
  iter <- 0
  while (dif > control$refine.tol & iter < max.iter & lam0 > control$lower & lam0 < control$upper) {
    iter <- iter + 1
    if (is.null(control$xmax)) {
      control$xmax <- 0
      while (ploggamma(control$xmax,lambda=lam0) < 1) {
        control$xmax <- control$xmax+0.01
      }
    }
    rmax <- max(control$xmax,max(rs0))+0.1
    cu <- tutl.loggamma(control$pcut,lam0)
    grid <- seq(from=cu,to=rmax,length=control$nTML)
    zap <- adapt.CutoffLG(rs0,delta,lam0,grid)
    cl <- zap$xl
    cu <- zap$xu+0.001
# TMLone 
    TMo <- TMLone.loggamma(IX=matrix(1, nrow=length(x)), y=x, delta=delta, Beta.t=mu0, sigma.t=sig0, lambda.t=lam0, cl=cl, cu=cu, step=control$step, reparam=FALSE)
    if (as.numeric(TMo$sigma) < 0 ) {
      TMo <- TMLone.loggamma(IX=matrix(1, nrow=length(x)), y=x, delta=delta, Beta.t=mu0, sigma.t=sig0, lambda.t=lam0, cl=cl, cu=cu, step=control$step, reparam=TRUE)
    }
    dif <- max(abs(c(mu0-drop(TMo$Beta[1]), sig0-TMo$sigma, lam0-TMo$lambda)))
    mu0 <- drop(TMo$Beta[1])
    sig0 <- TMo$sigma
    lam0 <- TMo$lambda
  }

  result <- list(mu=mu0, sigma=sig0, lambda=lam0, weights=TMo$weights, iterations=iter, error=TMo$err, TQTau=TQTau, cut.lower=cl, cut.upper=cu)
  return(result)
}

VCOV.TMLoneCensloggamma <- function(y,delta,w,mu,sigma,lambda,cl,cu) {  
  p <- 1
  n <- length(y)
  rs0 <- as.vector((y - mu))/sigma
  pos <- w > 0.001
  IX <- matrix(rep(1,n), ncol=1)
  IX <- IX[pos, , drop = FALSE]
  y <- y[pos]
  delta <- delta[pos]
  w <- w[pos]
  rs0 <- rs0[pos]
  kk <- -cl*dloggamma(cl,mu=0,sigma=1,lambda=lambda) + cu*dloggamma(cu,mu=0,sigma=1,lambda=lambda)-ploggamma(cu,mu=0,sigma=1,lambda=lambda)+ploggamma(cl,mu=0,sigma=1,lambda=lambda)
  gg <- d.FL.1(cl,lambda)-d.FL.1(cu,lambda)
  JAC <- JAC.TML.loggamma(IX, y, delta, mu, mu, sigma, sigma, lambda, rs0, w, cl, cu, kk, gg)
  Q <- QMatrix.LG(rep(0,p),1,lambda,rs0,delta,w,IX,cl,cu,kk,gg) 
  JIC <- solve(JAC)
  Cov <- t(JIC) %*% Q %*% JIC /n
  return(Cov)
}

TMLCensloggammaWeights.control <- function(pcut=0.997,nTML=2000,xmax=NULL,step=1) {
  control <- list(pcut=pcut,nTML=nTML,xmax=xmax,step=step)
  return(control)
}

TMLCensloggammaWeights <- function(y,delta,mu0,bet0,sig0,lam0,control) {
  rs0 <- (y-mu0)/sig0
  if (is.null(control$xmax)) {
    control$xmax <- 0
    while (ploggamma(control$xmax,lambda=lam0) < 1) {
      control$xmax <- control$xmax+0.01
    }
  }
  rmax <- max(control$xmax,max(rs0))+0.1
  cu <- tutl.loggamma(control$pcut,lam0)
  grid <- seq(from=cu,to=rmax,length=control$nTML)
  zap <- adapt.CutoffLG(rs0,delta,lam0,grid)
  cl <- zap$xl
  cu <- zap$xu+0.001
  wgt <- weightsHD(rs0, cl, cu)
  return(wgt)
}
