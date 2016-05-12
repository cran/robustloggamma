loggammacenslmrob.oneWL <- function(x,y,delta,init=NULL,w,control) {
  control$d <- 100
  x <- as.matrix(x)
  p <- NCOL(x)
  if (is.null(init)) {
    keep <- w >= control$minw
    MM <- loggammacenslmrob.MM(x=cbind(1,x[keep,]),y=y[keep],delta=delta[keep],control=loggammacenslmrob.MM.control())
    init <- MM$coefficients
    residuals <- drop(y - x%*%init[2:(p+1)])    
    init2 <- loggammacensrob.TQTau(x=residuals, delta=delta, w=w, control=control)
    mu0 <- init2$mu
    bet0 <- init[-1]
    sig0 <- init2$sigma
    lam0 <- init2$lambda
  } else {
    mu0 <- init[1]
    bet0 <- init[2:(p+1)]
    sig0 <- init[p+2]
    lam0 <- init[p+3]
  }
  
  v <- y - x%*%bet0
  qemp <- SemiQuant(y=v,delta=delta,mu0,sig0,lam0,nmod=control$subdivisions, tol=control$refine.tol)
  w <- weightsCensloggamma(y=v,delta=delta,qemp=qemp,
                            mu=mu0,sigma=sig0,lambda=lam0,
                            bw=control$bw,raf=control$raf,
                            tau=control$tau,nmod=control$subdivisions)
  w[w > 1-control$minw] <- 1  
  pos <- w >= control$minw
  ## y <- y[pos]
  ## delta <- delta[pos]
  ## w <- w[pos]
  ## x <- x[pos,,drop=FALSE]
  control$expJ <- FALSE
  if (control$expJ) {
    zWL <- WLoneCensloggammareg.expJ(y=y[pos],delta=delta[pos],
            x=x[pos,,drop=FALSE],w=w[pos],
            mu0=mu0,bet0=bet0,sig0=sig0,
            lam0=lam0,qemp=qemp,step=control$step,d=control$d,
            reparam=FALSE)
    if (zWL$sigma < 0) {
      zWL <- WLoneCensloggammareg.expJ(y=y[pos],delta=delta[pos],
            x=x[pos,,drop=FALSE],w=w[pos],
            mu0=mu0,bet0=bet0,sig0=sig0,
            lam0=lam0,qemp=qemp,step=control$step,d=control$d,reparam=TRUE)
      cat("Repar","\n")
    }
  } else {
    zWL <- WLoneCensloggammareg.empJ(y=y[pos],delta=delta[pos],
             x=x[pos,,drop=FALSE],w=w[pos],
             mu0=mu0,bet0=bet0,sig0=sig0,
             lam0=lam0,step=control$step,d=control$d,
             reparam=FALSE) # WLone
    if (zWL$sigma < 0) {
      zWL <- WLoneCensloggammareg.empJ(y=y[pos],delta=delta[pos],
             x=x[pos,,drop=FALSE],w=w[pos],
             mu0=mu0,bet0=bet0,sig0=sig0,
             lam0=lam0,step=control$step,d=control$d,reparam=TRUE)
      cat("Repar","\n")
    }
  }
  res <- list()
  res$coefficients <- c(zWL$mu,zWL$beta)
  res$mu <- zWL$mu
  res$sigma <- zWL$sigma
  res$lambda <- zWL$lambda
  res$fitted.values <- drop(x%*%res$coefficients[-1] + res$mu)
  res$residuals <- y - res$fitted.values
  res$iter <- control$step
  res$weights <- w
  res$control <- control
  res$converged <- TRUE
  return(res)
}

VCOV.WLoneCensloggammareg.empJ <- function(y, delta, X, w, mu, sigma, lambda, beta) {
  Jf <- Jacobian.full.r(y[delta==1],X[delta==1,,drop=FALSE],w[delta==1],mu,beta,sigma,lambda)
  Jc <- Jacobian.cens.r(y[delta==0],X[delta==0,,drop=FALSE],w[delta==0],mu,beta,sigma,lambda)
  J <- (Jf*sum(delta==1) + Jc*sum(delta==0))
  solve(J)
}

WLoneCensloggammareg.expJ <- function(y,delta,x,w,mu0,bet0,sig0,lam0,qemp,step=1,d=100,reparam=FALSE) {
  .NotYetImplemented()
}

WLoneCensloggammareg.empJ <- function(y,delta,x,w,mu0,bet0,sig0,lam0,step=1,d=100,reparam=FALSE) {
  p <- ncol(x)
  n <- length(y)
  yo <- y[delta==1]
  yc <- y[delta==0]
  no <- sum(delta)
  nc <- n-no
  xo <- as.matrix(x[delta==1,])
  xc <- as.matrix(x[delta==0,])
  Uf <- Uscore.full.r(yo,xo,w[delta==1],mu0,bet0,sig0,lam0)
  Uc <- Uscore.cens.r(yc,xc,w[delta==0],mu0,bet0,sig0,lam0)
#U  <- (Uf*(no/n) + Uc*(nc/n))
  U  <- Uf*no + Uc*nc
  Jf <- Jacobian.full.r(yo,xo,w[delta==1],mu0,bet0,sig0,lam0)
  Jc <- Jacobian.cens.r(yc,xc,w[delta==0],mu0,bet0,sig0,lam0)

  eJf <- eigen(Jf)
  epos <- min(eJf$values[which(eJf$values > 0)])
  eJf$values[eJf$values < epos] <- epos
  if (abs(eJf$values[1]/eJf$values[p+3])>d) Jf <- eJf$vectors%*%diag(eJf$values+(eJf$values[1]-eJf$values[p+3]*d)/(d-1))%*%t(eJf$vectors)

  eJc <- eigen(Jc)
  epos <- min(eJc$values[which(eJc$values > 0)])
  eJc$values[eJc$values < epos] <- epos
  if (abs(eJc$values[1]/eJc$values[p+3])>d) Jc <- eJc$vectors%*%diag(eJc$values+(eJc$values[1]-eJc$values[p+3]*d)/(d-1))%*%t(eJc$vectors)
  J  <- Jf*no + Jc*nc
  tht <- c(mu0,bet0,sig0,lam0)
  if (reparam) {
    kappa <- 2*sig0^0.5
    U[p+2]  <- kappa*U[p+2]
    H <- diag(p+3)
    H[p+2,p+2] <- kappa
    J <- H%*%J%*%H 
    tht <- c(mu0,bet0,sqrt(sig0),lam0)
  }
  rank <- p+3
  err <- 0
  rank <- try(qr(J)$rank)
  if (rank == (p+3)) {
    eJ <- eigen(J)
    epos <- min(eJ$values[which(eJ$values > 0)])
    eJ$values[eJ$values < epos] <- epos
    if (abs(eJ$values[1]/eJ$values[p+3])>d) J <- eJ$vectors%*%diag(eJ$values+(eJ$values[1]-eJ$values[p+3]*d)/(d-1))%*%t(eJ$vectors)
    JI  <- qr.solve(J)
    tht <- tht-step*JI%*%U
  } else {
    err <- 1
    tht <- c(mu0,bet0,ifelse(reparam,sqrt(sig0),sig0),lam0)  
  }
  if (reparam) tht[p+2] <- tht[p+2]^2
  res <- list(mu=tht[1],beta=tht[2:(p+1)],sigma=tht[p+2],lambda=tht[p+3])
  return(res)
}
