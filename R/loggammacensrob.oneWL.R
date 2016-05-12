loggammacensrob.oneWL <- function(x, delta, start=NULL, w=rep(1, length(x)), control, expJ=FALSE, d=100) {
  control$d <- d
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
  qemp <- SemiQuant(y=x,delta=delta,mu0,sig0,lam0,nmod=control$subdivisions,
           tol=control$refine.tol) # 1000 quantili della semi-empirica basata sul WQtau
  w <- wi <- weightsCensloggamma(y=x,delta=delta,qemp=qemp,
               mu=mu0,sigma=sig0, lambda=lam0, bw=control$bw,
               raf=control$raf, tau=control$tau,nmod=control$subdivisions) # pesi per il WLone
  w[w > 1-control$minw] <- 1
  pos <- w >= control$minw
  x <- x[pos]
  delta <- delta[pos]
  w <- w[pos]
  if (expJ) {
    zWL <- WLoneCensloggamma.expJ(y=x,delta=delta,wi=w,mu0=mu0,sig0=sig0,
            lam0=lam0,qemp=qemp,step=control$step,d=control$d,reparam=FALSE)
    if (zWL$sigma < 0) {
      zWL <- WLoneCensloggamma.expJ(y=x,delta=delta,wi=w,mu0=mu0,sig0=sig0,
            lam0=lam0,qemp=qemp,step=control$step,d=control$d,reparam=TRUE)
      cat("Repar","\n")
    }
  } else {
    zWL <- WLoneCensloggamma.empJ(y=x,delta=delta,wi=w,mu0=mu0,sig0=sig0,
             lam0=lam0,step=control$step,d=control$d,reparam=FALSE) # WLone
    if (zWL$sigma < 0) {
      zWL <- WLoneCensloggamma.empJ(y=x,delta=delta,wi=w,mu0=mu0,sig0=sig0,
             lam0=lam0,step=control$step,d=control$d,reparam=TRUE)
      cat("Repar","\n")      
    }
  }
  result <- list(mu=zWL$mu, sigma=zWL$sigma, lambda=zWL$lambda, weights=wi, iterations=1, error=NULL)
  return(result)
}

WLoneCensloggamma.empJ <- function(y,delta,wi,mu0,sig0,lam0,step=1,d=100,reparam=FALSE) {
# 1SWL for lambda, mu, sigma
# This function uses the empirical Jacobian 
  n <- length(y)
  yo <- y[delta==1]
  yc <- y[delta==0]
  no <- sum(delta)
  nc <- n-no
  Uf <- Uscore.full(yo,wi[delta==1],mu0,sig0,lam0)
  Uc <- Uscore.cens(yc,wi[delta==0],mu0,sig0,lam0)
### U  <- (Uf*(no/n) + Uc*(nc/n))
  U  <- (Uf*no + Uc*nc)                                 
  Jf <- Jacobian.full(yo,wi[delta==1],mu0,sig0,lam0)
  Jc <- Jacobian.cens(yc,wi[delta==0],mu0,sig0,lam0)
### J  <- (Jf*(no/n) + Jc*(nc/n))
  J <- (Jf*no + Jc*nc)
  tht <- c(mu0,sig0,lam0)
  if (reparam) {
    kappa <- 2*sig0^0.5
    U[2]  <- kappa*U[2]
    H <- matrix(c(1,0,0,0,kappa,0,0,0,1),nrow=3,byrow=TRUE)
    J <- H%*%J%*%H 
    tht <- c(mu0,sqrt(sig0),lam0)
  }
  rank <- 3
  err <- 0
  rank <- try(qr(J)$rank)
#
  if (rank == 3) {
    eJ <- eigen(J)
    if (abs(eJ$values[1]/eJ$values[3])>d) J=eJ$vectors%*%diag(eJ$values+(eJ$values[1]-eJ$values[3]*d)/(d-1))%*%t(eJ$vectors)
    JI <- qr.solve(J)
    tht <- tht-step*JI%*%U
  } else {
    err <- 1
    tht <- c(mu0,ifelse(reparam,sqrt(sig0),sig0),lam0)
  }
  if (reparam) tht[2] <- tht[2]^2
  res <- list(mu=tht[1], sigma=tht[2], lambda=tht[3])
  return(res)
}

WLoneCensloggamma.expJ <- function(y,delta,wi,mu0,sig0,lam0,qemp,step=1,d=100,reparam=FALSE) {
# 1SWL for lambda, mu, sigma
# This function uses the expected Jacobian 
  n <- length(y)
  yo <- y[delta==1]
  yc <- y[delta==0]
  no <- sum(delta)
  nc <- n-no
  Uf <- Uscore.full(yo,wi[delta==1],mu0,sig0,lam0)
  Uc <- Uscore.cens(yc,wi[delta==0],mu0,sig0,lam0)
  U <- (Uf*(no/n) + Uc*(nc/n))
  no1 <- no; nc1 <- nc; yo1 <- yo; yc1 <- yc
  if (delta[n]==1) {
    nc1 =nc1+1
    no1 <- n-nc1
    yc1[nc1] <- max(y)
  }
  deltac <- delta
  deltac[n] <- 0
  SG <- survfit(Surv(y,(1-deltac))~1,error="tsiatis")
  Gn <- summary(SG)$surv       # survival
  Jn <- 1-Gn                   # cumulative
  gn <- diff(c(0,Jn))          # jumps
  x <- sort(yc1)
  y1 <- Jn
  y2 <- Jn-gn/2
  y2[nc1] <- Jn[nc1]-0.5/n
  xmin <- min(y)-0.2
  xmax <- max(y)+0.2
  Jf <- Exp.Jacobian.full(mu=mu0,sigma=sig0,lambda=lam0,x=x,y2=y2)
  Jc <- Exp.Jacobian.cens(mu=mu0,sigma=sig0,lambda=lam0,xmin,xmax,qemp,x=x,y2=y2)

  eJf <- eigen(Jf)
  epos <- min(eJf$values[which(eJf$values > 0)])
  eJf$values[eJf$values < epos] <- epos
  if (abs(eJf$values[1]/eJf$values[3])>d) Jf=eJf$vectors%*%diag(eJf$values+(eJf$values[1]-eJf$values[3]*d)/(d-1))%*%t(eJf$vectors)

  eJc <- eigen(Jc)
  epos <- min(eJc$values[which(eJc$values > 0)])
  eJc$values[eJc$values < epos] <- epos
  if (abs(eJc$values[1]/eJc$values[3])>d) Jc=eJc$vectors%*%diag(eJc$values+(eJc$values[1]-eJc$values[3]*d)/(d-1))%*%t(eJc$vectors)

  J <- Jf+Jc
  tht <- c(mu0,sig0,lam0)
  if (reparam) {
    kappa <- 2*sig0^0.5
    U[2] <- kappa*U[2]
    H <- matrix(c(1,0,0,0,kappa,0,0,0,1),nrow=3,byrow=TRUE)
    J <- H%*%J%*%H 
    tht <- c(mu0,sqrt(sig0),lam0)
  }
  rank <- 3
  err <- 0
  rank <- try(qr(J)$rank)
  if (rank == 3) {
    eJ <- eigen(J)
    epos <- min(eJ$values[which(eJ$values > 0)])
    eJ$values[eJ$values < epos] <- epos
    if (abs(eJ$values[1]/eJ$values[3])>d)
      J <- eJ$vectors%*%diag(eJ$values+(eJ$values[1]-eJ$values[3]*d)/(d-1))%*%t(eJ$vectors)
    JI <- qr.solve(J)
    tht <- tht-step*JI%*%U
  } else {
    err <- 1
    tht <- c(mu0,ifelse(reparam,sqrt(sig0),sig0),lam0)
  }
  if (reparam)
    tht[2] <- tht[2]^2
  res <- list(mu=tht[1], sigma=tht[2], lambda=tht[3])
  return(res)
}

VCOV.WLoneCensloggamma.empJ <- function(y, delta, w, mu, sigma, lambda) {
  Jf <- Jacobian.full(y[delta==1],w[delta==1],mu,sigma,lambda)
  Jc <- Jacobian.cens(y[delta==0],w[delta==0],mu,sigma,lambda)
  J <- (Jf*sum(delta==1) + Jc*sum(delta==0))
  solve(J)
}

VCOV.WLoneCensloggamma.expJ <- function(y, delta, w, mu, sigma, lambda, nmod=1000, tol=1e-4) {
  qemp <- SemiQuant(y=y,delta=delta,mu,sigma,lambda,nmod=nmod,
                    tol=tol) # 1000 quantili della semi-empirica basata sul WQtau
  n <- length(y)
  yo <- y[delta==1]
  yc <- y[delta==0]
  no <- sum(delta)
  nc <- n-no
  no1 <- no; nc1 <- nc; yo1 <- yo; yc1 <- yc
  if (delta[n]==1) {
    nc1 =nc1+1
    no1 <- n-nc1
    yc1[nc1] <- max(y)
  }
  deltac <- delta
  deltac[n] <- 0
  SG <- survfit(Surv(y,(1-deltac))~1,error="tsiatis")
  Gn <- summary(SG)$surv       # survival
  Jn <- 1-Gn                   # cumulative
  gn <- diff(c(0,Jn))          # jumps
  x <- sort(yc1)
  y1 <- Jn
  y2 <- Jn-gn/2
  y2[nc1] <- Jn[nc1]-0.5/n
  xmin <- min(y)-0.2
  xmax <- max(y)+0.2
  Jf <- Exp.Jacobian.full(mu=mu,sigma=sigma,lambda=lambda,x=x,y2=y2)
  Jc <- Exp.Jacobian.cens(mu=mu,sigma=sigma,lambda=lambda,xmin,xmax,qemp,x=x,y2=y2)
  J <- Jf+Jc
  solve(J)
}


J.trovazero0 <- function(xx, xmin, xmax, p, tol=10^(-8),x=x,y2=y2) {
  Js0 <- splinefun(x,y2, method= "monoH.FC")
  Js  <- function(x){ js0 <- Js0(x); pmax(0, pmin(1,js0))} 
  I <- ceiling(-log(tol,2)); i <- 1; step <- xx
  differenza <- xold <- xx + 1 + tol
  intervallo <- xmax-xmin
  neverdown <- TRUE
  while (i <= I & differenza > tol) {
    xold <- xx
    loss <- Js(xx*intervallo+xmin)
    if (loss < p) {
      if (neverdown) {
        step <- xx <- min(2*xx,1)
      } else {
        xx <- xx+step/2^i
        i <- i + 1}
    } else {
      neverdown <- FALSE
      xx <- xx-step/2^i
      i <- i + 1}
    differenza <- max(abs(xx-xold))}
  xx <- xx*intervallo+xmin
  return(xx)
}

J.Quant0 <- function(nmod=1000,xmin,xmax,tol=1e-4,x=x,y2=y2) {  # Quantiles of J
  p <- qJemp <- ppoints(nmod)
  qJemp[1] <- J.trovazero0(xx=p[1],xmin=xmin,xmax=xmax,p=p[1],tol=10^(-8),x=x,y2=y2)
  qJemp[nmod] <- J.trovazero0(xx=p[nmod], xmin=xmin, xmax=xmax, p=p[nmod],tol=10^(-8),x=x,y2=y2)
  for (i in 2:(nmod-1)) {
    qJemp[i] <- J.trovazero0(xx=p[i],xmin=qJemp[i-1]-0.1,xmax=qJemp[nmod]+0.1,p=p[i],tol=10^(-8),x=x,y2=y2)
  }
  return(qJemp)
}

#####################################################
#####################################################

weightsCensloggamma <- function(y,delta,qemp,mu,sigma,lambda,bw=0.3,raf="NED",tau=0.5,nmod=1000) {
  ## delta non e' usato per ora!  
  r <- (y-mu)/sigma
  Semp  <- density(qemp, kernel="gaussian", bw=bw, cut=3,n=512)
  semp  <- approxfun(x=Semp$x, y=Semp$y, rule=2)
  pe    <- semp(qemp)
  qmodl <- qloggamma(ppoints(nmod),lambda=lambda)
  Smod  <- density(qmodl, kernel="gaussian", bw=bw, cut=3,n=512)
  smod  <- approxfun(x=Smod$x, y=Smod$y, rule=2)
  pm    <- smod(qemp)
  pe    <- semp(r)
  pm    <- smod(r)
  Delta <- pe/pm-1
  w     <- pesi(x=Delta, raf=raf, tau=tau)
  return(w)
}

SurvQuant <- function(y,delta,w) {
# KM survival
  S1  <- survfit(Surv(y,delta)~1,weights=w,error="tsiatis")# summary(S1)
  ync <- summary(S1)$time   # ordered noncensored obervations
  Sn  <- summary(S1)$surv   # survival
  Fn  <- 1-Sn               # cumulative
  pn  <- diff(c(0,Fn))      # jumps
  qi  <- Fn-pn/2            # quantiles of KM
  ans <- list(ri=ync,qi=qi,S1=S1)
  return(ans)
}

SemiQuant <- function(y,delta,mu,sigma,lambda,nmod=1000,tol=1e-4) {
  p <- qemp <- ppoints(nmod)
  r <- (y-mu)/sigma
  minr <- min(r) - 1
  maxr <- max(r) + 1
  qemp[1] <- QnsemiGLG(x=p[1],xmin=minr,xmax=maxr,p=p[1],y=r,delta=delta,mu=0,sigma=1,lambda=lambda,tol=tol)
  qemp[nmod] <- QnsemiGLG(x=p[nmod],xmin=minr,xmax=maxr,p=p[nmod],y=r,delta=delta,mu=0,sigma=1,lambda=lambda, tol=tol)
  for (i in 2:(nmod-1)) {
    qemp[i] <- QnsemiGLG(x=p[i],xmin=qemp[i-1]-0.1,xmax=qemp[nmod]+0.1,p=p[i],y=r,delta=delta,mu=0,sigma=1,lambda=lambda,tol=tol)
  }
  return(qemp)
}
