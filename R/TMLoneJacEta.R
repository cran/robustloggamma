# TMLone_loggamma_cens-JacEta.s
# =============================

# Note: rso is the standardozed residual vector wrt the inital values; 
#       rsd is the standardized residual vector wrt the current estimate

#--- jac1*

TMLjac11.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of Beta-equation wrt Beta.hat
  n <- length(rs0)
  p <- ncol(X) 
  zero <- 1e-9  
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(q=rsd,lambda=lambda)
  fo    <- dloggamma(x=rsd,lambda=lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl)-xbi)/d.sigma
  bi    <- (cu  -        xbi)/d.sigma
  foai  <- dloggamma(x=ai,lambda=lambda)
  fobi  <- dloggamma(x=bi,lambda=lambda)
  fopai <- p.loggamma(ai,lambda=lambda)  # foai*(-ps0LG(ai,lambda)) 
  fopbi <- p.loggamma(bi,lambda=lambda)  # fobi*(-ps0LG(bi,lambda))
  f1[ok]<- fo[ok]/(1-Fo[ok])^2
  f2[ok]<-      1/(1-Fo[ok])
  D1    <- - delta*wgt*csiLG.prime(rsd,lambda)
  D2    <- - (  (1-delta)*f1*(fobi-foai )    )
  D3    <- - (  (1-delta)*f2*( fopbi-fopai ) )
  D     <- D1 + D2 + D3
  Jac   <- t(X) %*% (as.vector(D)*X) /d.sigma/n
  return(Jac)
}

TMLjac12.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of Beta-equation wrt sigma.hat
  n <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(q=rsd,lambda=lambda)
  fo    <- dloggamma(x=rsd,lambda=lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -          xbi)/d.sigma
  foai  <- dloggamma(x=ai,lambda=lambda)
  fobi  <- dloggamma(x=bi,lambda=lambda)
  fopai <- p.loggamma(ai,lambda=lambda)
  fopbi <- p.loggamma(bi,lambda=lambda)
  f1[ok]<- fo[ok]/(1-Fo[ok])^2
  f2[ok]<-      1/(1-Fo[ok])
  D1    <- -  delta*wgt*csiLG.prime(rsd,lambda)*rsd
  D2    <- - (  (1-delta)*f1*(fobi-foai)*rsd  )
  D3    <- - (  (1-delta)*f2*(fopbi*bi-fopai*ai )  )
  D     <- D1 + D2 + D3
  Jac   <- t(X)%*%(as.vector(D))/d.sigma/n
  return(Jac)
}

TMLjac13.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of Beta-equation wrt lambda
  n <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(q=rsd,lambda=lambda)
  Fod   <- d.FL.1(rsd,lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -          xbi)/d.sigma
  foai  <- dloggamma(x=ai,lambda=lambda)
  fobi  <- dloggamma(x=bi,lambda=lambda)
  fodai <- d.fL(ai,lambda) 
  fodbi <- d.fL(bi,lambda)
  f1[ok]<- Fod[ok]/(1-Fo[ok])^2
  f2[ok]<-       1/(1-Fo[ok])
  D1    <- delta*wgt*csiLG.dot(rsd,lambda)
  D2    <- (1-delta)*f1*(fobi - foai)
  D3    <- (1-delta)*f2*(fodbi- fodai)
  D     <- D1 + D2 + D3
  Jac   <- t(X)%*%(as.vector(D))/d.sigma/n
  return(Jac)
}

#--- jac2*

TMLjac21.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of sigma-equation wrt Beta.hat
  n <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(q=rsd,lambda=lambda)
  fo    <- dloggamma(x=rsd,lambda=lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -        - xbi)/d.sigma
  foai  <- dloggamma(x=ai,lambda=lambda)
  fobi  <- dloggamma(x=bi,lambda=lambda)
  Foai  <- ploggamma(q=ai,lambda=lambda)
  Fobi  <- ploggamma(q=bi,lambda=lambda)
  fopai <- p.loggamma(ai,lambda=lambda)
  fopbi <- p.loggamma(bi,lambda=lambda)
  f1[ok]<- fo[ok]/(1-Fo[ok])^2
  f2[ok]<-      1/(1-Fo[ok])
  D1    <- - delta*wgt*phiLG.prime(rsd,lambda)
  D2    <- - (  (1-delta)*f1*(fobi*bi   - foai*ai + Foai - Fobi))
  D3    <- - (  (1-delta)*f2*( fopbi*bi - fopai*ai ) )
  D     <- D1 + D2 + D3
  Jac   <- t(X)%*%(as.vector(D))/d.sigma/n
  return(Jac)
}

TMLjac22.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of sigma-equation wrt sigma.hat
  n     <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n) 
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(q=rsd,lambda=lambda)
  fo    <- dloggamma(x=rsd,lambda=lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -          xbi)/d.sigma
  foai  <- dloggamma(ai,lambda=lambda)
  fobi  <- dloggamma(bi,lambda=lambda)
  Foai  <- ploggamma(ai,lambda=lambda)
  Fobi  <- ploggamma(bi,lambda=lambda)
  fopai <- p.loggamma(ai,lambda=lambda)
  fopbi <- p.loggamma(bi,lambda=lambda)
  f1[ok]<- fo[ok]/(1-Fo[ok])^2
  f2[ok]<-      1/(1-Fo[ok])
  D1    <- - delta*wgt*phiLG.prime(rsd,lambda)*rsd
  D2    <- -(  (1-delta)*f1*(fobi*bi    - foai*ai + Foai - Fobi)*rsd  )
  D3    <- -(  (1-delta)*f2*(fopbi*bi^2 - fopai*ai^2 )                )
  D     <- D1 + D2 + D3
  Jac   <- sum(D)/d.sigma/n
  return(Jac)
}

TMLjac23.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of sigma-equation wrt lambda
  n     <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9  
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(rsd,lambda=lambda)
  fo    <- dloggamma(rsd,lambda=lambda)
  Fod   <- d.FL.1(rsd,lambda)
  fod   <- d.fL(rsd,lambda) 
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -          xbi)/d.sigma
  foai  <- dloggamma(ai,lambda=lambda)
  fobi  <- dloggamma(bi,lambda=lambda)
  Foai  <- ploggamma(ai,lambda=lambda)
  Fobi  <- ploggamma(bi,lambda=lambda)
  fodai <- d.fL(ai,lambda) 
  fodbi <- d.fL(bi,lambda)
  Fodai <- d.FL.1(ai,lambda)
  Fodbi <- d.FL.1(bi,lambda)
  f1[ok]<- Fod[ok]/(1-Fo[ok])^2
  f2[ok]<-       1/(1-Fo[ok])
  D1    <- delta*wgt*phiLG.dot(rsd,lambda)
  D2    <- (1-delta)*f1*( bi*fobi - ai*foai -  Fobi  +  Foai)
  D3    <- (1-delta)*f2*( bi*fodbi -ai*fodai - Fodbi + Fodai)
  D     <- D1 + D2 + D3
  Jac   <- sum(D)/d.sigma/n
  return(Jac)
}

#--- jac3*

TMLjac31.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of lambda-equation wrt Beta.hat
  n  <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(rsd,lambda=lambda)
  fo    <- dloggamma(rsd,lambda=lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -          xbi)/d.sigma
  fodai <- d.fL(ai,lambda) 
  fodbi <- d.fL(bi,lambda)
  Fodai <- d.FL.1(ai,lambda)
  Fodbi <- d.FL.1(bi,lambda)
  f1[ok]<- fo[ok]/(1-Fo[ok])^2
  f2[ok]<-      1/(1-Fo[ok])
  D1    <- - delta*wgt*psiLG.prime(rsd,lambda)
  D2    <- - (  (1-delta)*f1*(Fodai - Fodbi)  )
  D3    <- - (  (1-delta)*f2*(fodai - fodbi)  )
  D     <- D1 + D2 + D3
  Jac   <- t(X)%*%(as.vector(D))/d.sigma/n
  return(Jac)
}

TMLjac32.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of lambda-equation wrt sigma.hat
  n <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9  
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  xbi   <- X%*%d.Beta
  rsd   <- (rs0 - xbi)/d.sigma
  Fo    <- ploggamma(rsd,lambda=lambda)
  fo    <- dloggamma(rsd,lambda=lambda)
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0,cl) - xbi)/d.sigma
  bi    <- (cu  -          xbi)/d.sigma
  fodai <- d.fL(ai,lambda) 
  fodbi <- d.fL(bi,lambda)
  Fodai <- d.FL.1(ai,lambda)
  Fodbi <- d.FL.1(bi,lambda)
  f1[ok]<- fo[ok]/(1-Fo[ok])^2
  f2[ok]<-      1/(1-Fo[ok])
  D1    <- - delta*wgt*psiLG.prime(rsd,lambda)*rsd
  D2    <- -(  (1-delta)*f1*( Fodai   -  Fodbi  )*rsd  ) 
  D3    <- -(  (1-delta)*f2*(ai*fodai - bi*fodbi)      )
  D     <- D1 + D2 + D3
  Jac   <- sum(D)/d.sigma/n
  return(Jac)
}

int7  <- function(x,lambda) {
  (psiLG.dot(x,lambda)-psiLG(x,lambda)^2)*dloggamma(x=x,lambda=lambda)
}

F2dab <- function(i,lambda,ai,bi) {
  integrate(int7, ai[i], bi[i], lambda=lambda)$value
}

TMLjac33.LG <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) {
# Jacobian of lambda-equation wrt lambda
  n    <- length(rs0); D1 <- D2 <- D3 <- D <- 0;  zero <- 1e-9
  indu <- (1:n)[delta==1]; indc <- (1:n)[delta==0]; nc <- sum(delta==0); f1 <- f2 <- rep(0,nc)
  xbi  <- X%*%d.Beta
rsd  <- (rs0-xbi)/d.sigma
  if (nc < n)
    D1 <- sum( wgt[indu]*psiLG.dot(rsd[indu],lambda) )
  if (nc > 0) {
    rsdc  <- rsd[indc]
  Fo    <- ploggamma(q=rsdc,lambda=lambda)
  fo    <- dloggamma(x=rsdc,lambda=lambda)
  Fod   <- d.FL.1(rsdc,lambda)
  fod   <- d.fL(rsdc,lambda) 
  ok    <- (1-Fo) > zero
  ai    <- (pmax(rs0[indc],cl)-xbi[indc])/d.sigma
  bi    <- (cu  -              xbi[indc])/d.sigma
  Fodai <- d.FL.1(ai,lambda)
  Fodbi <- d.FL.1(bi,lambda)
  Fod2.diff <- apply( as.matrix(1:nc), MARGIN=1, FUN=F2dab, lambda=lambda, ai=ai, bi=bi)
  f1[ok]<- Fod[ok]/(1-Fo[ok])^2
  f2[ok]<-       1/(1-Fo[ok])
  D2  <- sum( f1*( Fodai  - Fodbi) )
  D3  <- sum( f2*( Fod2.diff )    ) }
  D  <- D1 + D2 + D3
  Jac   <- D/d.sigma/n
  return(Jac)
}

TMLjac33.LG.0 <- function(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu) { # this version is very slow !
# Jacobian of lambda-equation wrt lambda
  n <- length(rs0)
  p <- ncol(X)
  zero <- 1e-9  
  D1 <- D2 <- D3 <- D <- f1 <- f2 <- rep(0,n)
  rsd    <- (rs0-X%*%d.Beta)/d.sigma
  Fo     <- ploggamma(q=rsd,lambda=lambda)
  fo     <- dloggamma(x=rsd,lambda=lambda)
  Fod    <- d.FL.1(rsd,lambda)
  fod    <- d.fL(rsd,lambda) 
  ok     <- (1-Fo) > zero & Fod != 0
  ai     <- (pmax(rs0,cl)-X%*%d.Beta)/d.sigma
  bi     <- (cu  -        X%*%d.Beta)/d.sigma
  Fodai  <- d.FL.1(ai,lambda)
  Fodbi  <- d.FL.1(bi,lambda)
  Fod2ai <- apply(ai,1,d2.FL,lambda=lambda)
  Fod2bi <- apply(bi,1,d2.FL,lambda=lambda)
  f1[ok] <- Fod[ok]/(1-Fo[ok])^2
  f2[ok] <-       1/(1-Fo[ok])
  D1     <- delta*wgt*psiLG.dot(rsd,lambda)
  D2     <- (1-delta)*f1*( Fodai  - Fodbi  )
  D3     <- (1-delta)*f2*( Fod2ai - Fod2bi )
  D      <- D1 + D2 + D3
  Jac    <- sum(D)/d.sigma/n
  return(Jac)
}

JAC.TML.loggamma <- function(X,yo,delta,Beta,Beta.t,sigma,sigma.t,lambda,rs0,wgt,cl,cu,kk,gg) {
  p <- ncol(X); JAC <- matrix(0, nrow=(p+2), ncol=(p+2))
  d.Beta  <- (Beta - Beta.t)/sigma.t
  d.sigma <- sigma/sigma.t
  J11 <- TMLjac11.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J12 <- TMLjac12.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J13 <- TMLjac13.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J21 <- TMLjac21.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J22 <- TMLjac22.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J23 <- TMLjac23.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J31 <- TMLjac31.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J32 <- TMLjac32.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  J33 <- TMLjac33.LG(d.Beta,d.sigma,lambda,rs0,wgt,delta,X,cl,cu)
  JAC[1:p,1:p]     <- J11
  JAC[1:p,(p+1)]   <- J12
  JAC[1:p,(p+2)]   <- J13
  JAC[(p+1),1:p]   <- J21
  JAC[(p+1),(p+1)] <- J22
  JAC[(p+1),(p+2)] <- J23
  JAC[(p+2),1:p]   <- J31
  JAC[(p+2),(p+1)] <- J32
  JAC[(p+2),(p+2)] <- J33
  JAC <- JAC/sigma.t
  return(JAC)
}

#--- Eta ----------------------------------------------------------------------

# Eta.s

Eta.LG <- function(X, yo, delta, sigma, sigma.t, mui, mui.t, lambda, rs0, wgt, cl, cu, kk,gg) {
p <- ncol(X);  n <- length(rs0); zero <- 1e-9
iu  <- delta==1; ic <- delta==0;  nc <- sum(ic); nu <- n-nc
Etau <- Etac <- matrix(0,nrow=(p+2),ncol=1)
Hi <- EHi <- rep(0,p); ki=gi=eki=egi=0
#
rsd  <- (yo - mui)/sigma
#
if(nc < n) {
 Xu  <- X[iu,,drop=FALSE]
 hi  <- wgt[iu]*csiLG(rsd[iu],lambda)
 Hi  <- colSums( as.vector(hi)*Xu )/n
 ki  <- sum(wgt[iu]*phiLG(rsd[iu],lambda))/n
 gi  <- sum(wgt[iu]*psiLG(rsd[iu],lambda))/n}
#
if (nc > 0) {
 ehi <- eki <- egi <- rep(0,nc)
 Xc <- X[ic,,drop=FALSE]
 muic <- mui[ic]
 muit <- mui.t[ic]
 rci <- rsd[ic]
 ai   <- pmax( rci, (sigma.t*cl-muic+muit)/sigma )
 bi    <-           (sigma.t*cu-muic+muit)/sigma
 Fo   <- ploggamma(q=rci,lambda=lambda)
 ok   <- (1-Fo) > zero  
 Ihi  <- dloggamma(x=bi,lambda=lambda) - dloggamma(x=ai,lambda=lambda)
 Iki  <- bi * dloggamma(x=bi,lambda=lambda) - ai * dloggamma(x=ai,lambda=lambda) + ploggamma(q=ai,lambda=lambda) - ploggamma(q=bi,lambda=lambda)
 Igi  <- d.FL.1(ai,lambda) - d.FL.1(bi,lambda)
 ehi[bi>ai & ok]  <- Ihi[bi>ai & ok]/(1-Fo[bi>ai & ok])
 eki[bi>ai & ok]  <- Iki[bi>ai & ok]/(1-Fo[bi>ai & ok])
 egi[bi>ai & ok]  <- Igi[bi>ai & ok]/(1-Fo[bi>ai & ok])
 EHi <- colSums ( as.vector(ehi)*Xc )/n
 eki <- sum(eki)/n
 egi <- sum(egi)/n }
Eta =c( (Hi+EHi), (ki+eki - kk), (gi+egi - gg) )
  return(Eta)
}

#--- TMLone -------------------------------------------------------------------

TMLone.loggamma <- function(IX,y,delta,Beta.t,sigma.t,lambda.t,cl,cu,step=1,reparam=FALSE) {
  if (missing(cl) | is.na(cl) | missing(cu) | is.na(cu)) {
    err <- 2
    n.ret <- NA
    theta <- c( Beta.t, sigma.t, lambda.t)
  } else {
    p <- ncol(IX) # Beta.ini   <- Beta.t; sigma.ini  <- sigma.t
    mui <- mui.t  <- IX %*% as.matrix(Beta.t)
    rs0    <- (y-mui.t)/sigma.t
# right hand sides
    hh <-  0
    kk <- -cl*dloggamma(x=cl,lambda=lambda.t) + cu*dloggamma(x=cu,lambda=lambda.t)-ploggamma(q=cu,lambda=lambda.t)+ploggamma(q=cl,lambda=lambda.t)
gg <-  d.FL.1(cl,lambda.t)-d.FL.1(cu,lambda.t)
# weights
    w <- wgt <- weightsHD(rs0, cl, cu)
#
# estimate
#
    pos <- wgt > 0.001
    IX <- IX[pos,,drop=FALSE]
    y <- y[pos]
    delta <- delta[pos]
    mui <- mui[pos]
    mui.t <- mui.t[pos]
    wgt <- wgt[pos]
    rs0 <- rs0[pos]

    Eta <- Eta.LG(IX, y, delta, sigma.t, sigma.t, mui, mui.t, lambda.t, rs0, wgt, cl, cu, kk, gg)
   JAC <- JAC.TML.loggamma(IX, y, delta, Beta.t, Beta.t, sigma.t, sigma.t, lambda.t, rs0, wgt, cl, cu, kk, gg)
   theta.ini  <- c(Beta.t,sigma.t,lambda.t)
#
    if (reparam==TRUE) {
      theta.ini  <- c(Beta.t,sqrt(sigma.t),lambda.t)
      kappa <- 2*sigma.t^0.5
      Eta[p+1] <- Eta[p+1]*kappa  # sigma component
      G <- diag(p+2)
      G[p+1,p+1] <- kappa
      JAC <- G%*%JAC%*%G
    }
#
    J=JAC; d <- 100; rank=p+2; err=0
    rank <- try(qr(J,LAPACK=TRUE)$rank)
    if (rank == (p+2)) {
      eJ <- eigen(J)
      if (abs(eJ$values[1]/eJ$values[p+2])>d) J=eJ$vectors%*%diag(eJ$values+(eJ$values[1]-eJ$values[p+2]*d)/(d-1))%*%t(eJ$vectors)
      JIC <- QRsolve(J)
      theta <- theta.ini-step*JIC%*%Eta
    } else {
      err <- 1
      theta <- c( Beta.t, sigma.t, lambda.t)
    }
    if (reparam==TRUE) theta[p+1] <- theta[p+1]^2
    n.ret <- sum(wgt)
  }
  res <- list(Beta=Re(theta[1:p]), sigma=Re(theta[p+1]), lambda=Re(theta[p+2]), n.ret=n.ret, error=err, weights=w)
  return(res)
}

QRsolve <- function (a, tol=1e-07) {
    if (!is.qr(a)) 
        a <- qr(a, tol <- tol, LAPACK=TRUE)
    nc <- ncol(a$qr)
    nr <- nrow(a$qr)
    if (a$rank != min(nc, nr)) 
        stop("input singular matrix in solve")
    if (nc != nr) stop("only square matrices can be inverted")
    b <- diag(1, nc)
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    return(res)
}

#--- Matrix Q for covariance --------------------------------------------------------

QMatrix.LG <- function(d.Beta,d.sigma,lambda,rs0,delta,wgt,X,cl,cu,kk,gg) {
p <- ncol(X);  n  <- length(rs0); zero <- 1e-9
Hu <- Hc <- matrix(0,nrow=(p+2),ncol=(p+2))
iu  <- delta==1; ic <- delta==0;  nc <- sum(ic); nu <- n-nc
if(nc < n) {
 Xu  <- X[iu,,drop=FALSE]
 r0u <- rs0[iu]
 ru  <- (r0u - Xu%*%d.Beta)/d.sigma
 hi  <- wgt[iu]*csiLG(ru,lambda)
 ki  <- wgt[iu]*phiLG(ru,lambda) - kk
 gi  <- wgt[iu]*psiLG(ru,lambda) - gg
 Hu[1:p,1:p]      <- t(Xu)%*% (as.vector(hi^2) * Xu)
 Hu[1:p,(p+1)]    <- t(Xu)%*%as.matrix(hi*ki)
 Hu[1:p,(p+2)]    <- t(Xu)%*%as.matrix(hi*gi)
 Hu[(p+1),(p+1)]  <- sum(ki^2)
 Hu[(p+1),(p+2)]  <- sum(ki*gi)
 Hu[(p+2),(p+2)]  <- sum(gi^2)
 Hu[(p+1),1:p]    <- Hu[1:p,(p+1)] 
 Hu[(p+2),1:p]    <- Hu[1:p,(p+2)] 
 Hu[(p+2),(p+1)]  <- Hu[(p+1),(p+2)] }
if(nc > 0) {
 Xc <- X[ic,,drop=FALSE]
 r0c <- rs0[ic]
 rc  <- (r0c - Xc%*%d.Beta)/d.sigma
 Fo  <- ploggamma(q=rc,lambda=lambda)
 ok  <- (1-Fo) > zero
 ai  <- (pmax(r0c,cl)- Xc%*%d.Beta)/d.sigma
 bi  <- (         cu - Xc%*%d.Beta)/d.sigma
 Ihi <- dloggamma(x=bi,lambda=lambda) - dloggamma(x=ai,lambda=lambda)
 Iki <- bi * dloggamma(x=bi,lambda=lambda) - ai * dloggamma(x=ai,lambda=lambda) + ploggamma(q=ai,lambda=lambda) - ploggamma(q=bi,lambda=lambda)
 Igi <- d.FL.1(ai,lambda) - d.FL.1(bi,lambda)
 ehi <- eki <- egi <- rep(0,nc)
 ehi[bi>ai & ok]  <- Ihi[bi>ai & ok]/(1-Fo[bi>ai & ok])
 eki[bi>ai & ok]  <- Iki[bi>ai & ok]/(1-Fo[bi>ai & ok]) - kk
 egi[bi>ai & ok]  <- Igi[bi>ai & ok]/(1-Fo[bi>ai & ok]) - gg
 Hc[1:p,1:p]      <- t(Xc)%*% (as.vector(ehi^2) * Xc)
 Hc[1:p,(p+1)]    <- t(Xc)%*%as.matrix(ehi*eki)
 Hc[1:p,(p+2)]    <- t(Xc)%*%as.matrix(ehi*egi)
 Hc[(p+1),(p+1)]  <- sum(eki^2)
 Hc[(p+1),(p+2)]  <- sum(eki*egi)
 Hc[(p+2),(p+2)]  <- sum(egi^2)
 Hc[(p+1),1:p]    <- Hc[1:p,(p+1)] 
 Hc[(p+2),1:p]    <- Hc[1:p,(p+2)] 
 Hc[(p+2),(p+1)]  <- Hc[(p+1),(p+2)] }
 res <- (Hu+Hc)/n
 return(res)
}

#--- just for checking --------------------------------------------------------

TML.Ave2LG.lambda <- function(y,delta,sigma, sigma.t, mui, mui.t,wgt, cl,cu, lambda) {
  n <- length(y); ku <- kc <- 0; zero <- 1.e-9
  ku   <- kc <- 0; zero <- 1.e-9
  indu <- (1:n)[delta==1]; indc <- (1:n)[delta==0]; nc <- sum(delta==0); Iki <- eki <- rep(0,nc)
  rs <- (y-mui)/sigma
  if (nc < n) {
    ki  <- wgt[indu]*psiLG(rs[indu],lambda)
    ku  <- sum(ki)
  }
  if (nc > 0) {
    muic <- mui[indc]
    muit <- mui.t[indc]
    rci <- rs[indc]
    ai   <- pmax( rci, (sigma.t*cl-muic+muit)/sigma )
    bi    <-           (sigma.t*cu-muic+muit)/sigma 
    Fo   <- ploggamma(q=rci,lambda=lambda)
    ok   <- (1-Fo) > zero 
    for (i in 1:nc)
      Iki[i] <- d.FL.1(ai[i],lambda) - d.FL.1(bi[i],lambda)
    eki  <- rep(0,nc)
    # check
    eki[bi > ai & ok] <- Iki[bi > ai & ok]/(1-Fo[bi > ai & ok])
    # check
    kc <- sum(eki)
  }
  res <- (ku+kc)/n
  return(res)
}

TML.Ave2LG.sigma <- function(X, y, delta, sigma, sigma.t, mui, mui.t, wgt, cl, cu, lambda){
# equation left-hand sides for delta-sigma using a rectangular weight function
  n <- length(y); p <- ncol(X); ku <- kc <- 0; zero <- 1.e-9
  indu <- (1:n)[delta==1]; indc <- (1:n)[delta==0]; nc <- sum(delta==0)
  rs <- (y-mui)/sigma
  if (nc < n) {
    ki <- wgt[indu]* phiLG(rs[indu],lambda)
    ku <- sum(ki)
  }
  if (nc > 0) {
    muic <- mui[indc]
    muit <- mui.t[indc]
    rci  <- rs[indc]
    ai   <- pmax( rci, (sigma.t*cl-muic+muit)/sigma )
    bi   <-            (sigma.t*cu-muic+muit)/sigma 
    Fo   <- ploggamma(q=rci,lambda=lambda)
    ok   <- (1-Fo) > zero 
    Iki  <- -ai*dloggamma(x=ai,lambda=lambda)+bi*dloggamma(x=bi,lambda=lambda)-ploggamma(q=bi,lambda=lambda)+ploggamma(q=ai,lambda=lambda)
    eki  <- rep(0,nc)
  # check
    eki[bi > ai & ok] <- Iki[bi > ai & ok]/(1-Fo[bi > ai & ok])
  # check
    kc <- sum(eki)
  }
  res <- (ku+kc)/(n-p)
  return(res)
}
