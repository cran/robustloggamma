# FunctionsNP.R

# lint: Gauss: lint=0 (nonpar) 
# ialg: 3:bisection 1:fixed point 2:regula falsi
# meth: 1: S  4: sn

#######################
ResExpG <- function(r, tol=1e-5) {
  num <- dnorm(r)
  den <- 1-pnorm(r)
  val <- num/den
  val <- ifelse(den < tol, r, val)
  return(val)
}

ResExpNP <- function(r, delta) {
  n <- length(r)
  delta[which.max(r)] <- 1
  nu <- sum(delta)
  if (nu < n) {
    rc <- r[delta==0]
    km <- kpmr(r,delta)
    pt <- km$pt
    tu <- km$yi
    rc.exp <- apply(as.matrix(rc), 1, SuKM, tu=tu, pt=pt)
    r[delta==0] <- rc.exp
  }
  return(r)
}

SuKM <- function(r, tu, pt) {
  u <- tu
  isum <- tu > r
  s <- sum(u[isum]*pt[isum])
  den <- 1 - F.KM(r,tu,pt)
  s/den
} 

BtamatG <- function(X, y, delta, N, q, seed=153, cint=c("Param","Nonpar"), maxit=100, tol=0.001, ialg=3, bb=0.5, xk=1.5477) {     
# Function for determining the matrix of betas
  cint <- match.arg(cint)
  n <- length(y)
  p <- ncol(X)
  if (q < p)
    q <- p
  set.seed(seed)
  dlt1 <- (1:n)[delta == 1]
  ind <- apply(matrix(rep(dlt1, N), nrow = N, byrow = TRUE), 1, sample, size = q)
  intcp <- any(X[,1,drop=TRUE]!= 1)
##  if (intcp)
##    X <- cbind(1,X)
  beta <- apply(ind, 2, CandidateG, X, y, delta, cint, maxit, tol, ialg, bb, xk)
  if (p == 1)
    beta <- matrix(beta, ncol = 1, nrow = N)
  else
    beta <- t(beta)
  list(beta = beta)
}

CandidateG <- function(ind, X, y, delta, cint=c("Param","Nonpar"), maxit=100, tol=0.001, ialg=3, bb=0.5, xk=1.5477) {
  cint <- match.arg(cint)  
  sigma0 <- 1
  s0 <- 1
  n <- length(y)
  p <- ncol(X)
  options(warn=1)
  beta0 <- lsfit(X[ind,],y[ind],intercept=FALSE)$coef
  options(warn=0)
  if (cint=="Param") {
    sini <- s.eq.Gauss(X,y,N=1,delta,sigma0,bb,beta0,maxit,tol,s0=s0,ipsi=4,xk=xk,lint=1,ialg=ialg,meth=4)$S
    rs <- as.vector((y-X%*%as.matrix(beta0))/sini)
    Ind <- (1:n)[delta==0]
    rs0 <- rs[Ind]
    ur0   <- unique(rs0)
    if (length(ur0)>0) {
      uIxfd <- ResExpG(ur0)
      jjj <- match(rs0,ur0,nomatch=0)
      rs[Ind] <- uIxfd[jjj]
    }
    cu <- 0.6745
    input <- list(lambda=beta0,sigma=sini)
    yy <- X%*%as.matrix(beta0) + rs*sini
    z <- TML.gauss(X,yy,cu=cu,initial="input",otp="fixed",cov ="no",input=input,iv=1,nrep=0,tol=0.0001,seed=1313)
    coef  <- z$th1
  } else {
    r <- as.vector( (y-X%*%as.matrix(beta0)) )
    rs <- ResExpNP(r,delta)
    Ind <- order(abs(rs))[1:(n/2)]
    if (p==1)
      sres <- survreg(Surv(y,delta)~1, dist="gaussian", subset=Ind, control=list(maxiter=100, outer.max=20))
    else
      sres <- survreg(Surv(y,delta)~X[,-1], dist="gaussian", subset=Ind, control=list(maxiter=100, outer.max=20))
    coef <- sres$coef
  }
  return(coef)
}

# Auxiliary function for refinement
# ---------------------------------

ChiSG <- function(x, k=1.5477) {
  z <- x/k
  (3*z^2 - 3*z^4 +z^6)*(abs(x) <= k) + 1*(abs(x) > k)
}

ChiSMM <- function(x, k=3.444) {
  z <- x/k
  (3*z^2 - 3*z^4 +z^6)*(abs(x) <= k) + 1*(abs(x) > k)
}

PsiSG  <- function(x, k=1.5477) {
  z <- x/k
  (6*z-12*z^3+6*z^5)/k*(abs(x) <= k) + 0*(abs(x) > k)
}

kpmr <- function(y, d, eps=1.e-5) { 
# right-continuous Kaplan Meier estimate
# y: duration data
# d: censoring indicator (1=full, 0=censored)
# Attention: returns NA if the largest observation is censored
  oy <- order(y)
  do <- d[oy]
  yi <- yo <- y[oy]
  ny <- length(y)
  yi <- yi[do==1]
  e <- n <- rep(0, length(yi))
  e[1] <- sum(do[abs(yo - yi[1]) <= eps])
  n[1] <- ny - 1
  nu <- 1
  for (i in 2:ny) {
    if (abs(yo[i]-yi[nu]) <= eps | do[i]==0)
      next
    nu <- nu + 1
    yi[nu] <- yo[i]
    e[nu] <- sum(do[abs(yo - yi[nu]) <= eps])
    n[nu] <- ny - sum(yo < yi[nu])
  }
  yi <- yi[1L:nu]
  e <- e[1L:nu]
  n <- n[1L:nu]
  Shat <- cumprod(1-e/n)
  Fhat <- 1-Shat
  pt   <- as.numeric(diff(c(0, Fhat)))
  ans <- list(yi=yi, Shat=Shat, Fhat=Fhat, pt=pt)
  return(ans)
}

F.KM <- function(x, tu, pt) { # KM cdf at x
  sum(pt[tu <= x])
}

## SPsiKMR <- function(r, sigma, tu, pt, k=1.5477) {
##   Psiu <- PsiSG(tu/sigma, k=k)
##   isum <- tu > r
##   s <- sum(Psiu[isum]*pt[isum])
##   den  <- 1-F.KM(r,tu,pt)
##   s/den
## }

## vi era un errore nella funzione in robustAFTtemp che e' stato qui sistemato
SPsiKM <- function(r, sigma, tu, pt, k=1.5477) {
  nr <- length(r)
  nt <- length(tu)
  res <- .Fortran("SPsiKM",
     as.double(r/sigma),
     as.integer(nr),
     as.double(tu/sigma),
     as.double(pt),
     as.integer(nt),
     as.double(k),
     s=double(nr),
     PACKAGE="robustloggamma"
   )
   return(res$s)
}

SChiKM <- function(r, sigma, tu, pt, k=1.5477) {
  Chiu <- ChiSG(tu/sigma, k=k)
  isum <- tu > r
  s <- sum(Chiu[isum]*pt[isum])
  den  <- 1-F.KM(r,tu,pt)
  s/den
}

Refeqn1NP <- function(Beta, sigma, X, y, delta, xk=1.5477) {
  n <- length(y)
  p <- ncol(X)
  hu <- hc <- rep(0,p)
  r  <- as.vector(y-X%*%Beta)
  delta[(r==max(r))] <- 1
  nu <- sum(delta)  
  if (nu > 0) {
    Xu <- X[delta==1,]
    ru <- r[delta==1]
    hu <- t(Xu)%*%as.matrix(PsiSG(ru/sigma))
  }
  if (nu < n) {
    Xc <- X[delta==0,]
    rc <- r[delta==0]
    km <- kpmr(r,delta)
    hc <- t(Xc)%*%SPsiKM(rc,sigma=sigma,tu=km$yi,pt=km$pt,k=xk)
  }
  (hu+hc)/n
}

Refeqn2NP <- function(sigma, Beta, X, y, delta, xk=1.5477) {
  n <- length(y)
  p <- ncol(X)
  ku <- kc <- 0
  r  <- as.vector(y-X%*%Beta)
  delta[(r==max(r))] <- 1
  nu <- sum(delta)   
  if (nu > 0) {
    ru <- r[delta==1]
    ku <- sum(ChiSG(ru/sigma))
  }
  if (nu < n) { 
    rc <- r[delta==0]
    km <- kpmr(r,delta)
    kc <- sum(apply(as.matrix(rc),1,SChiKM,sigma=sigma,tu=km$yi,pt=km$pt))
  }
  (ku+kc)/(n-p)-0.5
}

Refopt1NP <- function(Beta, sigma, X, y, delta) {
  z <- Refeqn1NP(Beta,sigma,X,y,delta)
  as.numeric(Nrm2(z))
}

# MM optimization
# ---------------

MChi <- function(rG, sigma, tu, pt, k=3.444) {
  r <- rG[1]
  G <- rG[2]
  Chiu <- ChiSMM((tu-G)/sigma, k=k)
  isum <- tu >= r
  s <- sum(Chiu[isum]*pt[isum])
  den  <- 1-F.KM(r,tu,pt)
  s/den
}

MMobj <- function(Gam,sigma,r0,X,delta, xk=3.444) {
  n  <- length(r0)
  p <- ncol(X)
##  dfcomn(ipsi=4,xk=xk)
  ku <- kc <- 0
  Gx <- as.vector(X%*%Gam)
  delta[r0 == max(r0)] <- 1 
  nu <- sum(delta)   
  if (nu > 0) {
    ru <- r0[delta==1]
    Gu <- Gx[delta==1]
    ku <- sum(ChiSMM((ru-Gu)/sigma), k=xk)
  }
  if (nu < n) { 
    rc <- r0[delta==0]
    Gc <- Gx[delta==0]
    km <- kpmr(r0,delta)
    rcGc <- cbind(rc,Gc)
    kc <- sum(apply(as.matrix(rcGc),1,MChi,sigma=sigma,tu=km$yi,pt=km$pt))
  }
  (ku+kc)/n
}
