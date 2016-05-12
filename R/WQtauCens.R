# WQtauCens.R
WQtauCens <- function(ri,qi,w,lgrid,c1=1.547647,c2=6.08,N=100,maxit=750,tol=1e-6,Qrefine=TRUE) {
# Weighted Qtau estimate (grid optimization)
  m <- length(ri)
  X <- matrix(0,ncol=2,nrow=m)
  nl <- length(lgrid)
  sL1 <- sL2 <- aL1 <- aL2 <- bL1 <- bL2 <- lgrid
  qq <- list(beta = 0, Tscale = 0, nit=0)
  b1 <- integrate(rhoBWdnorm,lower=-10,upper=10,k=c1)$value
  b2 <- integrate(rhoBWdnorm,lower=-10,upper=10,k=c2)$value
  for (i in 1:nl) {
    lam <- lgrid[i]
    ql <- qloggamma(qi,lambda=lam)
    z <- RegtauW.f(x=ql,y=ri,w=w,
           b1=b1,c1=c1,b2=b2,c2=c2,N=N,tol=tol,seed=567) # Fortran
    # z  <- RegtauW(ql,ri,w,b1,c1,b2,c2,N) # S
    X <- cbind(1,ql)
    B0 <- c(z$ao,z$bo)
    s0 <- z$to
    if (Qrefine) {
      qq <- IRLStauW(X=X,y=ri,w=w,inib=B0,iniscale=s0,
              b1=b1,c1=c1,b2=b2,c2=c2,maxit=maxit,tol=tol)
    } else {
      qq$Tscale <- s0
      qq$beta <- B0
    } 
#cat(i,lgrid[i],z$to,qq$Tscale,qq$nit,"\n")
    sL1[i] <- z$to; sL2[i] <- qq$Tscale
    aL1[i] <- z$ao; aL2[i] <- qq$beta[1]
    bL1[i] <- z$bo; bL2[i] <- qq$beta[2]
  }  
  io1 <- (1:nl)[sL1 == min(sL1)]
  io1 <- min(io1)
  lam.est1 <- lgrid[io1]
  s.est1   <- sL1[io1]
  a.est1   <- aL1[io1]
  b.est1   <- bL1[io1]
  io2 <- (1:nl)[sL2 == min(sL2)]
  io2 <- min(io2)
  lam.est2 <- lgrid[io2]
  s.est2   <- sL2[io2]
  a.est2   <- aL2[io2]
  b.est2   <- bL2[io2]
  res <- list(lam1=lam.est1,Stau1=s.est1,mu1=a.est1,sig1=b.est1,
    lam2=lam.est2,Stau2=s.est2,mu2=a.est2,sig2=b.est2,sL1=sL1,sL2=sL2)
  return(res)
}

#######################
# auxiliary functions #
#######################

KpMrG <- function(t,x) {
# Kaplan-Meier with Greenwood's estimate
# t: sorted survival times
# x: censored times
if (length(x)==0) {
 k  <- length(t)
 tu <- unique(t)
 pt <- rep(1/k,k)
 pt[k+1] <- 0
 Ft <- cumsum(pt[1:k]) 
 Vt <- Ft*(1-Ft)/k
 if (is.na(Vt[k])) Vt[k] =Vt[k-1]}
else {
 n0 <- length(t)+length(x)
 n1 <- n0-sum(x<min(t))
 tu <- unique(t)
 d  <- table(t)
 k  <- length(tu)
 m  <- rep(0,k-1)
 for (i in 1:(k-1)) m[i] <- sum(tu[i] <= x & x < tu[i+1])
 nj <- rep(n1,k)
 ds <- cumsum(c(0,d)[1:k])
 nj <- nj-ds
 cs <- cumsum(c(0,m)[1:k])
 nj <- nj-cs
 Ft <- cumprod((nj-d)/nj)
 pt <- as.numeric(diff(c(0,1-Ft,1))) 
# Greenwood's estimate                                       # da controllare
 gt <- rep(0,k)
 gt <- (d/(nj*(nj-d)))
 Vt <- (Ft^2)*cumsum( gt )
 if (is.na(Vt[k])) Vt[k] =Vt[k-1] }                          # ???
# k        : number of different survival times
# tu[1:k]  : unique survival times
# Ft[1:k]  : Kaplan Meier survival function at tu (length=k)
# pt[1:k+1]: point masses at tu; 
#            the support of pt[k+1] is undefined
# Vt[1:k]  : variance of survival function
list(k=k,tu=tu,Ft=Ft,pt=pt,Vt=Vt)}

F.KM <- function(x,tu,wkm){ # KM cdf at x
sum(wkm[tu <= x])}

# Biweight functions

rhoBW <- function(x,k){
k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2; x6 <- x4*x2
(3*x2/k2-3*x4/k4+x6/k6)*(abs(x)<k)+(abs(x)>=k)}

psiBW <- function(x,k){
(6/k)*(x/k)*(1-(x/k)^2)^2*(abs(x)<k)}

pspBW <- function(x,k){
k2  <- k*k; k4 <- k2*k2; k6 <- k2*k4
x2  <- x*x; x4 <- x2*x2
(6/k2-36*x2/k4+30*x4/k6)*(abs(x)<k)}

## rhoBW <- function(x,k){
##   Mchi(x,k,psi="optimal")
## }

## psiBW <- function(x,k){
##   Mpsi(x,k,psi="optimal")
## }

## pspBW <- function(x,k){
## NA
## }

rhoBWdnorm <- function(x,k) {
  rhoBW(x,k)*dnorm(x)
}

MscaleW <- function(u,w,b1,c1,tol) {
  h <- 1
  it <- 0
  s0 <- median(abs(u*w))/.6745
  if (s0 > tol) {
    while((h > tol) & (it < 50)) {
      it <- it+1
      s1 <- (s0^2)*mean(rhoBW((u*w/s0),c1)) / b1
      s1 <- s1^(1/2)
      h  <- abs(s1-s0)/s0
      s0 <- s1
    }
  }
  return(s0)
}

TauscaleW <- function(u,w,b1,c1,b2,c2,tol) {
  tau <- tol
  s0  <- MscaleW(u,w,b1,c1,tol)
  if (s0 > tol)
    tau <- sqrt(s0^2*mean(rhoBW(u*w/s0,c2)) / b2)
  return(tau)
}
