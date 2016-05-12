TMLCensloggammaCutoff.control <- function(pcut=0.997,gridsize=2000,xmax=NULL,step=1) {
  control <- list(pcut=pcut,gridsize=gridsize,xmax=xmax,step=step)
  return(control)
}

TMLCensloggammaCutoff <- function(y,delta,mu0,bet0,sig0,lam0,control) {
  rs0 <- (y-mu0)/sig0
  if (is.null(control$xmax)) {
    control$xmax <- 0
    while (ploggamma(control$xmax,lambda=lam0) < 1) {
      control$xmax <- control$xmax+0.01
    }
  }
  rmax <- max(control$xmax,max(rs0))+0.1
  cu <- tutl.loggamma(control$pcut,lam0)
  grid <- seq(from=cu,to=rmax,length=control$gridsize)
  zap <- adapt.CutoffLG(rs0,delta,lam0,grid)
  cl <- zap$xl
  cu <- zap$xu+0.001
  ans <- c(cl, cu)
  return(ans)
}

# Functions to compute the truncation points
# ------------------------------------------

tutl.loggamma <- function(p,lambda) {
# p-quantile of standard distribution
  if (p>1 | p<0) {
    warning("Argument out of range\n")
    return(0)
  }
  qloggamma(p=p,lambda=lambda)
}

rho.loggamma <- function(uu,lambda,zero,xmax=1e100){
  if (missing(zero))
    zero <- 0.0001
# negative loglikelihood
  if (abs(lambda) > zero) {
    alpha <- 1/lambda^2 
    llam  <- log(abs(lambda))
    lamu  <- lambda*uu
    ka    <- llam-2*alpha*llam-lgamma(alpha)
    res   <- -ka-alpha*(lamu-exp(lamu))
  } else {
    res <- log(sqrt(2*pi))+(uu^2)/2
  }
  res <- min(res, xmax)
  return(res)
} 

F0.loggamma <- function(u,lambda,zero) {
# Fo+(u) : model distribution of the negative loglikelihood
  if (missing(zero))
    zero <- 0.0001
  if (abs(lambda) > zero) {
    alpha <- 1/lambda^2 
    llam  <- log(abs(lambda))
    lamu  <- lambda*u
    ka    <- llam-2*alpha*llam-lgamma(alpha)
    if (u <= alpha-ka)
      return(0)
  }
  tt <- solverho.loggamma.eqn(lambda,u)
  Tu <- tt[2]
  Tl <- tt[1]
  z  <- ploggamma(q=Tu,lambda=lambda)-ploggamma(q=Tl,lambda=lambda)
  return(z)
}

solverho.loggamma.eqn <- function(lambda,k,pos=TRUE,neg=TRUE,tol=1e-9,maxit=300,zero) {
  if (missing(zero))
    zero <- 0.0001
# to solve negative log-lik = k
# note: if lambda < -tol, negative and positive solutions are exchanged
  lam <- lambda
  nitup <- nitun <- NA
  if (abs(lam) > zero) {
    alpha <- 1/lam^2 
    llam  <- log(abs(lam))
    ka    <- llam-2*alpha*llam-lgamma(alpha)
    un <- up <- nitun <- nitup <- NA
    if (is.na(k)  | is.na(alpha) | is.na(ka) )
      return(c(NA,NA,NA,NA))
    if (k  < alpha-ka) {
      warning(paste("No solution"))
      return(c(NA,NA,NA,NA))
    }
    if (k == alpha-ka)
      return(c(0,0,NA,NA))
# positive solution
    if (pos) {
      nit  <- 0
      conv <- FALSE
      up   <- log((k+ka)/alpha)/lam
      while(!conv){
        lamup <- lam*up
        f0 <- -ka - alpha*(lamup-exp(lamup)) - k
        f1 <- ( exp(lamup) - 1 )/lam
        delta <- -f0/f1
        conv  <- abs(f0) <= tol
        up    <- up+delta
        ## cat(up, delta, f0, f1, "\n")
        ## if (is.nan(f0)) browser()
        nit   <- nit+1
        if (nit==maxit | delta==0)
          break
      }
      nitup <- nit
    }
# negative solution
    if (neg) {
      nit  <- 0; conv <- F
      un   <- - lam*(k+ka)
      while(!conv){
        lamun <- lam*un
        f0 <- -ka - alpha*(lamun-exp(lamun)) - k
        f1 <-  ( exp(lamun) - 1 )/lam
        delta <- - f0/f1
        conv  <- abs(f0) <= tol
        un    <- un+delta
        nit   <- nit+1
        if (nit==maxit)
          break
      }
      nitun <- nit
    }
  } else {
    d <- 2*(k-log(sqrt(2*pi)))
    if (is.na(d))
      return(c(NA,NA,NA,NA))
    if (d > 0) { 
      absu  <- sqrt(d)
      un    <- -absu
      up    <-  absu
      nitun <- nitup <- 1
    } else {
      un <- NA
      up <- NA
    }
  }
  sol <- c(un,up)
  if (lambda < -zero)
    sol <- c(up,un)
  ans <- c(sol,nitun,nitup)
  return(ans)
}

# censoring case
# --------------

adapt.CutoffLG <- function(res,delta,lambda,grid) {
  n <- length(res)
  nu <- sum(delta)
  nc <- n-nu
  ngrid <- length(grid)
  rU <- res[delta==1]
  rC <- res[delta==0]
#
# calcolo della cdf e del rapporto
#
  GU   <- GC <- GN <- RP <- rep(0,ngrid)
  for (i in 1:ngrid) {
    x <- grid[i]
    rhox  <- rho.loggamma(x,lambda)
####    cat(rhox, "\n")
    if (is.infinite(rhox)) browser()
    kx    <- solverho.loggamma.eqn(lambda,rhox)[1]
    if ( nu>0)
      GU[i] <- sum( rU <= x & rU > kx )/n
    if ( nc>0) {
      num   <- (ploggamma(q=x,lambda=lambda) - ploggamma(q=pmax(rC,kx), lambda=lambda))[rC < x]
      den   <- (1 - ploggamma(q=rC,lambda=lambda))[rC < x]
      tmp   <- num[den>0]/den[den>0]
      GC[i] <- sum(tmp)/n
    }
    GN[i] <- GU[i]+GC[i]
    RP[i] <- GN[i]/(ploggamma(x,lambda=lambda) - ploggamma(kx,lambda=lambda))
  }
#
# calcolo di alpha e dei cutoff
#
  alpha <- min(min(RP),1)
  xu    <- max(grid[GN < alpha] )
  rhox  <- rho.loggamma(xu,lambda)
  xl    <- solverho.loggamma.eqn(lambda,rhox)[1]
  ans <- list(alpha=alpha,xl=xl,xu=xu,GN=GN,RP=RP)
  return(ans)
}

#---------

fast.adapt.CutoffLG <- function(res,delta,lambda,cu,xmax=NULL) { # ancora imperfetto
  if (is.null(xmax)) {
    xmax <- 0
    while (ploggamma(xmax,lambda=lambda) < 1) {
      xmax <- xmax+0.01
    }
  }
  n <- length(res) 
  nu <- sum(delta); nc <- n-nu
  rU <- res[delta==1]; rC <- res[delta==0]
#
# Punti di salto di Gn e del rapporto
#
  rho.rU <- IrU <- rep(0,nu)
  for (j in 1:nu) {
    rho.rU[j] <- rho.loggamma(rU[j], lambda)
    IrU[j]    <-  solverho.loggamma.eqn(lambda,rho.rU[j])[2] }
#
# cdf e rapporto sui valori corrispondenti ai salti di del rapporto
#
    Irr <- c(IrU[IrU>cu]-1e-100, xmax/3);  nrr <- length(Irr)
Gu <- Gc <- Gn <- Rp <- rep(0,nrr)
for (i in 1:nrr) {
  x <- Irr[i]
  rhox  <- rho.loggamma(x,lambda)
  kx    <- solverho.loggamma.eqn(lambda,rhox)[1]
if( nu>0)  Gu[i] <- sum( rU <= x & rU > kx )/n
if( nc>0) {
  num   <- (ploggamma(q=x,lambda=lambda) - ploggamma(q=pmax(rC,kx),lambda=lambda))[rC < x]
  den   <- (1-ploggamma(q=rC,lambda=lambda))[rC < x]
  tmp   <- num[den>0]/den[den>0]
  Gc[i] <- sum(tmp)/n}
Gn[i] <- Gu[i]+Gc[i]
Rp[i] <- Gn[i]/(ploggamma(q=x,lambda=lambda) - ploggamma(q=kx,lambda=lambda)) }
#
ord <- order(Irr)
Irr <- Irr[ord]
Gn  <- Gn[ord]
Rp  <- Rp[ord]
#
alpha1 <- min(min(Rp),1)                                 # non sempre coincide con la ricerca sul grid
xu1    <- max(Irr[ Gn <  alpha1 ])                       # spesso non coincide con la ricerca sul grid
rhox1  <- rho.loggamma(xu1,lambda)
xl1    <- solverho.loggamma.eqn(lambda,rhox1)[1]
res <- list(alpha=alpha1,xl=xl1,xu=xu1,Gn=Gn,Rp=Rp,Irr=Irr) 
res}

#--- Old functions --------------------------------------------------------------------------

RappLG  <- function(x,res,delta,lambda){ 
# computes semi-empirical cdf of rho_lambda(u) at rho_lambda(x) for positive x
# and ratio to compute cutoff
if (x <=0) cat("negative argument in RappLG", "\n")
n <- length(delta); FU <- FC <- 0; nu <- sum(delta); nc <- n-nu
xmax=0; while (ploggamma(q=xmax,lambda=lambda)< 1) {xmax=xmax+0.01}
rhox   <- rho.loggamma(x,lambda)
sol    <- izerox <- solverho.loggamma.eqn(lambda,rhox)
#
izerox <- sol[1]
rU <- res[delta==1]
rC <- res[delta==0]
if (nu > 0)            FU <- sum( rU <= x & rU > izerox )
#
if (nc > 0 & x > xmax) FC <- nc
if (nc > 0 & x<= xmax) { 
 tmp     <- rep(0,nc); i3 <- rC < xmax
 num     <- (ploggamma(q=x,lambda=lambda) - ploggamma(q=pmax(rC,izerox),lambda=lambda)  )*(rC <= x)
 den     <- 1-ploggamma(q=rC,lambda=lambda)
 i4 <- i3 & (den > 0)
 tmp[i4] <- num[i4]/den[i4]
 tmp[is.na(tmp)] <- 1
 FC <- sum(tmp) }
#
Fnz <- (FC+FU)/n
Rap <- Fnz/(ploggamma(q=x,lambda=lambda) - ploggamma(q=izerox,lambda=lambda))
c(Fnz, Rap, FU/n, FC/n) }

adapt.cutoff.loggamma.cens <- function(res,delta,lambda,p=0.995,gridsize=2000){
# computes cutoff on residual scale 
xmax=0; while (ploggamma(q=xmax,lambda=lambda)< 1) {xmax=xmax+0.01}
xmax=xmax+0.25
cu <- tutl.loggamma(p,lambda)
xs <- seq(from=cu,to=xmax,length=gridsize)
fs <- apply(as.matrix(xs),1,RappLG,res=res,delta=delta,lambda=lambda)
alpha <- min(min(fs[2,]),1)
xu    <- max(xs[fs[1,] < alpha])
rhox  <- rho.loggamma(xu,lambda)
xl    <- solverho.loggamma.eqn(lambda,rhox)
list(alpha=alpha,cl=xl[1],cu=xu) }

#------------------------------------------------------------------------------

# no censoring case
# -----------------

Discr <- function(yo,Po,cu){
# minimum ratio between Fn and Po (for argument > cu)
n  <- length(yo)
# if (cu > yo[n]) {cat("cu is too large ","\n"); return(list(j=NA))}
Fn <- (0:(n-1))/n; r  <- Fn/Po
j  <- (1:n)[yo > cu]
if (length(j)>0) {j <- j[1]; rj <- r[j:n]} else {j <- n+1; rj <- 1}
r.min <- min(rj,1)
if (r.min < 1) j.min <- (j:n)[rj==r.min] else j.min <- NULL
list(j=j,j.min=j.min,r.min=r.min)}

adaptCutoff.loggamma <- function(rso,pu=0.95,lam0,option,maxit,zero){
  if (missing(zero))
    zero <- 0.0001
# adaptive cut-off values
Lambda <- NA; n <- length(rso)  
  cu       <- tutl.loggamma(pu,lam0)
  rhocu    <- rho.loggamma(cu,lam0)
  cl       <- solverho.loggamma.eqn(lam0,rhocu)[1]
if (option=="adaptive"){
#
# compute rhotu (cutoff on rho scale)
#
  rhoi     <- rho.loggamma(rso,lam0)                          
  rhos     <- sort(rhoi)
  Po       <- apply(as.matrix(rhos),1,F0.loggamma,lam=lam0)
  dsr      <- Discr(rhos,Po,rhocu)
  Lambda   <- dsr$r.min
  tu2      <- as.numeric(quantile(rhos,probs=Lambda))
  tmp      <- rhos[rhos < tu2]
  if (length(tmp)==0) return(list(rhocu=rhocu,rhotu=NA,cl=cl,cu=cu,tl=NA,tu=NA,Lambda=Lambda))
  rhotu    <- max(tmp)
#
# compute cutoff on data scale
#
  zz       <- solverho.loggamma.eqn(lam0,rhotu,neg=FALSE)
  if(lam0 >= -zero) {
    tu <- zz[2]
    tu <- max(tu,cu)
    rhotu <- rho.loggamma(tu,lam0)
    zz    <- solverho.loggamma.eqn(lam0,rhotu,pos=FALSE)
    tl    <- zz[1] }                                                                          
  else {
    tl <- zz[1]
    tl <- min(tl,cl)
    rhotl <- rho.loggamma(tl,lam0)
    zz    <- solverho.loggamma.eqn(lam0,rhotl,pos=FALSE)
    tu    <- zz[2]}}
if (option=="fixed"){tu <- cu; tl <- cl; rhotu <- rhocu}
list(rhocu=rhocu,rhotu=rhotu,cl=cl,cu=cu,tl=tl,tu=tu,Lambda=Lambda)}

adaptCutoff_d.loggamma <- function(rso,pu=0.95,lam0,option,maxit,zero){
  if (missing(zero))
    zero <- 0.0001
# adaptive cut-off values based on difference
Lambda <- NA; n <- length(rso)  
  cu       <- tutl.loggamma(pu,lam0)
  rhocu    <- rho.loggamma(cu,lam0)
  cl       <- solverho.loggamma.eqn(lam0,rhocu)[1]
if (option=="adaptive"){
#
# compute rhotu (cutoff on rho scale)
#
  rhoi     <- rho.loggamma(rso,lam0)                          
  rhos     <- sort(rhoi)
  Po       <- apply(as.matrix(rhos),1,F0.loggamma,lam=lam0)
  tmp      <- (1:n)[rhos < rhocu]
  if (length(tmp)==0) return(list(rhocu=rhocu,rhotu=NA,cl=cl,cu=cu,tl=NA,tu=NA,Lambda=NA))
  ic  <- max( tmp )
  Fn  <- (0:(n-1))/n
  Fnc <- Fn[(ic+1):n]
  Poc <- Po[(ic+1):n]
  dn  <- max( Poc - Fnc, 0)
  ni <- n-floor(n*dn)
  rhotu <- rhos[ni]
#
# compute cutoff on data scale
#
  zz       <- solverho.loggamma.eqn(lam0,rhotu,neg=FALSE)
  if(lam0 >= -zero) {
    tu <- zz[2]
    tu <- max(tu,cu)
    rhotu <- rho.loggamma(tu,lam0)
    zz    <- solverho.loggamma.eqn(lam0,rhotu,pos=FALSE)
    tl    <- zz[1]
  } else {
    tl <- zz[1]
    tl <- min(tl,cl)
    rhotl <- rho.loggamma(tl,lam0)
    zz    <- solverho.loggamma.eqn(lam0,rhotl,pos=FALSE)
    tu    <- zz[2]}}
if (option=="fixed"){tu <- cu; tl <- cl; rhotu <- rhocu}
list(rhocu=rhocu,rhotu=rhotu,cl=cl,cu=cu,tl=tl,tu=tu,Lambda=1-dn)
}

