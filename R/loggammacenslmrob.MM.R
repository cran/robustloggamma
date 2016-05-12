loggammacenslmrob.MM.control <- function(method="MM", Srefine=FALSE,N=150,q=8,sigma0=1,max.it=100,tol=0.001,ialg=3,seed=123,bb=0.5, tuning.rho=1.547647, tuning.psi=3.444, cint="Nonpar",tolB=0.001,tolS=0.01,maxit.Srefine=20,d=NULL,lower=NULL,upper=NULL,verbose=FALSE) {
  control <- list(method=method,Srefine=Srefine,N=N,q=q,sigma0=sigma0,max.it=max.it,tol=tol,ialg=ialg,seed=seed,bb=bb,tuning.rho=tuning.rho,tuning.psi=tuning.psi,cint=cint,tolB=tolB,tolS=tolS,maxit.Srefine=maxit.Srefine,d=d,lower=lower,upper=upper,verbose=verbose)
  return(control)
}

loggammacenslmrob.MM <- function(x,y,delta,init=NULL,control) {    
  SNP <- loggammacenslmrob.S(x=x,y=y,delta=delta,control=control)
  init <- SNP$coefficients
  scale <- SNP$scale
  if (control$Srefine) {
    RNP <- RefineSNP(X=x,y=y,delta=delta,Bmin0=init,Smin0=scale,
                       tolB=control$tolB,tolS=control$tolS,
                       maxit=control$maxit.Srefine,d=control$d,
                       lower=control$lower,upper=control$upper,
                       verbose=control$verbose)
    init <- RNP$coefficients
    scale <- RNP$scale
  } else {
    RNP <- NULL
  }
  MM <- MMest(x=x, y=y, delta=delta, init=init, scale=scale, tuning.psi=control$tuning.psi)
  fitted.values <- drop(x%*%MM$coefficients)
  residuals <- y-fitted.values
  res <- list(coefficients=MM$coefficients, scale=scale, fitted.values=fitted.values, residuals=residuals, init=init, init.S=SNP, init.Sr=RNP, init.M=MM, control=control, converged=TRUE)
  return(res)
}

RefineSNP <- function(X, y, delta, Bmin0, Smin0, tolB=0.001, tolS=0.001, maxit=10, d=NULL, lower=NULL, upper=NULL, verbose=FALSE) {
# Refinement for non parametric S, Gaussian case
  p <- ncol(X)
  Nit <- 1
  Mxt <- 10
  Bmin <- Bmin0
  Smin <- Smin0
  DeltaB <- DeltaS <- 100000
  while(( DeltaS > tolS | DeltaB > tolB ) & Nit <= maxit) {
# determine lower and upper values for sigma
    smu  <- sml  <- Smin
    nit  <- 1 
    funu <- funl <- Refeqn2NP(smu,Bmin,X,y,delta)
    while(funu > 0 & nit <= Mxt) {
      smu  <- smu*1.5
      funu <- Refeqn2NP(smu,Bmin,X,y,delta)
      nit  <- nit+1
    }
    if (verbose)
      cat("smu, funu, nit",smu,funu,nit,"\n")
    nit <- 1
    while(funl < 0 & nit <= Mxt) {
      sml  <- sml*0.5
      funl <- Refeqn2NP(sml,Bmin,X,y,delta)
      nit  <- nit+1
    }
    if (verbose)
      cat("sml, funl, nit",sml,funl,nit,"\n")
    zS <- regfal(Refeqn2NP,cc=0,lower=sml,upper=smu,nint=5,tol=0.001,maxit=5,Beta=Bmin,X=X,y=y,delta=delta)
    Smin  <- zS$solution
    nitS <- zS$nit
    if (verbose)
      cat("Smin", Smin, nitS, "\n")
    if (verbose)
      cat("optRefopt1", Refopt1NP(Bmin,Smin,X,y,delta),"\n")
# zB   <- nlregb(nres=p,start=Bmin,res=Refeqn1NP,control=nlregb.control(abs.tol=10^(-5)),sigma=Smin,X=X,y=y,delta=delta)
    if (!is.null(d)) {
      lower <- Bmin-rep(d,length(Bmin))
      upper <- Bmin+rep(d,length(Bmin))
    }
    if (is.null(lower) | is.null(upper))
      zB <- nlminb(start=Bmin,objective=Refopt1NP,sigma=Smin,X=X,y=y,delta=delta)
    else
      zB <- nlminb(start=Bmin,objective=Refopt1NP,sigma=Smin,X=X,y=y,delta=delta, lower=lower, upper=upper)  
    
    Bmin <- zB$par
    if (verbose)
      cat("optim", zB$objective," Bmin",Bmin, " convergence:", zB$convergence, "\n")
    DeltaB <- max(abs(Bmin - Bmin0))
    DeltaS <- abs(Smin-Smin0)
    if (verbose)
      cat("Nit=",Nit,"DeltaB",DeltaB,"DeltaS",DeltaS,"\n")
    Bmin0 <- Bmin
    Smin0 <- Smin
    Nit <- Nit+1
  }
  zres <- list(coefficients=Bmin,scale=Smin,iter=Nit)
  return(zres)
}

## RefineSNP <- function(x, y, delta, init, scale, tolB=0.001, tolS=0.001, maxit=10, d=NULL, lower=NULL, upper=NULL, verbose=FALSE) {
## # Refinement for non parametric S, Gaussian case
##   p <- ncol(x)
##   Nit <- 1
##   Mxt <- 10
##   init0 <- init
##   scale0 <- scale
##   DeltaB <- DeltaS <- max(tolB,tolS)+1
##   while(( DeltaS > tolS | DeltaB > tolB ) & Nit <= maxit) {
## # determine lower and upper values for sigma
##     smu  <- sml  <- scale
##     nit  <- 1 
##     funu <- funl <- Refeqn2NP(sigma=smu, Beta=init, X=x, y=y, delta=delta, xk=control$tuning.rho)
##     while(funu > 0 & nit <= Mxt) {
##       smu  <- smu*1.5
##       funu <- Refeqn2NP(smu,init,x,y,delta)
##       nit  <- nit+1
##     }
##     if (verbose)
##       cat("smu, funu, nit",smu,funu,nit,"\n")
##     nit <- 1
##     while(funl < 0 & nit <= Mxt) {
##       sml  <- sml*0.5
##       funl <- Refeqn2NP(sml,init,x,y,delta)
##       nit  <- nit+1
##     }
##     if (verbose)
##       cat("sml, funl, nit",sml,funl,nit,"\n")
##     zS <- regfal(Refeqn2NP,cc=0,lower=sml,upper=smu,nint=5,tol=0.001,maxit=5,Beta=init,X=x,y=y,delta=delta)
##     scale  <- zS$solution
##     nitS <- zS$nit
##     if (verbose)
##       cat("scale", scale, nitS, "\n")
##     if (verbose)
##       cat("optRefopt1", Refopt1NP(init,scale,x,y,delta),"\n")
##     if (!is.null(d)) {
##       lower <- init-rep(d,length(init))
##       upper <- init+rep(d,length(init))
##     }
##     if (is.null(lower) | is.null(upper))
##       zB <- nlminb(start=init,objective=Refopt1NP,sigma=scale,X=x,y=y,delta=delta)
##     else
##       zB <- nlminb(start=init,objective=Refopt1NP,sigma=scale,X=x,y=y,delta=delta, lower=lower, upper=upper)
##     init <- zB$par
##     if (verbose)
##       cat("optim", zB$objective," init",init, " convergence:", zB$convergence, "\n")
##     DeltaB <- max(abs(init - init0))
##     DeltaS <- abs(scale-scale0)
##     if (verbose)
##       cat("Nit=",Nit,"DeltaB",DeltaB,"DeltaS",DeltaS,"\n")
##     init0 <- init
##     scale0 <- scale
##     Nit <- Nit+1
##   }
##   res <- list(coefficients=init,scale=scale,iter=Nit)
##   return(res)
## }

MMest <- function(x, y, delta, init, scale, tuning.psi) {
  p <- ncol(x)
  gamma <- rep(0, p)
  r0 <- drop(y-x%*%init)
  obj.ini <- MMobj(Gam=rep(0,p),scale,r0,x,delta,xk=tuning.psi)
  dB <- nlminb(start=rep(0,p),objective=MMobj,sigma=scale,r0=r0,X=x,delta=delta, xk=tuning.psi)
  mess    <- dB$message
  obj.fin <- dB$objective
  iter     <- dB$iterations
  if (dB$obj < obj.ini)
    gamma <- dB$par
  coefficients <- init+gamma
  list(coefficients=coefficients, gamma=gamma, init.MM=obj.ini, final.MM=obj.fin,iter=iter,message=mess)
}
