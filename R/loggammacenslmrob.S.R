loggammacenslmrob.S.control <- function(method="S",N=150, q=8, sigma0=1, max.it=100, tol=0.001, ialg=3, seed=153, bb=0.5, tuning.rho=1.547647, cint="Nonpar") {
  control <- list(method=method, N=N, q=q, sigma0=sigma0, max.it=max.it, tol=tol, ialg=ialg, seed=seed, bb=bb, tuning.rho=tuning.rho, cint=cint)
  return(control)
}

loggammacenslmrob.S <- function(x, y, delta, control) {
# ialg: 3: bisection 1:fixed point 2:regula falsi
# meth: 1: S  4: sn
  x <- as.matrix(x)
  n <- length(y)
  p <- ncol(x)
# generate matrix of Beta-s
  zbet <- BtamatG(X=x, y=y, delta=delta, N=control$N, q=control$q, seed=control$seed, cint=control$cint, maxit=control$max.it, tol=control$tol, ialg=control$ialg, bb=control$bb, xk=control$tuning.rho)
  beta <- matrix(zbet$beta, ncol=p)
  beta <- beta[!duplicated(beta),,drop=FALSE]
  N <- NROW(beta)
  gamok <- 0
  mu0 <- 0
  s0 <- 1
  kappa <- 9.E9 
  smink <- 9.E9
  mes2 <- integer(4)
  for (j in 1:N) {
    betj <- matrix(beta[j,],nrow=N,ncol=p,byrow=TRUE)
    if (gamok==0) {
      gamj <- sweep(beta,2,beta[j,])
      nu <- apply(gamj,1,Nrm2)
    }
    B <- which(nu <= kappa)
    ni <- length(B)
    if (ni==0)
      next
    betgamb <- betj[B,]
    sj <- rep(control$sigma0, ni)
    Sj <- rep(Inf, N) 
    gamjb <- gamj[B,]
    betjb <- betj[B,]
    tmp  <- S.eq.Gauss(x,y,ni,delta,sj,control$sigma0,control$bb,betjb,gamjb,control$max.it,control$tol,ipsi=4,xk=control$tuning.rho,lint=0,ialg=control$ialg)
    Sj[B] <- tmp$S
    mes2 <- mes2+tmp$mes2
    omega <- min(Sj)
    kStar <- (1:N)[Sj==omega]
    numin <- nu[kStar]
    nustar <- min(numin)
    kstar <- (kStar)[numin==nustar]
# Check for equal norm
    if (abs(nustar-kappa)<1e-6 & omega>smink)
      next
    Bbar <- setdiff(1:N,B)
    truemin <- TRUE
    if (length(Bbar)>0)
      for (k in Bbar) {
        sk <- control$sigma0 
        tmp <- S.eq.Gauss(x,y,1,delta,sk,control$sigma0,control$bb,betj[k,],gamj[k,],control$max.it,control$tol,ipsi=4,xk=control$tuning.rho,lint=0,ialg=control$ialg)
        tst <- tmp$S
        mes2 <- mes2+tmp$mes2
        if (tst < omega) {
          truemin <- FALSE
          break
        }
      }
    if (truemin) {
      smink <- omega
      Jstar <- j
      Kstar <- kstar
      smin <- Sj[kstar]
      gamin <- gamj
      kappa <- nustar
    } 
  } # end for j
  cmes2 <- c(as.character(mes2)," ")
  Eqexit <- paste(c("Normal (","Iter=Maxit (","f(a)*f(b)>0 (","|f(a)-f(b)|<tl ("," "), cmes2, sep="",collapse=")  ")
  fitted.values <- drop(x%*%beta[Jstar,])
  residuals <- y - fitted.values
  control$fast.s.large.n <- Inf
  res <- list(coefficients=beta[Jstar,], scale=smin,
    fitted.values=fitted.values, residuals=residuals,
    init=beta, gama=gamin, kappa=kappa, Jstar=Jstar,
    Kstar=Kstar, message=Eqexit, control=control, converged=TRUE)
  return(res)
}
