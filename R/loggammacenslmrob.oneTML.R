loggammacenslmrob.oneTML <- function(x,y,delta,init=NULL,w,control) {
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
  
  mu <- x%*%bet0 + mu0
  rs0 <- (y-mu)/sig0
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

    TMo <- TMLone.loggamma(IX=cbind(1,x), y=y, delta=delta, Beta.t=c(mu0,bet0), sigma.t=sig0, lambda.t=lam0, cl=cl, cu=cu, step=control$step, reparam=FALSE)
    if (as.numeric(TMo$sigma) < 0 ) {
      TMo <- TMLone.loggamma(IX=cbind(1,x), y=y, delta=delta, Beta.t=c(mu0,bet0), sigma.t=sig0, lambda.t=lam0, cl=cl, cu=cu, step=control$step, reparam=TRUE)
    }
    dif <- max(abs(c(mu0-drop(TMo$Beta[1]), bet0-TMo$Beta[-1], sig0-TMo$sigma, lam0-TMo$lambda)))
    ## if (control$verbose)
    ##   cat("iteration: ", iter, " convergence tolerance: ", dif, "\n")
    mu0 <- drop(TMo$Beta[1])
    bet0 <- TMo$Beta[-1]
    sig0 <- TMo$sigma
    lam0 <- TMo$lambda
  }
  res <- list()
  res$coefficients <- TMo$Beta
  res$mu <- TMo$Beta[1]
  res$sigma <- TMo$sigma
  res$lambda <- TMo$lambda
  res$fitted.values <- drop(x%*%res$coefficients[-1] + res$mu)
  res$residuals <- y - res$fitted.values
  res$cut.lower <- cl
  res$cut.upper <- cu
  res$iter <- iter
  res$weights <- TMo$weights
  res$errors <- TMo$err
  res$n.ret <- TMo$n.ret
  res$control <- control
  res$converged <- TRUE
  return(res)
}

VCOV.TMLoneCensloggammareg <- function(y,delta,IX,w,coefficients,sigma,lambda,cl,cu) {  
  p <- length(coefficients)
  n <- length(y)
  mu <- drop(IX %*% as.matrix(coefficients))
  rs0 <- drop((y - mu)/sigma)
  pos <- w > 0.001
  IX <- IX[pos, , drop = FALSE]
  y <- y[pos]
  delta <- delta[pos]
  w <- w[pos]
  rs0 <- rs0[pos]
  kk <- -cl*dloggamma(cl,mu=0,sigma=1,lambda=lambda) + cu*dloggamma(cu,mu=0,sigma=1,lambda=lambda)-ploggamma(cu,mu=0,sigma=1,lambda=lambda)+ploggamma(cl,mu=0,sigma=1,lambda=lambda)
  gg <- d.FL.1(cl,lambda)-d.FL.1(cu,lambda)
  JAC <- JAC.TML.loggamma(IX, y, delta, coefficients, coefficients, sigma, sigma, lambda, rs0, w, cl, cu, kk, gg)
  Q <- QMatrix.LG(rep(0,p),1,lambda,rs0,delta,w,IX,cl,cu,kk,gg) 
  JIC <- solve(JAC)
  Cov <- t(JIC) %*% Q %*% JIC /n
  return(Cov)
}
