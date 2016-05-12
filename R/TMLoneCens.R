TMLCensloggammaWeights.control <- function(pcut=0.997,gridsize=2000,xmax=NULL,step=1) {
  control <- list(pcut=pcut,gridsize=gridsize,xmax=xmax,step=step)
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
  grid <- seq(from=cu,to=rmax,length=control$gridsize)
  zap <- adapt.CutoffLG(rs0,delta,lam0,grid)
  cl <- zap$xl
  cu <- zap$xu+0.001
  wgt <- weightsHD(rs0, cl, cu)
  return(wgt)
}

