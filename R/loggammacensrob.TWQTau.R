#############################################################
#	loggammacensrob.TWQTau function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 16, 2015
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################
####  TWQtauCens

loggammacensrob.TWQTau <- function(x, delta, w=rep(1, length(x)), control=loggammarob.control()) {
  lgrid <- seq(control$lower, control$upper, length.out=control$n)
  n <- length(x)
  # TWQtau for lambda, mu, sigma based on residuals w.r.t. an initial regression.
  SQ <- SurvQuant(x,delta,w) # Kaplan Meier
  ri <- SQ$ri # ri=dati non censurati ordinati
  qi <- SQ$qi # qi=quantili di kaplan meier
  pos <- which(qi <= control$qthreshold)
  ri <- ri[pos]
  qi <- qi[pos]
  QT <- WQtauCens(ri=ri,qi=qi,w=rep(1,length(ri)),lgrid=lgrid,
          c1=control$tuning.rho,c2=control$tuning.psi,
          N=control$nResample, maxit=control$maxit,
          tol=control$refine.tol,Qrefine=FALSE) # TQtau
  muTQT <- QT$mu1
  sigTQT <- QT$sig1
  lamTQT <- QT$lam1
  vg0 <- summary(survfit(Surv(x,delta)~1,weights=w))$std.err^2 # Varianze di Greenwood
  vg <- vg0
  if (is.na(vg[length(vg)])) {
    vg[length(vg)] = vg0[length(vg0)-1]
  } 
  vg <- vg[pos]
  ql <- qloggamma(qi,lambda=lamTQT)
  dl <- dloggamma(ql,lambda=lamTQT)
  vg1 <- vg/(dl^2)
  wl <- 1/sqrt(vg1) # pesi per il WQtau
  WQT <- WQtauCens(ri=ri,qi=qi,w=wl,lgrid=lgrid,
           c1=control$tuning.rho,c2=control$tuning.psi,
           N=control$nResample, maxit=control$maxit,
           tol=control$refine.tol,Qrefine=FALSE) # WQtau
  result <- list(mu=WQT$mu2, sigma=WQT$sig2, lambda=WQT$lam2, weights=w, iterations=NULL, error=NULL)
  return(result)
}
