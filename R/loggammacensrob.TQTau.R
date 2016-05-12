#############################################################
#	loggammacensrob.TQTau function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 16, 2015
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################
####  TWQtauCens

loggammacensrob.TQTau <- function(x, delta, w=rep(1, length(x)), control=loggammarob.control()) {
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
          tol=control$refine.tol,Qrefine=FALSE)
  result <- list(mu=QT$mu1, sigma=QT$sig1, lambda=QT$lam1, weights=w, iterations=NULL, error=NULL)  
  return(result)
}
