loggammacenslmrob.TQTau <- function(x,y,delta,init=NULL,w,control) {
  x <- as.matrix(x)
  p <- NCOL(x)
  if (is.null(init)) {
    keep <- w >= control$minw
    MM <- loggammacenslmrob.MM(x=cbind(1,x[keep,]),y=y[keep],delta=delta[keep],control=loggammacenslmrob.MM.control())
    init <- MM$coefficients
  }
  residuals <- drop(y - x%*%init[2:(p+1)])
  res <- loggammacensrob.TQTau(x=residuals, delta=delta, w=w, control=control)
  res$coefficients <- init[1:(p+1)]
  res$coefficients[1] <- res$mu
  res$fitted.values <- drop(x%*%res$coefficients[-1] + res$mu)
  res$residuals <- y - res$fitted.values
  res$control <- control  
  res$converged <- TRUE
  return(res)
}

