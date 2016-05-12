QnsemiGLG <- function(x, xmin, xmax, p, y, delta, mu, sigma, lambda, tol=10^(-8)) {
  n <- length(y)
  res <- .Fortran("QnsemGLG",
    as.double(x),
    as.double(xmin),
    as.double(xmax),
    as.double(p),                  
    as.double(y),
    as.integer(delta),
    as.integer(n),
    as.double(mu),
    as.double(sigma),
    as.double(lambda),
    as.double(tol),                  
    xout = double(1),
    PACKAGE="robustloggamma"
  )
  return(res$xout)
}
