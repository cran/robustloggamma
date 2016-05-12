FnsemiGLG <- function(x, y, delta, mu, sigma, lambda) {
# Semiparametric cdf using GLG
  n <- length(y)
  nc <- sum(delta==0 & y <= x)
  res <- .Fortran("FnsemGLG",
    as.double(x),
    as.double(y),
    as.integer(delta),
    as.integer(n),
    as.integer(nc),
    as.double(mu),
    as.double(sigma),
    as.double(lambda),
    p = double(1),
    PACKAGE="robustloggamma"
  )
  return(res$p)
}
