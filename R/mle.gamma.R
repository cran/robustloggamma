mle.gamma <- function(x, tol=1e-8) {
  n <- length(x)
  shape1 <- mean(x)^2/(mean(x^2)-mean(x)^2)
  sumlogx <- sum(log(x))
  logmeanx <- log(mean(x))

  g <- function(x) sumlogx + n*(log(x) - logmeanx) - n*digamma(x)
  gprime <- function(x) n/x - n*trigamma(x) 
  nr <- function(x) x - g(x)/gprime(x)
  diffshape <-  tol + 1
  while (abs(diffshape) > tol) {
    shape2 <- nr(shape1)
    diffshape <- shape2 - shape1
    shape1 <- shape2
  }
  rate <- shape1/mean(x)
  ans <- list(shape=shape2, rate=rate, scale=1/rate)
  return(ans)
}
