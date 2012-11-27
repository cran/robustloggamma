### THIS FUNCTION IS NOT RELEASED YET TO THE END USER!

wWilksTest <- function(x, thetainit=NULL, ...) {
  y <- x$data
  if (is.null(x$weights)) 
    weights <- rep(1, length(y))
  else
    weights <- x$weights
  
  minusloglikH0 <- function(theta, x, weights) {
    if (is.null(dim(theta)))
      res <- drop(log(dloggamma(x, theta[1], theta[2], theta[2]))%*%weights)
    else {
      res <- rep(NA, NROW(theta))    
      for (i in 1:NROW(theta)) {
        res[i] <- drop(log(dloggamma(x, theta[i,1], theta[i,2], theta[i,2]))%*%weights)
      }
    }
    res <- -res
    return(res)
  }
  if (is.null(thetainit))
    thetainit <- c(x$mu, x$sigma)

  theta0 <- optim(par=thetainit, fn=minusloglikH0, x=y, weights=weights, ...)$par
  
  statistic <- drop(2*(log(dloggamma(y, x$mu, x$sigma, x$lambda))-log(dloggamma(y, theta0[1], theta0[2], theta0[2])))%*%weights)
  result <- list()
  result$statistic <- statistic
  result$pvalue <- pchisq(statistic, df=1, lower.tail=FALSE)
  result$mu0 <- theta0[1]
  result$sigma0 <- theta0[2]
  result$shape0 <- theta0[2]  
  return(result)
}
