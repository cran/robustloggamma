#############################################################
#	loggammacensrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 16, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

loggammacensrob <- function(x, delta, start=NULL, weights=rep(1, length(x)), method=c("oneTML", "oneWL", "TWQTau", "TQTau", "ML"), control, ...) {
  method <- match.arg(method)
  x <- na.omit(x)
  if (!is.null(na <- attr(x, "na.action")))
    weights <- weights[-na]
  
  if (missing(control)) 
    control <- if (missing(method))
                 loggammarob.control(...)
               else
                 loggammarob.control(method = method, ...)
    cl <- match.call()

  if (!missing(control) && !missing(method) && method != control$method) {
    warning("Methods argument set by method is different from method in control\n", "Using method = ", method)
    control$method <- method
  }
  ## x.orig <- x
  ## mx <- median(x)
  ## sx <- mad(x)
  ## x <- (x - mx)/sx
  
  if (control$method=="TQTau")
    result <- loggammacensrob.TQTau(x, delta, weights, control)
  else if (method=="TWQTau")
    result <- loggammacensrob.TWQTau(x, delta, weights, control)
  else if (method=="oneWL")
    result <- loggammacensrob.oneWL(x, delta, start, weights, control)
  else if (method=="oneTML")
    result <- loggammacensrob.oneTML(x, delta, start, weights, control)
  else if (method=="ML")
    result <- loggammacensrob.ML(x, delta, start, weights, control)
  result$eta <- Exp.response(result$mu, result$sigma, result$lambda)
  result$data <- x
  result$delta <- delta
  result$method <- method
  result$call <- cl
  class(result) <- c("loggammacensrob", "loggammarob")
  return(result)
}

#############################################################
#	summary.loggammacensrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 16, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

########### SUMMARY
summary.loggammacensrob <- function(object, p=NULL, conf.level=0.95, prob=0.00001, ...) {
  if ((object$method %in% c("oneTML", "oneWL", "ML"))) {
    p <- p[ifelse(p > 0 & p < 1, TRUE, FALSE)]
    nrep <- 4+length(p)
    if (length(conf.level) == 1)
      conf.level <- rep(conf.level, nrep)  
    if (length(conf.level) != nrep)
      stop("'conf.level' must have length equal to 1 or 4+length(p)")

    z <- qnorm(1-(1-conf.level)/2)
    if (is.finite(object$eta)) {
      if ((object$method %in% c("oneWL", "ML"))) {
        asd <- sqrt(diag(cov <- VCOV.WMLCensloggamma(par=c(object$mu,object$sigma,object$lambda), y=object$data, delta=object$delta, w=object$weights)))
      } else {
        asd <- sqrt(diag(cov <- VCOV.TMLoneCensloggamma(y=object$data,delta=object$delta,w=object$weights,mu=object$mu,sigma=object$sigma,lambda=object$lambda,cl=object$cut.lower,cu=object$cut.upper)))
      }
      asdeta <- sqrt(VAR.cens.eta.loggamma(cov=cov, mu=object$mu, sigma=object$sigma, lambda=object$lambda))
    } else {
      asd <- rep(NA,3)
      asdeta <- NA
    }
    object$muconf.int <- object$mu+c(-1,1)*z[1]*asd[1]/sqrt(sum(object$weights))
    object$sigmaconf.int <- object$sigma+c(-1,1)*z[2]*asd[2]/sqrt(sum(object$weights))
    object$lambdaconf.int <- object$lambda+c(-1,1)*z[3]*asd[3]/sqrt(sum(object$weights))
    object$etaconf.int <- object$eta+c(-1,1)*z[4]*asdeta/sqrt(sum(object$weights))

    attr(object$muconf.int, "conf.level") <- conf.level[1]
    attr(object$sigmaconf.int, "conf.level") <- conf.level[2]
    attr(object$lambdaconf.int, "conf.level") <- conf.level[3]
    attr(object$etaconf.int, "conf.level") <- conf.level[4]

    object$muse <- asd[1]/sqrt(sum(object$weights))
    object$sigmase <- asd[2]/sqrt(sum(object$weights))
    object$lambdase <- asd[3]/sqrt(sum(object$weights))
    object$etase <- asdeta/sqrt(sum(object$weights))

    object$p <- p
    
    if (!is.null(p)) {
      object$q <- asdq <- rep(NA, length(p))
      object$qconf.int <- matrix(NA, nrow=length(p), ncol=2)
      for (i in 1:length(p)) {
        object$q[i] <- qloggamma(p=p[i], mu=object$mu, sigma=object$sigma, lambda=object$lambda)
        if (is.finite(object$eta))
          asdq[i] <- sqrt(VAR.cens.quantile.loggamma(p=p[i], cov=cov, mu=object$mu, sigma=object$sigma, lambda=object$lambda))
        else
          asdq[i] <- NA
        object$qconf.int[i,] <- object$q[i]+c(-1,1)*z[4+i]*asdq[i]/sqrt(sum(object$weights))
      }
#####      object$qconf.int <- drop(qconf.int)
      attr(object$qconf.int, "conf.level") <- conf.level[4+(1:length(p))]
      object$qse <- asdq/sqrt(sum(object$weights))
    }
  }

  object$call <- match.call()
  class(object) <- c("summary.loggammacensrob", "summary.loggammarob")
  return(object)
}

#############################################################
#	print.summary.loggammacensrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 16, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

#PRINT.SUMMARY
print.summary.loggammacensrob <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("Location: ", format(x$mu, digits=digits))
  if ((x$method %in% c("oneTML", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$muse, digits=digits), "\n")
    cat("(", format(x$muconf.int[1L], digits=digits), ", ",
        format(x$muconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$muconf.int, "conf.level")),
        "percent confidence interval\n\n")
  } else {
    cat("\n")
  }
  cat("Scale: ", format(x$sigma, digits=digits))
  if ((x$method %in% c("oneTML", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$sigmase, digits=digits), "\n")    
    cat("(", format(x$sigmaconf.int[1L], digits=digits), ", ",
        format(x$sigmaconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$sigmaconf.int, "conf.level")),
        "percent confidence interval\n\n")
  } else {
    cat("\n")
  }
  cat("Shape: ", format(x$lambda, digits=digits))
  if ((x$method %in% c("oneTML", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$lambdase, digits=digits), "\n")    
    cat("(", format(x$lambdaconf.int[1L], digits=digits), ", ",
        format(x$lambdaconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$lambdaconf.int, "conf.level")),
        "percent confidence interval\n\n")
  } else {
    cat("\n")
  }
  cat("Mean(exp(X)): ", format(x$eta, digits=digits))
  if ((x$method %in% c("oneTML", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$etase, digits=digits), "\n")    
    cat("(", format(x$etaconf.int[1L], digits=digits), ", ",
        format(x$etaconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$etaconf.int, "conf.level")),
        "percent confidence interval\n\n\n")
  } else {
    cat("\n")
  }

######### Quantiles estimation and confidence interval
  if (!is.null(x$p)) {
    for (i in 1:length(x$p)) {
      cat("Quantile of order ", format(x$p[i], digits=digits), ": ", format(x$q[i], digits=digits))
      if ((x$method %in% c("oneTML", "oneWL", "ML"))) {
        cat(" s.e. ", format(x$qse[i], digits=digits), "\n")    
        cat("(", format(x$qconf.int[i,1L], digits=digits), ", ",
        format(x$qconf.int[i,2L], digits=digits), ") \n")
        cat(format(100 * attr(x$qconf.int, "conf.level")[i]),
          "percent confidence interval\n\n\n")
      } else {
        cat("\n")
      }
    }
  }

  if ((x$method %in% c("oneTML", "oneWL", "ML")))
    cat("\n\n")
  else
    cat("\n\nConfidence intervals and Standard Errors are available only for 'oneTML', 'oneWL' and 'ML'\n\n\n") 
  
  if (x$method %in% c("oneTML", "oneWL", "TQTau", "TWQTau"))
    summarizeRobWeights(x$weights, digits = digits, ...)
  invisible(x)
}

#############################################################
#	VAR.cens.eta.loggamma function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: December, 04, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

VAR.cens.eta.loggamma <- function(cov, mu=0, sigma=1, lambda=0, eps=0.0001, npoints=100000) {
  gradient <- c(Exp.response1(mu, sigma, lambda, eps, npoints), Exp.response2(mu, sigma, lambda, eps, npoints), Exp.response3(mu, sigma, lambda, eps, npoints))
  res <- drop(t(gradient)%*%cov%*%gradient)
  return(res)
}

#############################################################
#	VAR.cens.quantile.loggamma function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: December, 04, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

VAR.cens.quantile.loggamma <- function(p, cov, mu=0, sigma=1, lambda=0) {
## p: quantile order  
  gradient <- c(d.Q1(p, mu, sigma, lambda), d.Q2(p, mu, sigma, lambda), d.Q3(p, lambda))
  res <- drop(t(gradient)%*%cov%*%gradient)
  return(res)
}
