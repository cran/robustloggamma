#############################################################
#	loggammacenslmrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 27, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

loggammacenslmrob <- function(formula, delta, data, subset, weights, na.action, method=c("oneTML", "oneWL", "TWQTau", "TQTau", "ML"), model=TRUE, x=!control$compute.rd, y=FALSE, singular.ok=TRUE, contrasts=NULL, offset=NULL, control=NULL, init=NULL, ...) {

  method <- match.arg(method)  
  chk.s <- function (...) {
    if (length(list(...))) 
      warning("arguments  ", sub(")$", "", sub("^list\\(", "", deparse(list(...), control = c()))), "  are disregarded in\n ", deparse(sys.call(-1), control = c()), call. = FALSE)
  }
  save.d <- delta
  if (miss.ctrl <- missing(control)) 
    control <- if (missing(method)) 
        loggammarob.control(...)
     else loggammarob.control(method = method, ...)
  else chk.s(...)
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset) && length(offset) != NROW(y)) 
    stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
      length(offset), NROW(y)), domain = NA)
  if (!miss.ctrl && !missing(method) && method != control$method) {
    warning("Methods argument set by method is different from method in control\n", "Using the former, method = ", method)
     control$method <- method
  }
  if (is.empty.model(mt)) {
    x <- NULL
    singular.fit <- FALSE
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
      3) else numeric(0), residuals = y, scale = NA,
      mu = NA, sigma = NA, lambda = NA, 
      fitted.values = 0 * y, cov = matrix(, 0, 0), weights = w,
      rweights = rep.int(as.numeric(NA), NROW(y)),              
      rank = 0, df.residual = NROW(y), 
      converged = TRUE, iter = 0)
      if (!is.null(offset)) {
        z$fitted.values <- offset
        z$residuals <- y - offset
        z$offset <- offset
      }
  } else {
    x <- model.matrix(mt, mf, contrasts)
    x <- x[ ,apply(x, 2, var)!=0, drop=FALSE]
    contrasts <- attr(x, "contrasts")
    assign <- attr(x, "assign")
    p <- ncol(x)
    if (!is.null(offset)) 
      y <- y - offset
    if (NCOL(y) > 1L)
      stop("y must be a vector")
    y <- as.vector(y)
    if (!is.null(w)) {
      n <- nrow(x)
      if (NROW(y) != n | length(w) != n) 
        stop("incompatible dimensions")
      if (any(w < 0 | is.na(w))) 
        stop("missing or negative weights not allowed")
      zero.weights <- any(w == 0)
      if (zero.weights) {
        save.r <- y
        save.w <- w
        save.f <- y
        save.d <- delta
        ok <- w != 0
        nok <- !ok
        w <- w[ok]
        x0 <- x[nok, , drop = FALSE]
        x <- x[ok, , drop = FALSE]
        n <- nrow(x)
        y0 <- y[nok]
        y <- y[ok]
        delta <- delta[ok]
        attr(mf, "zero.weights") <- which(nok)
      }
      wts <- sqrt(w)
      save.y <- y
    } else {
      wts <- rep(1,length(y))
    }
    ## if (getRversion() >= "3.1.0") {
    ##   z0 <- .lm.fit(wts*x, wts*y, tol = control$solve.tol)
    ##   piv <- z0$pivot
    ## } else {
      z0 <- lm.fit(wts*x, wts*y, tol = control$solve.tol)
      piv <- z0$qr$pivot
    ## }                                   
    rankQR <- z0$rank
    singular.fit <- rankQR < p
    if (rankQR > 0) {
      if (singular.fit) {
        if (!singular.ok) 
          stop("singular fit encountered")
        pivot <- piv
        p1 <- pivot[seq_len(rankQR)]
        p2 <- pivot[(rankQR + 1):p]
        dn <- dimnames(x)
        x <- x[, p1]
        attr(x, "assign") <- assign[p1]
      }
            
      z <- loggammacenslmrob.fit(x=x, y=y, delta=delta, init=init, w=w, method=method, control=control)
      z$rank <- rankQR
      z$qr <- z0$qr
      z$rweights <- z$weights
      df <- length(y) - z$rank
      z$degree.freedom <- z$df.residual <- df
      z$delta <- save.d
      if (singular.fit) {
        coef <- numeric(p)
        coef[p2] <- NA
        coef[p1] <- z$coefficients
        d.p <- p - rankQR       
        names(coef) <- dn[[2L]]
        z$coefficients <- coef
        n <- NROW(y)
        z$qr[c("qr", "qraux", "pivot")] <- list(matrix(c(z$qr$qr, 
          rep.int(0, d.p * n)), n, p, dimnames = list(dn[[1L]], 
          dn[[2L]][piv])), c(z$qr$qraux, rep.int(0, d.p)), piv)
      }
      z$df <- c(rankQR+1, NROW(y)-rankQR-1, p+1) 
    } else {
      z <- list(coefficients = if (is.matrix(y)) matrix(NA, 
        p, ncol(y)) else rep.int(as.numeric(NA), p), 
        residuals = y, scale = NA,
        mu = NA, sigma = NA, lambda = NA,       
        fitted.values = 0 * y, cov = matrix(, 0, 0),
        rweights = rep.int(as.numeric(NA), 
        NROW(y)), weights = w, rank = 0, df.residual = NROW(y), 
        converged = TRUE, iter = 0, control = control)
        if (is.matrix(y)) 
          colnames(z$coefficients) <- colnames(x)
        else names(z$coefficients) <- colnames(x)
        if (!is.null(offset)) 
          z$residuals <- y - offset
    }
    if (!is.null(w)) {
      z$residuals <- z$residuals/wts
      z$fitted.values <- save.y - z$residuals
      z$weights <- w
      if (zero.weights) {
        coef <- z$coefficients
        coef[is.na(coef)] <- 0
        f0 <- cbind(1,x0) %*% coef
        save.r[ok] <- z$residuals
        save.r[nok] <- y0 - f0
        save.f[ok] <- z$fitted.values
        save.f[nok] <- f0
        z$residuals <- save.r
        z$fitted.values <- save.f
        rw <- z$rweights
        z$rweights <- rep.int(0, length(save.w))
        z$rweights[ok] <- rw
        z$weights <- save.w
        z$delta <- save.d
      }
    } else {
      z$weights <- rep(1, length(y))
    }
  }
  if (!is.null(offset)) 
    z$fitted.values <- z$fitted.values + offset
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- contrasts
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$assign <- assign
  if (control$compute.rd && !is.null(x)) 
    z$MD <- robMD(x, attr(mt, "intercept"), wqr = z$qr)
  if (model) 
    z$model <- mf
  if (ret.x) 
    z$x <- if (singular.fit || (!is.null(w) && zero.weights)) 
    model.matrix(mt, mf, contrasts)
  else x
  if (ret.y) 
    z$y <- if (!is.null(w)) 
      model.response(mf, "numeric")
    else y
  z$scale <- sqrt(z$sigma)
  class(z) <- "loggammacenslmrob"
  return(z)
}

#############################################################
#	loggammacenslmrob.fit function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 26, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

loggammacenslmrob.fit <- function(x, y, delta, init=NULL, w=NULL, control, ...) {
  if (is.null(w))
    w <- rep(1, length(y))
  x <- as.matrix(x)
  y <- na.omit(y)  
  if (!is.null(na <- attr(y, "na.action"))) {
    w <- w[-na]
    x <- x[-na,]
  }  
  if (control$method=="TQTau")
    result <- loggammacenslmrob.TQTau(x, y, delta, init, w, control)
  else if (control$method=="TWQTau")
    result <- loggammacenslmrob.TWQTau(x, y, delta, init, w, control)
  else if (control$method=="oneWL")
    result <- loggammacenslmrob.oneWL(x, y, delta, init, w, control)
  else if (control$method=="oneTML")
    result <- loggammacenslmrob.oneTML(x, y, delta, init, w, control)
  else if (control$method=="ML")
    result <- loggammacenslmrob.ML(x, y, delta, init, w, control)
  else if (control$method=="MM")
    result <- loggammacenslmrob.MM(x, y, delta, init, control)
  else if (control$method=="S")
    result <- loggammacenslmrob.S(x, y, delta, control)
  result$method <- control$method
  class(result) <- c("loggammacenslmrob.fit")
  return(result)
}

#############################################################
#	summary.loggammacenslmrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 26, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

########### SUMMARY
summary.loggammacenslmrob <- function(object, correlation = FALSE, ...) {
  if (is.null(object$terms)) 
    stop("invalid 'loggammacenslmrob' object:  no terms component")
  p <- object$rank
  df <- object$df.residual
  sigma <- object[["sigma"]]
  aliased <- is.na(coef(object))
  cf.nms <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (p > 0) {
    resid <- object$residuals
    pred <- object$fitted.values
    resp <- if (is.null(object[["y"]])) 
        pred + resid
      else object$y
    if ((object$method %in% c("oneTML", "oneWL", "ML"))) {
      keep <- object$rweights > object$control$minw
      if (object$method=="oneWL") {
        object$cov <- VCOV.WLoneCensloggammareg.empJ(y=resp[keep], delta=object$delta[keep], X=object$x[keep,,drop=FALSE], w=object$rweights[keep], mu=object$mu, sigma=object$sigma, lambda=object$lambda, beta=object$coefficients[-1])
        ####object$cov <- VCOV.WMLCensloggammareg(x=object$x[keep,,drop=FALSE], y=resp[keep], delta=object$delta[keep], mu=object$mu, sigma=object$sigma, lambda=object$lambda, beta=object$coefficients[-1], w=object$rweights[keep], xmax=object$control$xmax, minw=object$control$minw)
      } else if (object$method=="ML"){
        object$cov <- VCOV.WMLCensloggammareg(X=object$x[keep,,drop=FALSE], y=resp[keep], delta=object$delta[keep], mu=object$mu, sigma=object$sigma, lambda=object$lambda, beta=object$coefficients[-1], w=object$weights[keep], xmax=object$control$xmax, minw=object$control$minw)
      } else {
        object$cov <- VCOV.TMLoneCensloggammareg(y=resp,delta=object$delta,IX=cbind(1,object$x),w=object$rweights,coefficients=object$coefficients,sigma=object$sigma,lambda=object$lambda,cl=object$cut.lower,cu=object$cut.upper)
      }
      negvar <- which(diag(object$cov) <= 0)
      object$cov[negvar,] <- object$cov[,negvar] <- NA
    } else {
      object$cov <- diag(rep(NA,p+3))
    }
    n <- p + df
    p1 <- seq_len(p)
    se <- sqrt(if (length(object$cov) == 1L) object$cov else diag(object$cov))
    est <- c(object$mu,(object$coefficients[-1])[object$qr$pivot[p1]],object$sigma,object$lambda)
    names(est) <- c("(Intercept)", colnames(object$x)[object$qr$pivot[p1]], "Sigma", "Lambda")
    pe <- length(est)
    zval <- est/se
    zval[pe-1] <- NA
    ans <- object[c("call", "terms", "residuals", "scale",
      "sigma", "lambda", "df.residuals", "df", "weights",
      "rweights", "converged", "iter", "control")]
    if (!is.null(ans$weights)) 
      ans$residuals <- ans$residuals * sqrt(object$weights)
    ans$coefficients <- if (ans$converged) 
      cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
      else cbind(est, if (sigma <= 0) 
          0
        else NA, NA, NA)
    dimnames(ans$coefficients) <- list(names(est), cf.nms)
    ans$cov <- object$cov
    if (length(object$cov) > 1L) 
      dimnames(ans$cov) <- list(names(est), names(est))
      if (correlation) {
        ans$correlation <- ans$cov/outer(se, se)
      }
    } else {
      ans <- object
      ans$df <- c(0L, df, length(aliased))
      ans$coefficients <- matrix(NA, 0L, 4L, dimnames = list(NULL, cf.nms))
      ans$cov <- object$cov
    }
    ans$aliased <- aliased
    structure(ans, class = "summary.loggammacenslmrob")
}

#############################################################
#	print.summary.loggammacenslmrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 30, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

#PRINT.SUMMARY
print.summary.loggammacenslmrob <- function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call, width.cutoff = 72), 
        sep = "\n", collapse = "\n"), "\n", sep = "")
####    control <- lmrob.control.neededOnly(x$control)
  cat(" \\--> method = \"", x$control$method, "\"\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat(if (!is.null(x$weights) && diff(range(x$weights))) 
        "Weighted ", "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (NCOL(resid) > 1) 
        structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
                dimnames(resid)[[2]]))
      else setNames(quantile(resid), nam)
        print(rq, digits = digits, ...)
  } else print(resid, digits = digits, ...)
  if (length(x$aliased)) {
    if (!(x$converged)) {
      if (x$scale == 0) {
        cat("\nExact fit detected\n\nCoefficients:\n")
      } else {
        cat("\nAlgorithm did not converge\n\nCoefficients:\n")
      }
      printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
    } else {
      if (nsingular <- df[3L] - df[1L]) 
        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
      else cat("\nCoefficients:\n")
      coefs <- x$coefficients
      if (!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                  colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }
      printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "  ", ...)
      if (x$control$method!="ML")
        cat("\nRobust residual standard error:", format(signif(x$scale, digits)), "\n")
      else
        cat("\nResidual standard error:", format(signif(x$scale, digits)), "\n")
      
      correl <- x$correlation
      if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
          cat("\nCorrelation of Coefficients:\n")
          correl <- format(round(correl, 2), nsmall = 2, digits = digits)
          correl[!lower.tri(correl)] <- ""
          print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
      }
####      cat("Convergence in", x$iter, "IRWLS iterations\n")
    }
    if (x$control$method!="ML") {
      cat("\n")
      if (!is.null(rw <- x$rweights)) {
        if (any(zero.w <- x$weights == 0)) 
          rw <- rw[!zero.w]
        eps.outlier <- if (is.function(EO <- x$control$eps.outlier)) 
          EO(nobs.loggammacenslmrob(x))
        else EO
        summarizeRobWeights(rw, digits = digits, eps = eps.outlier, ...)
      }
    }
  } else cat("\nNo Coefficients\n")
####  if (showAlgo && !is.null(x$control))
####    robustbase:::printControl(x$control, digits = digits, drop. = "method")
  invisible(x)
}

#############################################################
#	print.loggammacenslmrob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and I. Locatelli
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: November, 30, 2015
#	Version: 0.1
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and I. Locatelli
#############################################################

print.loggammacenslmrob <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n", cl <- deparse(x$call, width.cutoff = 72), "\n", sep = "")
  if (!any(grepl("method *= *['\"]", cl))) 
    cat(" \\--> method = \"", x$control$method, "\"\n", sep = "")
  else cat("\n")
  if (length((cf <- coef(x)))) {
    if (x$converged) 
      cat("Coefficients:\n")
    else {
      if (x$scale == 0) {
        cat("Exact fit detected\n\nCoefficients:\n")
      } else {
        cat("Algorithm did not converge\n\n")
      }
    }
    print(format(cf, digits = digits), print.gap = 2, quote = FALSE)
    cat("Sigma:\n")
    print(format(x$sigma, digits = digits), print.gap = 2, quote = FALSE)
    cat("Lambda:\n")
    print(format(x$lambda, digits = digits), print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

nobs.loggammacenslmrob <- function (object, ...) 
if (!is.null(w <- object$weights)) sum(w != 0) else NROW(object$residuals)

## robMD <- function (x, intercept, wqr, ...) {
##   if (intercept == 1) 
##     x <- x[, -1, drop = FALSE]
##   if (ncol(x) >= 1) {
##     rob <- tryCatch(covMcd(x, ...), warning = function(w) structure("covMcd failed with a warning", 
##     class = "try-error", condition = w), error = function(e) structure("covMcd failed with an error", 
##     class = "try-error", condition = e))
##     if (inherits(rob, "try-error")) {
##       warning("Failed to compute robust Mahalanobis distances, reverting to robust leverages.")
##       return(lmrob.leverages(wqr = wqr))
##     }
##     sqrt(mahalanobis(x, rob$center, rob$cov))
##   }
## }

lmrob.leverages <- function (x, w = rep(1, NROW(x)), wqr = qr(sqrt(w) * x)) {
  if (missing(wqr) && !is.matrix(x)) 
    x <- as.matrix(x)
  pmin(1, rowSums(qr.Q(wqr)^2))
}
