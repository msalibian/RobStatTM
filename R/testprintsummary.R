print.f <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nClass: ", class(x), "\n", sep='')
  x <- x$DCML
  cat("\nCall:\n", cl <- deparse(x$call, width.cutoff=72), "\n", sep = "")
  control <- x$control
  if(length((cf <- coef(x)))) {
    if( x$converged )
      cat("Coefficients:\n")
    else {
      if (x$scale == 0) {
        cat("Exact fit detected\n\nCoefficients:\n")
      } else {
        cat("Algorithm did not converge\n\n")
        if (control$method == "S")
          cat("Coefficients of the *initial* S-estimator:\n")
        else
          cat(sprintf("Coefficients of the %s-estimator:\n",
                      control$method))
      }
    }
    print(format(cf, digits = digits), print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

summary.f <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
  cat("\nIn summary.f\n")
  object <- object$DCML
  if (is.null(object$terms))
    stop("invalid 'lmrobdet' object:  no terms component")
  p <- object$rank
  df <- object$df.residual #was $degree.freedom
  sigma <- object[["scale"]]
  aliased <- is.na(coef(object))
  cf.nms <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  if (p > 0) {
    n <- p + df
    p1 <- seq_len(p)
    se <- sqrt(if(length(object$cov) == 1L) object$cov else diag(object$cov))
    est <- object$coefficients[object$qr$pivot[p1]]
    tval <- est/se
    ans <- object[c("call", "terms", "residuals", "scale",
                    "converged", "iter", "control")]
    ## 'df' vector, modeled after summary.lm() : ans$df <- c(p, rdf, NCOL(Qr$qr))
    ## where  p <- z$rank ; rdf <- z$df.residual ; Qr <- qr.lm(object)
    ans$df <- c(p, df, NCOL(object$qr$qr))
    ans$coefficients <-
      if( ans$converged)
        cbind(est, se, tval, 2 * pt(abs(tval), df, lower.tail = FALSE))
    else
      cbind(est, if(sigma <= 0) 0 else NA, NA, NA)
    dimnames(ans$coefficients) <- list(names(est), cf.nms)
    ans$r.squared <- ans$adj.r.squared <- NULL
    ans$cov <- object$cov
    if(length(object$cov) > 1L)
      dimnames(ans$cov) <- dimnames(ans$coefficients)[c(1,1)]
    if (correlation) {
      ans$correlation <- ans$cov / outer(se, se)
      ans$symbolic.cor <- symbolic.cor
    }
  } else { ## p = 0: "null model"
    ans <- object
    ans$df <- c(0L, df, length(aliased))
    ans$coefficients <- matrix(NA, 0L, 4L, dimnames = list(NULL, cf.nms))
    ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov <- object$cov
  }
  ans$aliased <- aliased # used in print method
  ans$sigma <- sigma # 'sigma': in summary.lm() & 'fit.models' pkg
  structure(ans, class = "summary.f")
}



print.summary.f <- function (x, digits = max(3, getOption("digits") - 3),
                                    symbolic.cor = x$symbolic.cor,
                                    signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nIn print.summary.f\n")
  cat("\nCall:\n",
      paste(deparse(x$call, width.cutoff=72), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  control <- x$control
  # cat(" \\--> method = \"", control$method, '"\n', sep = "")
  ## else cat("\n")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat(if (!is.null(x$weights) && diff(range(x$weights))) "Weighted ",
      "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <-
      if (NCOL(resid) > 1)
        structure(apply(t(resid), 1, quantile),
                  dimnames = list(nam, dimnames(resid)[[2]]))
    else setNames(quantile(resid), nam)
    print(rq, digits = digits, ...)
  }
  else print(resid, digits = digits, ...)
  ## FIXME: need to catch rdf == 0?
  if( length(x$aliased) ) {
    if( !(x$converged) ) {
      if (x$scale == 0) {
        cat("\nExact fit detected\n\nCoefficients:\n")
      } else {
        cat("\nAlgorithm did not converge\n")
        if (control$method == "S")
          cat("\nCoefficients of the *initial* S-estimator:\n")
        else
          cat(sprintf("\nCoefficients of the %s-estimator:\n",
                      control$method))
      }
      printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
                   ...)
    } else {
      if (nsingular <- df[3L] - df[1L])
        cat("\nCoefficients: (", nsingular,
            " not defined because of singularities)\n", sep = "")
      else cat("\nCoefficients:\n")
      coefs <- x$coefficients
      if(!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }
      
      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                   na.print="NA", ...)
      cat("\nRobust residual standard error:",
          format(signif(x$scale, digits)),"\n")
      if (!is.null(x$r.squared) && x$df[1] != attr(x$terms, "intercept")) {
        cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
        cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits),
            "\n")
      }
      ## FIXME: use naprint() here to list observations deleted due to missingness?
      correl <- x$correlation
      if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
          cat("\nCorrelation of Coefficients:\n")
          if (is.logical(symbolic.cor) && symbolic.cor) {
            print(symnum(correl), abbr.colnames = NULL)
          }
          else { correl <- format(round(correl, 2), nsmall = 2,
                                  digits = digits)
          correl[!lower.tri(correl)] <- ""
          print(correl[-1, -p, drop = FALSE], quote = FALSE)
          }
        }
      }
      cat("Convergence in", x$iter, "IRWLS iterations\n")
    }
    cat("\n")
    
    # if (!is.null(rw <- x$rweights)) {
    #   if (any(zero.w <- x$weights == 0))
    #     rw <- rw[!zero.w]
    #   eps.outlier <- if (is.function(EO <- control$eps.outlier))
    #     EO(nobs(x)) else EO
    #   summarizeRobWeights(rw, digits = digits, eps = eps.outlier, ...)
    # }
    
  } else cat("\nNo Coefficients\n")
  invisible(x)
}

