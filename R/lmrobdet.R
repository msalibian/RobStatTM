#' Robust linear regression estimators
#'
#' This function computes an MM-regression estimators for linear models
#' using deterministic starting points.
#'
#' This function computes MM-regression estimators
#' computed using Pen~a-Yohai candidates (instead of subsampling ones).
#' This function makes use of the functions \code{lmrob.fit},
#' \code{lmrob..M..fit}, \code{.vcov.avar1}, \code{lmrob.S} and
#' \code{lmrob.lar}, from robustbase,
#' along with utility functions used by these functions,
#' modified so as to include use of the analytic form of the
#' optimal psi and rho functions (for the optimal psi function , see
#' Section 5.8.1 of Maronna, Martin, Yohai and Salibian Barrera, 2019)
#'
#' @param formula a symbolic description of the model to be fit.
#' @param data an optional data frame, list or environment containing
#' the variables in the model. If not found in \code{data}, model variables
#' are taken from \code{environment(formula)}, which usually is the
#' root environment of the current R session.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param weights an optional vector of weights to be used in the fitting process.
#' @param na.action a function to indicates what should happen when the data contain NAs.
#' The default is set by the \link{na.action} setting of \code{\link[base]{options}}, and is
#' \code{na.fail} if that is unset.
#' @param model logical value indicating whether to return the model frame
#' @param x logical value indicating whether to return the model matrix
#' @param y logical value indicating whether to return the vector of responses
#' @param singular.ok logical value. If \code{FALSE} a singular fit produces an error.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \link{model.matrix.default}.
#' @param offset this can be used to specify an a priori known component to be included
#' in the linear predictor during fitting. An offset term can be included in the formula
#' instead or as well, and if both are specified their sum is used.
#' @param control a list specifying control parameters as returned by the function
#' \link{lmrobdet.control}.
#'
#' @return A list with the following components:
#' \item{coefficients}{The estimated vector of regression coefficients}
#' \item{scale}{The robust residual M-scale estimate using the final residuals from the converged iterated weighted least square (IRWLS) algorithm final estimate}
#' \item{residuals}{The vector of residuals associated with the robust fit}
#' \item{loss}{Value of the objective function at the final MM-estimator}
#' \item{converged}{Logical value indicating whether IRWLS iterations for the MM-estimator have converged}
#' \item{iter}{Number of IRWLS iterations for the MM-estimator}
#' \item{rweights}{Robustness weights for the MM-estimator}
#' \item{fitted.values}{Fitted values associated with the robust fit}
#' \item{rank}{Numeric rank of the fitted linear model}
#' \item{cov}{The estimated covariance matrix of the regression estimates}
#' \item{df.residual}{The residual degrees of freedom}
#' \item{degree.freedom}{The residual degrees of freedom}
#' \item{scale.S}{Minimum robust scale associated with the preliminary highly robust but inefficient S-estimator.}
#' \item{r.squared}{The robust multiple correlation coefficient}
#' \item{adj.r.squared}{The adjusted robust multiple correlation coefficient taking into account the degrees of freedom of each term}
#' \item{contrasts}{(only where relevant) the contrasts used}
#' \item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting}
#' \item{call}{the matched call}
#' \item{model}{if requested, the model frame used}
#' \item{x}{if requested, the model matrix used}
#' \item{y}{if requested, the response vector used}
#' \item{terms}{The \link{terms} object used.}
#' \item{assign}{Used to separate continuous from categorical columns in the design matrix}
#' \item{na.action}{(where relevant) information returned by model.frame on the special handling of NAs}
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, based on \code{lmrob} from package \code{robustbase}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \link{DCML}, \link{MMPY}, \link{SMPY}
#'
#' @examples
#' data(coleman, package='robustbase')
#' m2 <- lmrobdetMM(Y ~ ., data=coleman)
#' m2
#' summary(m2)
#'
#' @import stats
#' @useDynLib RobStatTM, .registration = TRUE
#' @export
lmrobdetMM <- function(formula, data, subset, weights, na.action,
                       model = TRUE, x = !control$compute.rd, y = FALSE,
                       singular.ok = TRUE, contrasts = NULL, offset = NULL,
                       control = lmrobdet.control())
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset) && length(offset) != NROW(y))
    stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                  length(offset), NROW(y)), domain = NA)

  if (is.empty.model(mt)) {
    x <- NULL
    singular.fit <- FALSE ## to avoid problems below
    z <- list(coefficients = if (is.matrix(y)) matrix(,0,3) else numeric(0),
              residuals = y, scale = NA, fitted.values = 0 * y,
              cov = matrix(,0,0), weights = w, rank = 0,
              df.residual = NROW(y), converged = TRUE, iter = 0)
    if(!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
      z$offset <- offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    contrasts <- attr(x, "contrasts")
    assign <- attr(x, "assign")
    p <- ncol(x)
    if(!is.null(offset))
      y <- y - offset
    if (!is.null(w)) {
      ## checks and code copied/modified from lm.wfit
      ny <- NCOL(y)
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
        ok <- w != 0
        nok <- !ok
        w <- w[ok]
        x0 <- x[nok, , drop = FALSE]
        x  <- x[ ok, , drop = FALSE]
        n <- nrow(x)
        y0 <- if (ny > 1L) y[nok, , drop = FALSE] else y[nok]
        y  <- if (ny > 1L) y[ ok, , drop = FALSE] else y[ok]
        ## add this information to model.frame as well
        ## need it in outlierStats()
        ## ?? could also add this to na.action, then
        ##    naresid() would pad these as well.
        attr(mf, "zero.weights") <- which(nok)
      }
      wts <- sqrt(w)
      save.y <- y
      x <- wts * x
      y <- wts * y
    }
    ## check for singular fit

    if(getRversion() >= "3.1.0") {
      z0 <- .lm.fit(x, y) #, tol = control$solve.tol)
      piv <- z0$pivot
    } else {
      z0 <- lm.fit(x, y) #, tol = control$solve.tol)
      piv <- z0$qr$pivot
    }
    rankQR <- z0$rank

    singular.fit <- rankQR < p
    if (rankQR > 0) {
      if (singular.fit) {
        if (!singular.ok) stop("singular fit encountered")
        pivot <- piv
        p1 <- pivot[seq_len(rankQR)]
        p2 <- pivot[(rankQR+1):p]
        ## to avoid problems in the internal fitting methods,
        ## split into singular and non-singular matrices,
        ## can still re-add singular part later
        dn <- dimnames(x)
        x <- x[,p1]
        attr(x, "assign") <- assign[p1] ## needed for splitFrame to work
      }
      # Check if there are factors
      if( control$initial=="SM" ) {
        split <- splitFrame(mf, x, control$split.type)
        if (ncol(split$x1) == 0) {
          control$initial <- 'S'
          warning("No categorical variables found in model. Reverting to an MM-estimator.")
        }
      }
      if( control$initial=="SM" ) {
        z <- SMPY(mf=mf, y=y, control=control, split=split)
      } else if( control$initial == "S" ) {
        z <- MMPY(X=x, y=y, control=control, mf=mf)
      } else stop('Unknown value for lmrobdet.control()$initial')
      # update residual scale estimator
      # re.dcml <- as.vector(y - x %*% beta.dcml)
      # si.dcml.final <- mscale(u=re.dcml, tol = control$mscale.tol, delta=dee, tuning.chi=control$tuning.chi)
      n <- length(z$resid)
      z$scale.S <- z$scale
      z$scale <- mscale(u=z$resid, tol = control$mscale_tol, delta=control$bb*(1-p/length(z$resid)), tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
      # compute robust R^2
      r.squared <- adj.r.squared <- NA
      if(z$scale > .Machine$double.eps) {
      # if(control$family == 'bisquare') {
      s2 <- mean(rho(z$resid/z$scale, family = control$family, cc=control$tuning.psi))
      if( p != attr(mt, "intercept") ) {
        df.int <- if (attr(mt, "intercept"))
          1L
        else 0L
        if(df.int == 1L) {
          tmp <- as.vector(refine.sm(x=matrix(rep(1,n), n, 1), y=y, initial.beta=median(y),
                                     initial.scale=z$scale, k=500,
                                     conv=1, family = control$family, cc = control$tuning.psi, step='M')$beta.rw)
          s02 <- mean(rho((y-tmp)/z$scale, family = control$family, cc=control$tuning.psi))
        } else {
          s02 <- mean(rho(y/z$scale, family = control$family, cc=control$tuning.psi))
        }
        r.squared <- INVTR2( (s02 - s2)/(s02*(1-s2)), control$family, control$tuning.psi)
        adj.r.squared <- ((n-1)/(n-p))*r.squared -(p-1)/(n-p) # ( s02/(n-1) - s2/(n-z$rank) ) / (s02/(n-1)) # n-p? p.193
      }
      }
      # }
      z$r.squared <- r.squared
      z$adj.r.squared <- adj.r.squared
      # DCML
      # LS is already computed in z0
      # z2 <- DCML(x=x, y=y, z=z, z0=z0, control=control)
      # z$MM <- z
      # # complete the MM object
      # z$MM$na.action <- attr(mf, "na.action")
      # z$MM$offset <- offset
      # z$MM$contrasts <- contrasts
      # z$MM$xlevels <- .getXlevels(mt, mf)
      # z$MM$call <- cl
      # z$MM$terms <- mt
      # z$MM$assign <- assign
      if(control$compute.rd && !is.null(x))
        z$MD <- robustbase::robMD(x, attr(mt, "intercept"), wqr=z$qr)
      if(model)
        z$model <- mf
      if(ret.x)
        z$x <- if (singular.fit || (!is.null(w) && zero.weights))
          model.matrix(mt, mf, contrasts) else x
      if (ret.y)
        z$y <- if (!is.null(w)) model.response(mf, "numeric") else y
      #
      #       z$coefficients <- z2$coefficients
      #       z$scale <- z2$scale
      #       z$residuals <- z2$residuals
      #       z$cov <- z2$cov
      #       z$fitted.values <- y - z2$residuals
      #       z$rweights <- z$loss <- NULL

      # z <- lmrob.fit(x, y, control, init=init, mf = mf) #-> ./lmrob.MM.R
      if (singular.fit) {
        coef <- numeric(p)
        coef[p2] <- NA
        coef[p1] <- z$coefficients
        names(coef) <- dn[[2L]]
        z$coefficients <- coef
        ## Update QR decomposition (z$qr)
        ## pad qr and qraux with zeroes (columns that were pivoted to the right in z0)
        d.p <- p-rankQR
        n <- NROW(y)
        z$qr[c("qr","qraux","pivot")] <-
          list(matrix(c(z$qr$qr, rep.int(0, d.p*n)), n, p,
                      dimnames = list(dn[[1L]], dn[[2L]][piv])),
               ## qraux:
               c(z$qr$qraux, rep.int(0, d.p)),
               ## pivot:
               piv)
      }
    } else { ## rank 0
      z <- list(coefficients = if (is.matrix(y)) matrix(NA,p,ncol(y))
                else rep.int(as.numeric(NA), p),
                residuals = y, scale = NA, fitted.values = 0 * y,
                cov = matrix(,0,0), rweights = rep.int(as.numeric(NA), NROW(y)),
                weights = w, rank = 0, df.residual = NROW(y),
                converged = TRUE, iter = 0, control=control)
      if (is.matrix(y)) colnames(z$coefficients) <- colnames(x)
      else names(z$coefficients) <- colnames(x)
      if(!is.null(offset)) z$residuals <- y - offset
    }
    if (!is.null(w)) {
      z$residuals <- z$residuals/wts
      z$fitted.values <- save.y - z$residuals
      z$weights <- w
      if (zero.weights) {
        coef <- z$coefficients
        coef[is.na(coef)] <- 0
        f0 <- x0 %*% coef
        if (ny > 1) {
          save.r[ok, ] <- z$residuals
          save.r[nok, ] <- y0 - f0
          save.f[ok, ] <- z$fitted.values
          save.f[nok, ] <- f0
        }
        else {
          save.r[ok] <- z$residuals
          save.r[nok] <- y0 - f0
          save.f[ok] <- z$fitted.values
          save.f[nok] <- f0
        }
        z$residuals <- save.r
        z$fitted.values <- save.f
        z$weights <- save.w
        rw <- z$rweights
        z$rweights <- rep.int(0, length(save.w))
        z$rweights[ok] <- rw
      }
    }
  }
  if(!is.null(offset))
    z$fitted.values <- z$fitted.values + offset

  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- contrasts
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$assign <- assign
  if(control$compute.rd && !is.null(x))
    z$MD <- robustbase::robMD(x, attr(mt, "intercept"), wqr=z$qr)
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- if (singular.fit || (!is.null(w) && zero.weights))
      model.matrix(mt, mf, contrasts) else x
  if (ret.y)
    z$y <- if (!is.null(w)) model.response(mf, "numeric") else y
  # class(z) <- c("lmrobdet", "lmrob")
  # z
  # tmp <- z
  # tmp$MM <- NULL
  # z2 <- list(DCML=tmp, MM=z$MM)
  # class(z2$DCML) <- c("DCML", "lmrob")
  # class(z2$MM) <- "lmrob"
  # class(z2) <- "lmrobdet"
  # z2
  class(z) <- c('lmrobdetMM', 'lmrob')
  z
}


#' Tuning parameters for lmrobdetMM and lmrobdetDCML
#'
#' This function sets tuning parameters for the MM estimator implemented in \code{lmrobdetMM} and
#' the Distance Constrained Maximum Likelihood regression estimators
#' computed by \code{lmrobdetDCML}.
#'
#' @rdname lmrobdet.control
#' @param tuning.chi tuning constant for the function used to compute the M-scale
#' used for the initial S-estimator. If missing, it is computed inside \code{lmrobdet.control} to match
#' the value of \code{bb} according to the family of rho functions specified in \code{family}.
#' @param bb tuning constant (between 0 and 1/2) for the M-scale used to compute the initial S-estimator. It
#' determines the robusness (breakdown point) of the resulting MM-estimator, which is
#' \code{bb}. Defaults to 0.5.
#' @param tuning.psi tuning parameters for the regression M-estimator computed with a rho function
#' as specified with argument \code{family}. If missing, it is computed inside \code{lmrobdet.control} to match
#' the value of \code{efficiency} according to the family of rho functions specified in \code{family}.
#' Appropriate values for \code{tuning.psi} for a given desired efficiency for Gaussian errors
#' can be constructed using the functions \link{bisquare}, \link{mopt} and \link{opt}.
#' @param efficiency desired asymptotic efficiency of the final regression M-estimator. Defaults to 0.95.
#' @param max.it maximum number of IRWLS iterations for the MM-estimator
#' @param refine.tol relative covergence tolerance for the S-estimator
#' @param rel.tol relative covergence tolerance for the IRWLS iterations for the MM-estimator
#' @param refine.PY number of refinement steps for the Pen~a-Yohai candidates
#' @param solve.tol (for the S algorithm): relative tolerance for matrix inversion. Hence, this corresponds to \code{\link{solve.default}}'s tol.
#' @param trace.lev positive values (increasingly) provide details on the progress of the MM-algorithm
#' @param compute.rd logical value indicating whether robust leverage distances need to be computed.
#' @param family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "opt" and "mopt"). Incomplete entries will be matched to the current valid options. Defaults to "mopt".
#' @param corr.b logical value indicating whether a finite-sample correction should be applied
#' to the M-scale parameter \code{bb}.
#' @param split.type determines how categorical and continuous variables are split. See
#' \code{\link[robustbase]{splitFrame}}.
#' @param initial string specifying the initial value for the M-step of the MM-estimator. Valid
#' options are \code{'S'}, for an S-estimator and \code{'MS'} for an M-S estimator which is
#' appropriate when there are categorical explanatory variables in the model.
#' @param psc_keep For \code{pyinit}, proportion of observations to remove based on PSCs. The effective proportion of removed
#' observations is adjusted according to the sample size to be \code{prosac*(1-p/n)}. See \code{\link{pyinit}}.
#' @param resid_keep_method For \code{pyinit}, how to clean the data based on large residuals. If
#' \code{"threshold"}, all observations with scaled residuals larger than \code{C.res} will
#' be removed, if \code{"proportion"}, observations with the largest \code{prop} residuals will
#' be removed. See \code{\link{pyinit}}.
#' @param resid_keep_thresh See parameter \code{resid_keep_method} above. See \code{\link{pyinit}}.
#' @param resid_keep_prop See parameter \code{resid_keep_method} above. See \code{\link{pyinit}}.
#' @param py_maxit Maximum number of iterations. See \code{\link{pyinit}}.
#' @param py_eps Relative tolerance for convergence.  See \code{\link{pyinit}}.
#' @param mscale_maxit Maximum number of iterations for the M-scale algorithm. See \code{\link{pyinit}} and \code{\link{mscale}}.
#' @param mscale_tol Convergence tolerance for the M-scale algorithm. See \code{\link{mscale}} and \code{\link{mscale}}.
#' @param mscale_rho_fun String indicating the loss function used for the M-scale. See \code{\link{pyinit}}.
#'
#' @return A list with the necessary tuning parameters.
#'
#' @details The argument \code{family} specifies the name of the family of loss function to be used. Current valid
#' options are "bisquare", "opt" and "mopt"--"opt" refers to the optimal psi function defined in Section 5.8.1. of the
#' book Robust Statistics: Theory and Methods (with R) by Maronna, Martin, Yohai and Salibian-Barrera,
#' "mopt" is a modified  version of the optimal psi function to make it
#' strictly increasing close to 0, and to make the corresponding weight function
#' non-increasing near 0.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @seealso \code{\link{pyinit}}, \code{\link{mscale}}.
#'
#' @examples
#' data(coleman, package='robustbase')
#' m2 <- lmrobdetMM(Y ~ ., data=coleman, control=lmrobdet.control(refine.PY=50))
#' m2
#' summary(m2)
#'
#' @export
lmrobdet.control <- function(bb = 0.5,
                             efficiency = 0.95,
                             family = 'mopt',
                             tuning.psi,
                             tuning.chi,
                             compute.rd = FALSE,
                             corr.b = TRUE,
                             split.type = "f",
                             initial='S',
                             max.it = 100, refine.tol = 1e-7, rel.tol = 1e-7,
                             refine.PY = 10,
                             solve.tol = 1e-7, trace.lev = 0,
                             psc_keep = 0.5, resid_keep_method = 'threshold',
                             resid_keep_thresh = 2, resid_keep_prop = .2, py_maxit = 20, py_eps = 1e-5,
                             mscale_maxit = 50, mscale_tol = 1e-06, mscale_rho_fun = 'bisquare')
{
  family <- match.arg(family, choices = FAMILY.NAMES)
  if( (efficiency > .9999 ) & ( (family=='mopt') | (family=='opt') ) ) {
    efficiency <- .9999
    warning("Current implementation of \'opt\' or \'mopt\' only allows efficiencies up to 99.99%. Efficiency set to 99.99% for this call.")
  }
  if(missing(tuning.psi))
    tuning.psi <- do.call(family, args=list(e=efficiency))
  if( (length(tuning.psi) == 1) & is.null(names(tuning.psi)) )
    tuning.psi <- c( 'c' = tuning.psi )
  if(missing(tuning.chi))
    tuning.chi <- adjustTuningVectorForBreakdownPoint(family=family, cc=tuning.psi, breakdown.point = bb)
  return(list(family=family, # psi=psi,
              tuning.chi=tuning.chi, bb=bb, tuning.psi=tuning.psi,
              max.it=max.it,
              refine.tol=refine.tol,
              corr.b = corr.b, refine.PY = refine.PY,
              rel.tol=rel.tol,
              solve.tol=solve.tol, trace.lev=trace.lev,
              compute.rd=compute.rd,
              split.type=split.type,
              initial=initial, # method=method, subsampling=subsampling,
              psc_keep=psc_keep, resid_keep_method=resid_keep_method, resid_keep_thresh=resid_keep_thresh,
              resid_keep_prop=resid_keep_prop, py_maxit=py_maxit, py_eps=py_eps, mscale_maxit=mscale_maxit,
              mscale_tol=mscale_tol, mscale_rho_fun='bisquare', mts=1000)) # mts maximum number of subsamples. Un-used, but passed (unnecessarily) to the function that performs M-iterations (lmrob..M..fit), so set here.
}


#' @export
print.lmrobdetMM <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  # x <- x$DCML
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
        # if (control$method == "S")
        # cat("Coefficients of the *initial* S-estimator:\n")
        # else
        cat(sprintf("Coefficients of the %s-estimator:\n",
                    control$method))
      }
    }
    print(format(cf, digits = digits), print.gap = 2, quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' @export
summary.lmrobdetMM <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
  # object <- object$DCML
  if (is.null(object$terms))
    stop("invalid 'lmrobdetMM' object:  no terms component")
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
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), df, lower.tail = FALSE))
    #   if( ans$converged)
    #    cbind(est, se, tval, 2 * pt(abs(tval), df, lower.tail = FALSE))
    # else
    #   cbind(est, if(sigma <= 0) 0 else NA, NA, NA)
    dimnames(ans$coefficients) <- list(names(est), cf.nms)
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
  ans$r.squared <- object$r.squared
  ans$adj.r.squared <- object$adj.r.squared
  structure(ans, class = "summary.lmrobdetMM")
}

#' @export
print.summary.lmrobdetMM <- function (x, digits = max(3, getOption("digits") - 3),
                                      symbolic.cor = x$symbolic.cor,
                                      signif.stars = getOption("show.signif.stars"), ...)
{
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
      if (!is.null(x$r.squared) && !is.na(x$r.squared) && x$df[1] != attr(x$terms, "intercept")) {
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


#' IRWLS iterations for S- or M-estimators
#'
#' This function performs iterative improvements for S- or
#' M-estimators.
#'
#' This function performs iterative improvements for S- or
#' M-estimators. Both iterations are formally the same, the
#' only difference is that for M-iterations the residual
#' scale estimate remains fixed, while for S-iterations
#' it is updated at each step. In this case, we follow
#' the Fast-S algorithm of Salibian-Barrera and Yohai
#' an use one step updates for the M-scale, as opposed
#' to a full computation. This as internal function.
#'
#' @param x design matrix
#' @param y vector of responses
#' @param initial.beta vector of initial regression estimates
#' @param initial.scale initial residual scale estimate. If missing the (scaled) median of
#' the absolute residuals is used.
#' @param k maximum number of refining steps to be performed
#' @param conv an integer indicating whether to check for convergence (1) at each step,
#' or to force running k steps (0)
#' @param b tuning constant for the M-scale estimator, used if iterations are for an S-estimator.
#' @param cc tuning constant for the rho function.
#' @param family string specifying the name of the family of loss function to be used (current
#' valid options are "bisquare", "opt" and "mopt")
#' @param step a string indicating whether the iterations are to compute an S-estiamator
#' ('S') or an M-estimator ('M')
#'
#' @return A list with the following components:
#' \item{beta.rw}{The updated vector of regression coefficients}
#' \item{scale.rw}{The corresponding estimated residual scale}
#' \item{converged}{A logical value indicating whether the algorithm
#' converged}
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}.
#' @export
refine.sm <- function(x, y, initial.beta, initial.scale, k=50,
                      conv=1, b, cc, family, step='M') {

  #refine.sm <- function(x, y, initial.beta, initial.scale, k=50,
  #                     conv=1, b, cc, step='M') {

  ## Weight function   # weight function = psi(u)/u
  #f.w <- function(u, cc) {
  #  tmp <- (1 - (u/cc)^2)^2
  #  tmp[abs(u/cc) > 1] <- 0
  #  return(tmp)
  #}
  f.w <- function(u, family, cc)
    Mwgt(x = u, cc = cc, psi = family)


  n <- dim(x)[1]
  # p <- dim(x)[2]

  res <- as.vector( y - x %*% initial.beta )

  if( missing( initial.scale ) ) {
    initial.scale <- scale <- median(abs(res))/.6745
  } else {
    scale <- initial.scale
  }

  beta <- initial.beta


  converged <- FALSE

  # lower.bound <- median(abs(res))/cc

  if( scale == 0) {
    beta.1 <- initial.beta
    scale <- initial.scale
  } else {
    for(i in 1:k) {
      # do one step of the iterations to solve for the scale
      scale.super.old <- scale
      #lower.bound <- median(abs(res))/1.56
      if(step=='S') {
        scale <- sqrt( scale^2 * mean( rho( res / scale, family = family, cc = cc ) ) / b     )
        # scale <- mscale(res, tol=1e-7, delta=b, max.it=500, tuning.chi=cc)
      }
      # now do one step of IRWLS with the "improved scale"
      weights <- f.w( res/scale, family = family, cc = cc )
      # W <- matrix(weights, n, p)
      xw <- x * sqrt(weights) # sqrt(W)
      yw <- y *   sqrt(weights)
      beta.1 <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
      if(any(is.na(beta.1))) {
        beta.1 <- initial.beta
        scale <- initial.scale
        break
      }
      if( (conv==1) ) {
        # check for convergence
        if( norm.sm( beta - beta.1 ) / norm.sm(beta) < 1e-7 ) { # magic number alert!!!
          converged <- TRUE
          break
        }
      }
      res <- as.vector( y - x %*% beta.1 )
      beta <- beta.1
      # print(as.vector(t(x) %*% rhoprime(res/scale, cc))/n)
      # print(scale)
    }
  }
  # res <- as.vector( y - x %*% beta )
  # get the residuals from the last beta
  return(list(beta.rw = beta.1, scale.rw = scale, converged=converged))
}

norm.sm <- function(x) sqrt(sum(x^2))



our.solve <- function(a,b) {
  a <- qr(a)
  da <- dim(a$qr)
  if(a$rank < (p <- da[2]))
    return(NA)
  else qr.coef(a, b)
}



#' Robust Distance Constrained Maximum Likelihood estimators for linear regression
#'
#' This function computes robust Distance Constrained Maximum Likelihood
#' estimators for linear models.
#'
#' This function computes Distance Constrained Maximum Likelihood regression estimators
#' computed using an MM-regression estimator based on Pen~a-Yohai
#' candidates (instead of subsampling ones).
#' This function makes use of the functions \code{lmrob.fit},
#' \code{lmrob..M..fit}, \code{.vcov.avar1}, \code{lmrob.S} and
#' \code{lmrob.lar}, from robustbase,
#' along with utility functions used by these functions,
#' modified so as to include use of the analytic form of the
#' optimal psi and rho functions (for the optimal psi function , see
#' Section 5.8.1 of Maronna, Martin, Yohai and Salibian Barrera, 2019)
#'
#' @param formula a symbolic description of the model to be fit.
#' @param data an optional data frame, list or environment containing
#' the variables in the model. If not found in \code{data}, model variables
#' are taken from \code{environment(formula)}, which usually is the
#' root environment of the current R session.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param weights an optional vector of weights to be used in the fitting process.
#' @param na.action a function to indicates what should happen when the data contain NAs.
#' The default is set by the \link{na.action} setting of \code{\link[base]{options}}, and is
#' \code{na.fail} if that is unset.
#' @param model logical value indicating whether to return the model frame
#' @param x logical value indicating whether to return the model matrix
#' @param y logical value indicating whether to return the vector of responses
#' @param singular.ok logical value. If \code{FALSE} a singular fit produces an error.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \link{model.matrix.default}.
#' @param offset this can be used to specify an a priori known component to be included
#' in the linear predictor during fitting. An offset term can be included in the formula
#' instead or as well, and if both are specified their sum is used.
#' @param control a list specifying control parameters as returned by the function
#' \link{lmrobdet.control}.
#'
#' @return A list with the following components:
#' \item{coefficients}{The estimated vector of regression coefficients}
#' \item{scale}{The estimated scale of the residuals}
#' \item{residuals}{The vector of residuals associated with the robust fit}
#' \item{converged}{Logical value indicating whether IRWLS iterations for the MM-estimator have converged}
#' \item{iter}{Number of IRWLS iterations for the MM-estimator}
#' \item{rweightsMM}{Robustness weights for the MM-estimator}
#' \item{fitted.values}{Fitted values associated with the robust fit}
#' \item{rank}{Numeric rank of the fitted linear model}
#' \item{cov}{The estimated covariance matrix of the regression estimates}
#' \item{df.residual}{The residual degrees of freedom}
#' \item{contrasts}{(only where relevant) the contrasts used}
#' \item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting}
#' \item{call}{the matched call}
#' \item{model}{if requested, the model frame used}
#' \item{x}{if requested, the model matrix used}
#' \item{y}{if requested, the response vector used}
#' \item{na.action}{(where relevant) information returned by model.frame on the special handling of NAs}
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, based on \code{lmrob}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @examples
#' data(coleman, package='robustbase')
#' m1 <- lmrobdetDCML(Y ~ ., data=coleman)
#' m1
#' summary(m1)
#'
#' @export
lmrobdetDCML <- function(formula, data, subset, weights, na.action,
                         model = TRUE, x = !control$compute.rd, y = FALSE,
                         singular.ok = TRUE, contrasts = NULL, offset = NULL,
                         control = lmrobdet.control())
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset) && length(offset) != NROW(y))
    stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                  length(offset), NROW(y)), domain = NA)

  if (is.empty.model(mt)) {
    x <- NULL
    singular.fit <- FALSE ## to avoid problems below
    z <- list(coefficients = if (is.matrix(y)) matrix(,0,3) else numeric(0),
              residuals = y, scale = NA, fitted.values = 0 * y,
              cov = matrix(,0,0), weights = w, rank = 0,
              df.residual = NROW(y), converged = TRUE, iter = 0)
    if(!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
      z$offset <- offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    contrasts <- attr(x, "contrasts")
    assign <- attr(x, "assign")
    p <- ncol(x)
    if(!is.null(offset))
      y <- y - offset
    if (!is.null(w)) {
      ## checks and code copied/modified from lm.wfit
      ny <- NCOL(y)
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
        ok <- w != 0
        nok <- !ok
        w <- w[ok]
        x0 <- x[nok, , drop = FALSE]
        x  <- x[ ok, , drop = FALSE]
        n <- nrow(x)
        y0 <- if (ny > 1L) y[nok, , drop = FALSE] else y[nok]
        y  <- if (ny > 1L) y[ ok, , drop = FALSE] else y[ok]
        ## add this information to model.frame as well
        ## need it in outlierStats()
        ## ?? could also add this to na.action, then
        ##    naresid() would pad these as well.
        attr(mf, "zero.weights") <- which(nok)
      }
      wts <- sqrt(w)
      save.y <- y
      x <- wts * x
      y <- wts * y
    }
    ## check for singular fit

    if(getRversion() >= "3.1.0") {
      z0 <- .lm.fit(x, y) #, tol = control$solve.tol)
      piv <- z0$pivot
    } else {
      z0 <- lm.fit(x, y) #, tol = control$solve.tol)
      piv <- z0$qr$pivot
    }
    rankQR <- z0$rank

    singular.fit <- rankQR < p
    if (rankQR > 0) {
      if (singular.fit) {
        if (!singular.ok) stop("singular fit encountered")
        pivot <- piv
        p1 <- pivot[seq_len(rankQR)]
        p2 <- pivot[(rankQR+1):p]
        ## to avoid problems in the internal fitting methods,
        ## split into singular and non-singular matrices,
        ## can still re-add singular part later
        dn <- dimnames(x)
        x <- x[,p1]
        attr(x, "assign") <- assign[p1] ## needed for splitFrame to work
      }
      # Check if there are factors
      if( control$initial=="SM" ) {
        split <- splitFrame(mf, x, control$split.type)
        if (ncol(split$x1) == 0) {
          control$initial <- 'S'
          warning("No categorical variables found in model. Reverting to an MM-estimator.")
        }
      }
      if( control$initial=="SM" ) {
        z <- SMPY(mf=mf, y=y, control=control, split=split)
      } else if( control$initial == "S" ) {
        z <- MMPY(X=x, y=y, control=control, mf=mf)
      } else stop('Unknown value for lmrobdet.control()$initial')
      # save to complete the DCML object just below
      z.tmp <- list(rank=z$rank, converged=z$converged, qr=z$qr,
                    df.residual=z$df.residual, iter=z$iter, rweights = z$rweights)
      # DCML
      # LS is already computed in z0
      # if MMPY or SMPY return an exact fit then do not run DCML
      if(z$scale > .Machine$double.eps)
        z <- DCML(x=x, y=y, z=z, z0=z0, control=control)
      z$rank <- z.tmp$rank
      z$converged <- z.tmp$converged
      z$qr <- z.tmp$qr
      z$df.residual <- z.tmp$df.residual
      z$iter <- z.tmp$iter
      z$rweightsMM <- z.tmp$rweights
      if(control$compute.rd && !is.null(x))
        z$MD <- robustbase::robMD(x, attr(mt, "intercept"), wqr=z$qr)
      if(model)
        z$model <- mf
      if(ret.x)
        z$x <- if (singular.fit || (!is.null(w) && zero.weights))
          model.matrix(mt, mf, contrasts) else x
      if (ret.y)
        z$y <- if (!is.null(w)) model.response(mf, "numeric") else y
      #
      #       z$coefficients <- z2$coefficients
      #       z$scale <- z2$scale
      #       z$residuals <- z2$residuals
      #       z$cov <- z2$cov
      #       z$fitted.values <- y - z2$residuals
      #       z$rweights <- z$loss <- NULL

      # z <- lmrob.fit(x, y, control, init=init, mf = mf) #-> ./lmrob.MM.R
      if (singular.fit) {
        coef <- numeric(p)
        coef[p2] <- NA
        coef[p1] <- z$coefficients
        names(coef) <- dn[[2L]]
        z$coefficients <- coef
        ## Update QR decomposition (z$qr)
        ## pad qr and qraux with zeroes (columns that were pivoted to the right in z0)
        d.p <- p-rankQR
        n <- NROW(y)
        z$qr[c("qr","qraux","pivot")] <-
          list(matrix(c(z$qr$qr, rep.int(0, d.p*n)), n, p,
                      dimnames = list(dn[[1L]], dn[[2L]][piv])),
               ## qraux:
               c(z$qr$qraux, rep.int(0, d.p)),
               ## pivot:
               piv)
      }
    } else { ## rank 0
      z <- list(coefficients = if (is.matrix(y)) matrix(NA,p,ncol(y))
                else rep.int(as.numeric(NA), p),
                residuals = y, scale = NA, fitted.values = 0 * y,
                cov = matrix(,0,0), rweights = rep.int(as.numeric(NA), NROW(y)),
                weights = w, rank = 0, df.residual = NROW(y),
                converged = TRUE, iter = 0, control=control)
      if (is.matrix(y)) colnames(z$coefficients) <- colnames(x)
      else names(z$coefficients) <- colnames(x)
      if(!is.null(offset)) z$residuals <- y - offset
    }
    if (!is.null(w)) {
      z$residuals <- z$residuals/wts
      z$fitted.values <- save.y - z$residuals
      z$weights <- w
      if (zero.weights) {
        coef <- z$coefficients
        coef[is.na(coef)] <- 0
        f0 <- x0 %*% coef
        if (ny > 1) {
          save.r[ok, ] <- z$residuals
          save.r[nok, ] <- y0 - f0
          save.f[ok, ] <- z$fitted.values
          save.f[nok, ] <- f0
        }
        else {
          save.r[ok] <- z$residuals
          save.r[nok] <- y0 - f0
          save.f[ok] <- z$fitted.values
          save.f[nok] <- f0
        }
        z$residuals <- save.r
        z$fitted.values <- save.f
        z$weights <- save.w
        rw <- z$rweights
        z$rweights <- rep.int(0, length(save.w))
        z$rweights[ok] <- rw
      }
    }
  }
  if(!is.null(offset))
    z$fitted.values <- z$fitted.values + offset

  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- contrasts
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$assign <- assign
  if(control$compute.rd && !is.null(x))
    z$MD <- robustbase::robMD(x, attr(mt, "intercept"), wqr=z$qr)
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- if (singular.fit || (!is.null(w) && zero.weights))
      model.matrix(mt, mf, contrasts) else x
  if (ret.y)
    z$y <- if (!is.null(w)) model.response(mf, "numeric") else y
  # class(z) <- c("lmrobdet", "lmrob")
  # z
  # tmp <- z
  # tmp$MM <- NULL
  # z2 <- list(DCML=tmp, MM=z$MM)
  # class(z2$DCML) <- c("DCML", "lmrob")
  # class(z2$MM) <- "lmrob"
  # class(z2) <- "lmrobdet"
  # z2
  class(z) <- c('lmrobdetDCML', 'lmrobdetMM', 'lmrob')
  z
}


#' Robust estimators for linear regression with fixed designs
#'
#' This function computes a robust regression
#' estimator for a linear models with fixed designs.
#'
#' This function computes robust regression estimators for linear
#' models with fixed designs. It computes an L1 estimator,
#' and uses it as a starting point to find a minimum of a
#' re-descending M estimator. The scale is set to a quantile of the
#' absolute residuals from the L1 estimator.
#' This function makes use of the functions \code{lmrob.fit},
#' \code{lmrob..M..fit}, \code{.vcov.avar1}, \code{lmrob.S} and
#' \code{lmrob.lar}, from robustbase,
#' along with utility functions used by these functions,
#' modified so as to include use of the analytic form of the
#' optimal psi and rho functions (for the optimal psi function , see
#' Section 5.8.1 of Maronna, Martin, Yohai and Salibian Barrera, 2019)
#'
#' @param formula a symbolic description of the model to be fit.
#' @param data an optional data frame, list or environment containing
#' the variables in the model. If not found in \code{data}, model variables
#' are taken from \code{environment(formula)}, which usually is the
#' root environment of the current R session.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param weights an optional vector of weights to be used in the fitting process.
#' @param na.action a function to indicates what should happen when the data contain NAs.
#' The default is set by the \link{na.action} setting of \code{\link[base]{options}}, and is
#' \code{na.fail} if that is unset.
#' @param model logical value indicating whether to return the model frame
#' @param x logical value indicating whether to return the model matrix
#' @param y logical value indicating whether to return the vector of responses
#' @param singular.ok logical value. If \code{FALSE} a singular fit produces an error.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \link{model.matrix.default}.
#' @param offset this can be used to specify an a priori known component to be included
#' in the linear predictor during fitting. An offset term can be included in the formula
#' instead or as well, and if both are specified their sum is used.
#' @param control a list specifying control parameters as returned by the function
#' \link{lmrobM.control}.
#'
#' @return A list with the following components:
#' \item{coefficients}{The estimated vector of regression coefficients}
#' \item{scale}{The estimated scale of the residuals}
#' \item{residuals}{The vector of residuals associated with the robust fit}
#' \item{converged}{Logical value indicating whether IRWLS iterations for the MM-estimator have converged}
#' \item{iter}{Number of IRWLS iterations for the MM-estimator}
#' \item{rweights}{Robustness weights for the MM-estimator}
#' \item{fitted.values}{Fitted values associated with the robust fit}
#' \item{rank}{Numeric rank of the fitted linear model}
#' \item{cov}{The estimated covariance matrix of the regression estimates}
#' \item{df.residual}{The residual degrees of freedom}
#' \item{contrasts}{(only where relevant) the contrasts used}
#' \item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting}
#' \item{call}{the matched call}
#' \item{model}{if requested, the model frame used}
#' \item{x}{if requested, the model matrix used}
#' \item{y}{if requested, the response vector used}
#' \item{na.action}{(where relevant) information returned by model.frame on the special handling of NAs}
#'
#' @author Victor Yohai, \email{vyohai@gmail.com}, based on \code{lmrob}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(shock)
#' cont <- lmrobM.control(bb = 0.5, efficiency = 0.85, family = "bisquare")
#' shockrob <- lmrobM(time ~ n.shocks, data = shock, control=cont)
#' shockrob
#' summary(shockrob)
#'
#' @export
lmrobM <- function(formula, data, subset, weights, na.action,
                   model = TRUE, x = FALSE, y = FALSE,
                   singular.ok = TRUE, contrasts = NULL, offset = NULL,
                   control = lmrobM.control()) {
  # tuning.psi = 3.4434 # 85% efficiency
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset) && length(offset) != NROW(y))
    stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                  length(offset), NROW(y)), domain = NA)
  x <- model.matrix(mt, mf, contrasts)
  contrasts <- attr(x, "contrasts")
  assign <- attr(x, "assign")
  p <- ncol(x)
  if(!is.null(offset))
    y <- y - offset
  if (!is.null(w)) {
    ## checks and code copied/modified from lm.wfit
    ny <- NCOL(y)
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
      ok <- w != 0
      nok <- !ok
      w <- w[ok]
      x0 <- x[nok, , drop = FALSE]
      x  <- x[ ok, , drop = FALSE]
      n <- nrow(x)
      y0 <- if (ny > 1L) y[nok, , drop = FALSE] else y[nok]
      y  <- if (ny > 1L) y[ ok, , drop = FALSE] else y[ok]
      ## add this information to model.frame as well
      ## need it in outlierStats()
      ## ?? could also add this to na.action, then
      ##    naresid() would pad these as well.
      attr(mf, "zero.weights") <- which(nok)
    }
    wts <- sqrt(w)
    save.y <- y
    x <- wts * x
    y <- wts * y
  }
  ## check for singular fit

  if(getRversion() >= "3.1.0") {
    z0 <- .lm.fit(x, y) #, tol = control$solve.tol)
    piv <- z0$pivot
  } else {
    z0 <- lm.fit(x, y) #, tol = control$solve.tol)
    piv <- z0$qr$pivot
  }
  rankQR <- z0$rank

  singular.fit <- rankQR < p
  if (rankQR > 0) {
    if (singular.fit) {
      if (!singular.ok) stop("singular fit encountered")
      pivot <- piv
      p1 <- pivot[seq_len(rankQR)]
      p2 <- pivot[(rankQR+1):p]
      ## to avoid problems in the internal fitting methods,
      ## split into singular and non-singular matrices,
      ## can still re-add singular part later
      dn <- dimnames(x)
      x <- x[,p1]
      attr(x, "assign") <- assign[p1] ## needed for splitFrame to work
    }
    outL <- lmrob.lar(x=x, y=y, control = control, mf = NULL)
    resL <- sort(abs(outL$resid))
    p <- length(outL$coef)
    n <- length(outL$resid)
    L1scale <- resL[(n+p)/2]/.6745
    initial <- list(coefficients=outL$coef, scale=L1scale)
    # create a call to lmrob.control() using the
    # values of the common arguments with the user's call to lmrobdet.control()
    #tmp <- as.call(control)
    #tmp[[1]] <- quote(lmrob.control)
    #my.control <- eval(tmp) # there are no environments to worry about here
    # build the call to lmrob, adding the lmrob control paramaters
    #tmp2 <- cl
    ##tmp2$control <- my.control
    #tmp2[[1]] <- quote(lmrob)
    #tmp2$init <- initial
    #z <- eval(expr=tmp2, envir=parent.frame()) # lmrob(formula, control=my.control, init=initial) #tuning.psi=tuning.psi ,init=initial) lmrob.control(tuning.psi=lmrobdet.control()$tuning.psi),
    rb.ctl <- lmrob.control(psi = control$family, #tuning.psi$name,
                            tuning.psi = control$tuning.psi, #$cc,
                            method = "M", max.it = control$max.it,
                            rel.tol = control$rel.tol,
                            trace.lev = control$trace.lev,
                            cov = ".vcov.w")
    z <- lmrob.fit(x, y, rb.ctl, initial, mf)
    oldz.control <- z$control
    z$control <- control
    z$control$method <- oldz.control$method
    z$scale <- mscale(u=z$resid, tol = control$mscale_tol, delta=control$bb*(1-p/length(z$resid)), tuning.chi=control$tuning.chi, family=control$family, max.it = control$mscale_maxit)
    r.squared <- adj.r.squared <- NA
    if(z$scale > .Machine$double.eps) {
      # if(control$family == 'bisquare') {
      s2 <- mean(rho(z$resid/z$scale, family = control$family, cc=control$tuning.psi))
      if( p != attr(mt, "intercept") ) {
        df.int <- if (attr(mt, "intercept"))
          1L
        else 0L
        if(df.int == 1L) {
          tmp <- as.vector(refine.sm(x=matrix(rep(1,n), n, 1), y=y, initial.beta=median(y),
                                     initial.scale=z$scale, k=500,
                                     conv=1, family = control$family, cc = control$tuning.psi, step='M')$beta.rw)
          s02 <- mean(rho((y-tmp)/z$scale, family = control$family, cc=control$tuning.psi))
        } else {
          s02 <- mean(rho(y/z$scale, family = control$family, cc=control$tuning.psi))
        }
        r.squared <- INVTR2( (s02 - s2)/(s02*(1-s2)), control$family, control$tuning.psi)
        adj.r.squared <- ((n-1)/(n-p))*r.squared -(p-1)/(n-p) # ( s02/(n-1) - s2/(n-z$rank) ) / (s02/(n-1)) # n-p? p.193
      }
      # }
    }
      z$r.squared <- r.squared
      z$adj.r.squared <- adj.r.squared
  } else { ## rank 0
    z <- list(coefficients = if (is.matrix(y)) matrix(NA,p,ncol(y))
              else rep.int(as.numeric(NA), p),
              residuals = y, scale = NA, fitted.values = 0 * y,
              cov = matrix(,0,0), rweights = rep.int(as.numeric(NA), NROW(y)),
              weights = w, rank = 0, df.residual = NROW(y),
              converged = TRUE, iter = 0, control=control)
    if (is.matrix(y)) colnames(z$coefficients) <- colnames(x)
    else names(z$coefficients) <- colnames(x)
    if(!is.null(offset)) z$residuals <- y - offset
  }

  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- contrasts
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$assign <- assign
  # if(control$compute.rd && !is.null(x))
  #   z$MD <- robMD(x, attr(mt, "intercept"), wqr=z$qr)
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- if (singular.fit || (!is.null(w) && zero.weights))
      model.matrix(mt, mf, contrasts) else x
  if (ret.y)
    z$y <- if (!is.null(w)) model.response(mf, "numeric") else y
  class(z) <- c('lmrobM', 'lmrobdetMM')
  z
}



#' Robust likelihood ratio test for linear hypotheses
#'
#' This function computes a robust likelihood ratio test for linear hypotheses.
#'
#' @export rob.linear.test lmrobdetLinTest
#' @aliases rob.linear.test lmrobdetLinTest
#' @rdname rob.linear.test
#'
#' @param object1 an \code{lmrobdetMM} or \code{lmrobM} object with the fit corresponding to the complete model
#' @param object2 an \code{lmrobdetMM} or \code{lmrobM} object with the fit corresponding to the model
#' restricted under the null linear hypothesis.
#'
#' @return A list with the following components: c("test","chisq.pvalue","f.pvalue","df")
#' \item{test}{The value of the F-statistic}
#' \item{f.pvalue}{p-value based on the F distribution}
#' \item{chisq.pvalue}{p-value based on the chi-squared distribution}
#' \item{df}{degrees of freedom}
#'
#' @author Victor Yohai, \email{vyohai@gmail.com}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#'
#' @examples
#' data(oats)
#' cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")
#' oats1M <- lmrobM(response1 ~ variety+block, control=cont, data=oats)
#' oats1M_var <- lmrobM(response1 ~ block, control=cont, data=oats)
#' ( anov1M_var <- rob.linear.test(oats1M, oats1M_var) )
#'
lmrobdetLinTest <- rob.linear.test <- function(object1, object2)
{
  p <- length(object1$coeff)
  q <- length(object1$coeff) - length(object2$coeff)
  n <- length(object1$resid)

  family1 <- object1$control$family
  family2 <- object2$control$family
  cc1 <- object1$control$tuning.psi
  cc2 <- object2$control$tuning.psi

  same.name <- (family1 == family2)
  same.cc <- isTRUE(all.equal(cc1, cc2))

  if(!(same.name && same.cc))
    stop("object1 and object2 do not have the same rho family")

  s <- object1$scale
  a <- sum(rho(object2$resid/s, family=family1, cc=cc1, standardize = TRUE))
  b <- sum(rho(object1$resid/s, family=family1, cc=cc1, standardize = TRUE))
  c <- sum(rhoprime2(object1$resid/s, family=family1, cc=cc1, standardize = TRUE))
  d <- sum(rhoprime(object1$resid/s, family=family1, cc=cc1, standardize = TRUE)^2)
  test <- (2 * (a - b) * c)/d
  # Sanity check: passed
  # a2 <- sum(rho(object2$resid/s, cc=cc))
  # b2 <- sum(rho(object1$resid/s, cc=cc))
  # c2 <- sum(rhoprime2(object1$resid/s, cc=cc))*6/cc
  # d2 <- sum((rhoprime(object1$resid/s, cc=cc)*6/cc)^2)
  # test2 <- (2 * (a2 - b2) * c2)/d2
  # print(c(test, test2))
  chisq.pvalue <- 1- pchisq(test, q)
  F.pvalue <- 1-pf(test/q, q, n - p)
  df <- c(q,n-p)
  output <- list(test=test, chisq.pvalue=chisq.pvalue, F.pvalue=F.pvalue, df=df)
  output
}


#' Tuning parameters for lmrobM
#'
#' This function sets tuning parameters for the M estimators of regression implemented
#' in \code{\link{lmrobM}}.
#'
#' @rdname lmrobM.control
#' @param bb tuning constant (between 0 and 1/2) for the M-scale used to compute the residual
#' scale estimator. Defaults to 0.5.
#' @param tuning.chi tuning constant for the function used to compute the M-scale
#' used for the residual scale estimator. If missing, it is computed inside \code{lmrobdet.control} to match
#' the value of \code{bb} according to the family of rho functions specified in \code{family}.
#' @param tuning.psi tuning parameters for the regression M-estimator computed with a rho function
#' as specified with argument \code{family}. If missing, it is computed inside \code{lmrobdet.control} to match
#' the value of \code{efficiency} according to the family of rho functions specified in \code{family}.
#' Appropriate values for \code{tuning.psi} for a given desired efficiency for Gaussian errors
#' can be constructed using the functions \link{bisquare}, \link{mopt} and \link{opt}.
#' @param efficiency desired asymptotic efficiency of the final regression M-estimator. Defaults to 0.85.
#' @param max.it maximum number of IRWLS iterations for the M-estimator
#' @param mscale_tol Convergence tolerance for the M-scale algorithm. See \code{\link{mscale}}.
#' @param mscale_maxit Maximum number of iterations for the M-scale algorithm. See \code{\link{mscale}}.
#' @param rel.tol relative covergence tolerance for the IRWLS iterations for the M-estimator
#' @param trace.lev positive values (increasingly) provide details on the progress of the M-algorithm
#' @param family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "opt" and "mopt"). Incomplete entries will be matched to
#' the current valid options.
#'
#' @return A list with the necessary tuning parameters.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @examples
#' data(coleman, package='robustbase')
#' m2 <- lmrobM(Y ~ ., data=coleman, control=lmrobM.control())
#' m2
#' summary(m2)
#'
#' @export
lmrobM.control <- function(bb = 0.5,
                           efficiency = 0.99,
                           family = 'opt',
                           tuning.chi, tuning.psi,
                           max.it = 100, rel.tol = 1e-7,
                           mscale_tol = 1e-06,
                           mscale_maxit = 50,
                           trace.lev = 0)
{

  family <- match.arg(family, choices = FAMILY.NAMES)
  if( (efficiency > .9999 ) & ( (family=='mopt') | (family=='opt') ) ) {
    efficiency <- .9999
    warning("Current implementation of \'opt\' or \'mopt\' only allows efficiencies up to 99.99%. Efficiency set to 99.99% for this call.")
  }
  if(missing(tuning.psi))
    tuning.psi <- do.call(family, args=list(e=efficiency))
  if( (length(tuning.psi) == 1) & is.null(names(tuning.psi)) )
    tuning.psi <- c( 'c' = tuning.psi )
  if(missing(tuning.chi))
    tuning.chi <- adjustTuningVectorForBreakdownPoint(family=family, cc=tuning.psi, breakdown.point = bb)
  return(list(family=family, # psi=psi,
              bb=bb, tuning.psi=tuning.psi,
              tuning.chi = tuning.chi,
              max.it=max.it,
              rel.tol=rel.tol, mscale_tol = mscale_tol, mscale_maxit = mscale_maxit,
              trace.lev=trace.lev, mts=1000)) # mts maximum number of subsamples. Un-used, but passed (unnecessarily) to the function that performs M-iterations (lmrob..M..fit), so set here.
}

