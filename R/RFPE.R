#' Robust Final Prediction Error
#'
#' This function computes the robust Final Prediction Errors (RFPE) for a robust regression fit using M-estimates.
#' It is used internally by \code{\link{step.lmrobdetMM}} and not meant to be used
#' directly.
#'
#' @param object the \code{MM} element (of class \code{\link{lmrob}}) in an object of class \code{\link{lmrobdetMM}}.
#' @param scale a numeric value specifying the scale estimate used to compute the RFPE. Usually this 
#' should be the scale estimate from an encompassing model. If \code{NULL}, the scale estimate in 
#' \code{object} is used.
#'
#' @return the robust final prediction error (numeric).
#'
#' @rdname lmrobdetMM.RFPE
#' @author Victor Yohai, \email{victoryohai@gmail.com}, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{lmrobdetMM}}
#'
#' @export
lmrobdetMM.RFPE <- function(object, scale = NULL)
{
  if (!object$converged)
    warning("The algorithm did not converge, inference is not recommended.")
  # ocm <- tolower(object$control$method)
  # nocm <- nchar(ocm)
  # if(substr(ocm, nocm, nocm) != "m")
  #   stop("RFPE is only available for MM-estimates.")
  # if ( (casefold(object$control$method) != "sm") ) { # & (casefold(object$control$method) != "m-sm") )
  #   stop("RFPE is only available for MM-estimates.")
  p <- length(object$coef)
  if (is.null(scale))
    scale <- object$scale
  res <- residuals(object)/scale
  # psif <- object$control$psi
  # tun <- object$control$tuning.psi
  # efficiency <- object$robust.control$efficiency
  # if (casefold(psif[2]) == "optimal")
  #   ipsi <- 1
  # else ipsi <- 2
  # yc <- object$yc
  # wf <- .Mwgt.psi1(psi=psif, cc=tun)
  # a <- sum(Mpsi(x=res, cc=tun, psi=psif, deriv=-1))
  # b <- p * sum(Mpsi(x=res, cc=tun, psi=psif, deriv=0)^2)
  # d <- sum(Mpsi(x=res, cc=tun, psi=psif, deriv=1))
  # a <- mean(rho(u=res, family=object$control$tuning.psi))
  # b <- p * mean(rhoprime(u=res, family=object$control$tuning.psi)^2)
  # d <- mean(rhoprime2(u=res, family=object$control$tuning.psi))
  # tun <- object$control$tuning.psi$cc
  a2 <- mean(rho(u=res, family = object$control$family, cc = object$control$tuning.psi, standardize=TRUE))
  b2 <- p * mean(rhoprime(u=res, family = object$control$family, cc = object$control$tuning.psi, standardize=TRUE)^2)
  d2 <- mean(rhoprime2(u=res, family = object$control$family, cc = object$control$tuning.psi, standardize=TRUE))
  if (d2 <= 0)
    return(NA)
  return( (a2 + b2/d2/length(res)) ) # (a + b/d)*6 / tun^2 )
}



#' RFPE of submodels of an \code{\link{lmrobdetMM}} fit
#'
#' This function computes the RFPE for the MM-estimators obtained with \code{\link{lmrobdetMM}} by
#' recomputing it, successively removing each of a number of specified terms. 
#' It is used internally by \code{\link{step.lmrobdetMM}} and not meant to be used
#' directly.
#'
#' @param object the \code{MM} element (of class \code{\link{lmrob}}) in an object of class \code{\link{lmrobdetMM}}.
#' @param scope an optional \code{formula} giving the terms to be considered for dropping. Typically 
#' this argument is omitted, in which case all possible terms are dropped (without breaking hierarchy 
#' rules). The \code{scope} can also be a character vector of term labels. If the argument is supplied as a 
#' formula, any \code{.} is interpreted relative to the formula implied by the \code{object} argument.
#' @param scale an optional residual scale estimate. If missing the residual
#' scale estimate in \code{object} is used.
#' @param keep a character vector of names of components that should be saved for each subset model. 
#' Only names from the set \code{"coefficients"}, \code{"fitted"} and \code{"residuals"}
#' are allowed. If \code{keep == TRUE}, the complete set is saved. The default behavior is 
#' not to keep anything.
#' @param \dots additional parameters to match generic method \code{drop1}
#'
#' @return An anova object consisting of the term labels, the degrees of freedom, and Robust Final 
#' Prediction Errors (RFPE) for each subset model. If \code{keep} is missing, the anova object is 
#' returned. If \code{keep} is present, a list with components \code{"anova"} and \code{"keep"} is returned. 
#' In this case, the \code{"keep"} component is a matrix of mode \code{"list"}, with a column for each 
#' subset model, and a row for each component kept.
#'
#' @rdname drop1.lmrobdetMM
#' @author Victor Yohai, \email{victoryohai@gmail.com},  Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{lmrobdetMM}}
#'
#' @export
drop1.lmrobdetMM <- function (object, scope, scale, keep, ...)
{
  # if ( (casefold(object$control$method) != "sm") ) # & (casefold(object$control$method) != "m-sm") )
  #   stop("drop1 is only available for MM-estimates.")
  # ocm <- tolower(object$control$method) #MM$control$method)
  # nocm <- nchar(ocm)
  # if(substr(ocm, nocm, nocm) != "m")
  #   stop("RFPE is only available for MM-estimates.")
  if (!object$converged) #MM$converged)
    stop("The algorithm did not converge, inference is not recommended.")
  x <- model.matrix(object) #$MM)
  asgn <- attr(x, "assign")
  term.labels <- attr(object$terms, "term.labels") #MM$terms
  dfs <- table(asgn[asgn > 0])
  names(dfs) <- term.labels
  # psif <- object$control$psi # psif <- object$robust.control$weight
  if (missing(scope))
    scope <- drop.scope(object)
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, term.labels, FALSE)))
      stop("scope is not a subset of term labels")
  }
  dfs <- dfs[scope]
  k <- length(scope)
  if (missing(scale))
    scale <- object$scale #MM$scale
  if (!missing(keep)) {
    max.keep <- c("coefficients", "fitted", "residuals")
    if (is.logical(keep) && keep)
      keep <- max.keep
    else {
      if (!all(match(keep, max.keep, FALSE)))
        stop(paste("Can only keep one or more of: \"",
                   paste(max.keep, collapse = "\", \""), "\"",
                   sep = ""))
    }
    value <- array(vector("list", 3 * k), c(k, 3), list(scope, c("coefficients", "fitted", "residuals")))
  }
  else keep <- character(0)
  rfpe <- vector('numeric', k)
  for (i in 1:k) {
      curfrm <- as.formula(paste(".~.-", scope[[i]]))
      curobj <- update(object, curfrm)
      rfpe[i] <- lmrobdetMM.RFPE(curobj, scale) #$MM, scale)
      if (length(keep)) {
        value[i, 1] <- list(curobj$coefficients) #MM$coefficients)
        value[i, 2] <- list(curobj$fitted) #MM$fitted)
        value[i, 3] <- list(curobj$residuals) #MM$residuals)
      }
    }
  scope <- c("<none>", scope)
  dfs <- c(0, dfs)
  rfpe <- c(lmrobdetMM.RFPE(object, scale), rfpe) #MM, scale), rfpe)
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, RFPE = rfpe, row.names = scope,
                    check.names = FALSE)
  head <- c("\nSingle term deletions", "\nModel:", deparse(as.vector(formula(object)))) #$MM))))
  if (!missing(scale))
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  oldClass(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  if (length(keep)) {
    value <- value[, keep, drop = FALSE]
    oldClass(value) <- "matrix"
    list(anova = aod, keep = value)
  }
  else aod
}

#' Robust stepwise using RFPE
#'
#' This function performs stepwise model selection on a robustly fitted
#' linear model using the RFPE
#' criterion and the robust regression estimators computed with
#' \code{\link{lmrobdetMM}}. Only backwards stepwise is currently implemented.
#'
#' Presently only backward stepwise selection is supported. During each step the 
#' Robust Final Prediction Error (as computed by the function \code{lmrobdetMM.RFPE}) is 
#' calculated for the current model and for each sub-model achievable by deleting a 
#' single term. If the argument \code{whole.path} is \code{FALSE}, the function steps 
#' to the sub-model with the lowest 
#' Robust Final Prediction Error or, if the current model has the lowest Robust Final 
#' Prediction Error, terminates. If the argument \code{whole.path} is \code{TRUE}, the 
#' function steps through all smaller submodels removing, at each step, the variable 
#' that most reduces the Robust Final Prediction Error. The scale estimate from \code{object} 
#' is used to compute the Robust Final Prediction Error throughout the procedure.
#'
#' @param object a robust fit as returned by \code{\link{lmrobdetMM}}
#' @param scope either a formula or a list with elements \code{lower} and \code{upper} each of 
#' which is a formula. The terms in the right-hand-side of \code{lower} are always included 
#' in the model and the additional terms in the right-hand-side of \code{upper} are the
#' candidates for inclusion/exclusion from the model. If a single formula is given, it is 
#' taken to be \code{upper}, and \code{lower} is set to the empty model. The \code{.} operator 
#' is interpreted in the context of the formula in \code{object}.
#' @param direction the direction of stepwise search. Currenly only \code{backward} stepwise 
#' searches are implemented.
#' @param trace logical. If \code{TRUE} information about each step is printed on the screen.
#' @param keep a filter function whose input is a fitted model object and the associated AIC statistic, and whose output is arbitrary. Typically keep will select a subset of the components of the object and return them. The default is not to keep anything.
#' @param steps maximum number of steps to be performed. Defaults to 1000, which should mean as many as needed.
#' @param whole.path if \code{FALSE} (default) variables are dropped until the RFPE fails to improve. If \code{TRUE} the best variable to be dropped is removed, even if this does not improve the RFPE.
#'
#' @return If \code{whole.path == FALSE} the function returns the robust fit as obtained by \code{lmrobdetMM} using the final model. 
#' If \code{whole.path == TRUE} a list is returned containing the RFPE of each model on the sequence
#' of submodels. The names of the components of this list are the formulas that correspods to each model. 
#'
#' @aliases step.lmrobdetMM step.lmrobdet
#' @rdname step.lmrobdetMM
#' @author Victor Yohai, \email{victoryohai@gmail.com}, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @examples
#' cont <- lmrobdet.control(bb = 0.5, efficiency = 0.85, family = "bisquare")
#' set.seed(300)
#' X <- matrix(rnorm(50*6), 50, 6)
#' beta <- c(1,1,1,0,0,0)
#' y <- as.vector(X %*% beta) + 1 + rnorm(50)
#' y[1:6] <- seq(30, 55, 5)
#' for (i in 1:6) X[i,] <- c(X[i,1:3],i/2,i/2,i/2)
#' Z <- cbind(y,X)
#' Z <- as.data.frame(Z)
#' obj <- lmrobdetMM(y ~ ., data=Z, control=cont)
#' out <- step.lmrobdetMM(obj)
#'
#' @export step.lmrobdetMM step.lmrobdet
step.lmrobdet <- step.lmrobdetMM <- function (object, scope, direction = c("both", "backward", "forward"), trace = TRUE,
                        keep = NULL, steps = 1000, whole.path=FALSE)
{
  # object.MM <- object$MM
  if (missing(direction))
    direction <- "backward"
  else direction <- match.arg(direction)
  if (direction != "backward")
    stop("Presently step.lmrobdetMM only supports backward model selection.")
  if(whole.path) keep <- function(a, b) a
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, namc))
  }
  make.step <- function(models, fit, scale, object) {
    change <- sapply(models, "[[", "change")
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    RFPE <- sapply(models, "[[", "RFPE")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                 "\nInitial Model:", deparse(as.vector(formula(object))),
                 "\nFinal Model:", deparse(as.vector(formula(fit))),
                 "\n")
    aod <- data.frame(Step = change, Df = ddf, `Resid. Df` = rdf,
                      RFPE = RFPE, check.names = FALSE)
    attr(aod, "heading") <- heading
    oldClass(aod) <- c("anova", "data.frame")
    fit$anova <- aod
    fit
  }
  backward <- direction == "both" || direction == "backward"
  forward <- direction == "both" || direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric(0)
    fadd <- NULL
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower))
        attr(terms(update.formula(object, fdrop)), "factor")
      else numeric(0)
      fadd <- if (!is.null(fadd <- scope$upper))
        attr(terms(update.formula(object, fadd)), "factor")
    }
    else {
      fadd <- if (!is.null(fadd <- scope))
        attr(terms(update.formula(object, scope)), "factor")
      fdrop <- numeric(0)
    }
  }
  if (is.null(fadd)) {
    backward <- TRUE
    forward <- FALSE
  }
  m <- model.frame(object) #object$MM)
  obconts <- object$contrasts #$MM$contrasts
  objectcall <- object$call #MM$call
  control <- object$control #MM$control
  # if (forward) {
  #   add.rhs <- paste(dimnames(fadd)[[2]], collapse = "+")
  #   add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
  #   new.form <- update.formula(object, add.rhs, evaluate = FALSE)
  #   fc <- objectcall
  #   Terms <- terms(new.form)
  #   fc$formula <- Terms
  #   fobject <- list(call = fc)
  #   oldClass(fobject) <- oldClass(object)
  #   m <- model.frame(fobject)
  #   x <- model.matrix(Terms, m, contrasts = obconts)
  # }
  # else {
    Terms <- object$terms
    x <- model.matrix(Terms, m, contrasts = obconts)
  # }
  Asgn <- attr(x, "assign")
  term.labels <- attr(Terms, "term.labels")
  a <- attributes(m)
  y <- model.extract(m, "response")
  w <- model.extract(m, "weights")
  if (is.null(w))
    w <- rep(1, nrow(m))
  models <- vector("list", steps)
  if (!is.null(keep)) {
    keep.list <- vector("list", steps)
    nv <- 1
  }
  n <- length(object$fitted) #MM$fitted)
  scale <- object$scale #MM$scale
  fit <- object
  bRFPE <- lmrobdetMM.RFPE(fit) #$MM)
  nm <- 1
  Terms <- fit$terms
  if (trace)
    cat("Start:  RFPE=", format(round(bRFPE, 4)), "\n", deparse(as.vector(formula(fit))),
        "\n\n")
  models[[nm]] <- list(df.resid = fit$df.resid, change = "",
                       RFPE = bRFPE)
  if (!is.null(keep))
    keep.list[[nm]] <- keep(fit, bRFPE)
  RFPE <- bRFPE + 1
  while (bRFPE < RFPE & steps > 0) {
    steps <- steps - 1
    if(!whole.path) RFPE <- bRFPE
    bfit <- fit
    ffac <- attr(Terms, "factor")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && (ndrop <- length(scope$drop))) {
      aod <- drop1.lmrobdetMM(fit, scope$drop, scale)
      if (trace)
        print(aod)
      change <- rep("-", ndrop + 1)
    }
    # if (forward && (nadd <- length(scope$add))) {
    #   aodf <- add1.lmrobdet(fit, scope$add, scale, x = x)
    #   if (trace)
    #     print(aodf)
    #   change <- c(change, rep("+", nadd + 1))
    #   if (is.null(aod))
    #     aod <- aodf
    #   else {
    #     ncaod <- dim(aod)[1]
    #     aod[seq(ncaod + 1, ncaod + nadd + 1), ] <- aodf
    #   }
    # }
    if (is.null(aod))
      break
    o <- (oo <- order(aod[, "RFPE"]))[1]
    if(o[1] == 1) {
      if(!whole.path) break
      o[1] <- oo[2]
    }
    # if (!whole.path & (o[1] == 1) ) {
    #   break
    # }
    change <- paste(change[o], dimnames(aod)[[1]][o])
    Terms <- terms(update(formula(fit), eval(parse(text = paste("~ .", change)))))
    attr(Terms, "formula") <- new.formula <- formula(Terms)
    # control$method <- 'MM'
    newfit <- lmrobdetMM(new.formula, data = m, control = control) #, init=object$init$control$method)
    bRFPE <- aod[, "RFPE"][o]
    if (trace)
      cat("\nStep:  RFPE =", format(round(bRFPE, 4)), "\n",
          deparse(as.vector(formula(Terms))), "\n\n")
    if(!whole.path & bRFPE >= RFPE)
      break
    nm <- nm + 1
    models[[nm]] <- list(df.resid = newfit$df.resid, change = change,
                         RFPE = bRFPE)
    fit <- c(newfit, list(formula = new.formula))
    oc <- objectcall
    oc$formula <- as.vector(fit$formula)
    fit$call <- oc
    oldClass(fit) <- oldClass(object)
    if (!is.null(keep))
      keep.list[[nm]] <- keep(fit, bRFPE)
    if(whole.path) RFPE <- bRFPE + 1
  }
  if (!is.null(keep))
    fit$keep <- keep.list[seq(nm)] # re.arrange(keep.list[seq(nm)])
  if(!whole.path) {
    return(make.step(models = models[seq(nm)], fit, scale, object))
  } else {
    nms <- lapply(fit$keep, function(a) a$call$formula)
    RFPEs <- sapply(models, "[[", "RFPE")
    names(RFPEs) <- nms
    return(RFPEs[seq(nm)])
  }
}


