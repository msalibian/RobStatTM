#' RFPE of a model fit with \code{\link{lmrobdet}}
#'
#' This function computes the RFPE for the MM-estimator obtained using \code{\link{lmrobdet}}.
#'
#' @param object the \code{MM} element (of class \code{\link{lmrob}}) in an object of class \code{\link{lmrobdet}}.
#' @param scale an optional residual scale estimator. If \code{NULL} the residual
#' scale estimator that corresponds to the model in \code{object} is used.
#'
#' @return the RFPE corresponding to the fit model.
#'
#' @rdname lmrobdet.RFPE
#' @author Victor Yohai, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://thebook}
#' @seealso \code{\link{lmrobdet}}
#'
#' @export
lmrobdet.RFPE <- function (object, scale = NULL)
{
  if (!object$converged)
    warning("The algorithm did not converge, inference is not recommended.")
  ocm <- tolower(object$control$method)
  nocm <- nchar(ocm)
  if(substr(ocm, nocm, nocm) != "m")
    stop("RFPE is only available for MM-estimates.")
  # if ( (casefold(object$control$method) != "sm") ) { # & (casefold(object$control$method) != "m-sm") )
  #   stop("RFPE is only available for MM-estimates.")
  p <- length(object$coef)
  if (is.null(scale))
    scale <- object$scale
  res <- residuals(object)/scale
  psif <- object$control$psi
  tun <- object$control$tuning.psi
  # efficiency <- object$robust.control$efficiency
  # if (casefold(psif[2]) == "optimal")
  #   ipsi <- 1
  # else ipsi <- 2
  # yc <- object$yc
  # wf <- .Mwgt.psi1(psi=psif, cc=tun)
  a <- sum(Mpsi(x=res, cc=tun, psi=psif, deriv=-1))
  b <- p * sum(Mpsi(x=res, cc=tun, psi=psif, deriv=0)^2)
  d <- sum(Mpsi(x=res, cc=tun, psi=psif, deriv=1))
  if (d <= 0)
    return(NA)
  (a + b/d)*6 / tun^2
}



#' RFPE of submodels of an \code{\link{lmrobdet}} fit
#'
#' This function computes the RFPE for the MM-estimators obtained with \code{\link{lmrobdet}} by dropping each single variable in the model.
#'
#' @param object the \code{MM} element (of class \code{\link{lmrob}}) in an object of class \code{\link{lmrobdet}}.
#' @param scope a formula giving the terms to be considered for dropping.
#' @param scale an optional residual scale estimator. If \code{NULL} the residual
#' scale estimator that corresponds to the model in \code{object} is used.
#' @param keep a character vector of names of components that should be saved for each subset model. 
#' Only names from the set \code{"coefficients"}, \code{"fitted"} and \code{"residuals"}
#' are allowed. If \code{keep == TRUE}, the complete set is saved. The default behavior is 
#' not to keep anything.
#'
#' @return An anova object consisting of the term labels, the degrees of freedom, and Robust Final 
#' Prediction Errors (RFPE) for each subset model. If \code{keep} is missing, the anova object is 
#' returned. If \code{keep} is present, a list with components \code{"anova"} and \code{"keep"} is returned. 
#' In this case, the \code{"keep"} component is a matrix of mode \code{"list"}, with a column for each 
#' subset model, and a row for each component kept.
#'
#' @rdname drop1.lmrobdet
#' @author Victor Yohai, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://thebook}
#' @seealso \code{\link{lmrobdet}}
#'
#' @export
drop1.lmrobdet <- function (object, scope, scale, keep)
{
  # if ( (casefold(object$control$method) != "sm") ) # & (casefold(object$control$method) != "m-sm") )
  #   stop("drop1 is only available for MM-estimates.")
  ocm <- tolower(object$MM$control$method)
  nocm <- nchar(ocm)
  if(substr(ocm, nocm, nocm) != "m")
    stop("RFPE is only available for MM-estimates.")
  if (!object$MM$converged)
    stop("The algorithm did not converge, inference is not recommended.")
  x <- model.matrix(object$MM)
  asgn <- attr(x, "assign")
  term.labels <- attr(object$MM$terms, "term.labels")
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
    scale <- object$MM$scale
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
      rfpe[i] <- lmrobdet.RFPE(curobj$MM, scale)
      if (length(keep)) {
        value[i, 1] <- list(curobj$MM$coefficients)
        value[i, 2] <- list(curobj$MM$fitted)
        value[i, 3] <- list(curobj$MM$residuals)
      }
    }
  scope <- c("<none>", scope)
  dfs <- c(0, dfs)
  rfpe <- c(lmrobdet.RFPE(object$MM, scale), rfpe)
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, RFPE = rfpe, row.names = scope,
                    check.names = FALSE)
  head <- c("\nSingle term deletions", "\nModel:", deparse(as.vector(formula(object$MM))))
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
#' This function performs stepwise variable selection using the RFPE
#' criterion and the robust regression estimators computed with
#' \code{\link{lmrobdet}}.
#'
#' @param object a robust fit as returned by \code{\link{lmrobdet}}
#' @param scope defines the range of models to be examined. This should be either a single formula, or a list containing components upper and lower, both formulae. See the details below.
#' @param direction the direction of stepwise search. Currenly only \code{backward} is implemented.
#' @param trace if \code{TRUE} information about each step is printed on the screen.
#' @param keep a filter function whose input is a fitted model object and the associated AIC statistic, and whose output is arbitrary. Typically keep will select a subset of the components of the object and return them. The default is not to keep anything.
#' @param steps maximum number of steps to be performed. Defaults to 1000, which should mean as many as needed.
#' @param whole.path if \code{FALSE} (default) variables are dropped until the RFPE fails to improve. If \code{TRUE} the best variable to be dropped is removed, even if this does not improve the RFPE.
#'
#' @return either a robust fit as obtained by \code{lmrobdet} using the final model, or a list of fits one for each step in the process.
#'
#' @rdname step.lmrob
#' @author Victor Yohai, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://thebook}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @export
step.lmrobdet <- function (object, scope, direction = c("both", "backward", "forward"), trace = TRUE,
                        keep = NULL, steps = 1000, whole.path=FALSE)
{
  # object.MM <- object$MM
  if (missing(direction))
    direction <- "backward"
  else direction <- match.arg(direction)
  if (direction != "backward")
    stop("Presently step.lmrobdet only supports backward model selection.")
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
  m <- model.frame(object$MM)
  obconts <- object$MM$contrasts
  objectcall <- object$MM$call
  control <- object$MM$control
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
  n <- length(object$MM$fitted)
  scale <- object$MM$scale
  fit <- object
  bRFPE <- lmrobdet.RFPE(fit$MM)
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
      aod <- drop1.lmrobdet(fit, scope$drop, scale)
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
    newfit <- lmrobdet(new.formula, data = m, control = control) #, init=object$init$control$method)
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


