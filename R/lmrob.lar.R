lmrob.lar <- function(x, y, control = lmrob.control(), mf = NULL)
{
  ## LAR : Least Absolute Residuals -- i.e. L_1  M-estimate
  ## this function is identical to lmRob.lar of the robust package

  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(p > 0, n >= p, length(y) == n, is.numeric(control$rel.tol))
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  bet0 <- 0.773372647623  ## bet0 = pnorm(0.75)
  tmpn <- double(n)
  tmpp <- double(p)

  z1 <- .Fortran(rslarsbi, ##-> ../src/rllarsbi.f
                 x,
                 y,
                 as.integer(n),
                 as.integer(p),
                 as.integer(n),
                 as.integer(n),
                 as.double(control$rel.tol),
                 NIT=integer(1),
                 K=integer(1),
                 KODE=integer(1),
                 SIGMA=double(1),
                 THETA=tmpn,
                 RS=tmpn,
                 SC1=tmpn,
                 SC2=tmpp,
                 SC3=tmpp,
                 SC4=tmpp,
                 BET0=as.double(bet0))[c("THETA","SIGMA","RS","NIT","KODE")]
  if (z1[5] > 1)
      stop("calculations stopped prematurely in rllarsbi\n",
           "(probably because of rounding errors).")
  names(z1) <- c("coefficients", "scale", "residuals", "iter", "status")
  ##           c("THETA",        "SIGMA", "RS",        "NIT",  "KODE")
  z1$converged <- TRUE
  length(z1$coefficients) <- p
  z1
}

splitFrame <- function(mf, x = model.matrix(mt, mf),
                       type = c("f","fi", "fii"))
{
    mt <- attr(mf, "terms")
    type <- match.arg(type)
    x <- as.matrix(x)
    p <- ncol(x)

    ## --- split categorical and interactions of categorical vars.
    ##     from continuous variables
    factors <- attr(mt, "factors")
    factor.idx <- attr(mt, "dataClasses") == "factor"
    if (!any(factor.idx)) ## There are no factors
        return(list(x1.idx = rep.int(FALSE, p), x1=matrix(,nrow(x),0L), x2=x))
    switch(type,
           ## --- include interactions cat * cont in x1:
           fi = { factor.asgn <- which(factor.idx %*% factors > 0) },
           ## --- include also continuous variables that interact with factors in x1:
           ##     make sure to include interactions of continuous variables
           ##     interacting with categorical variables, too
           fii = { factor.asgn <- numeric(0)
                   factors.cat <- factors
                   factors.cat[factors.cat > 0] <- 1L ## fix triple+ interactions
                   factors.cat[, factor.idx %*% factors == 0] <- 0L
                   for (i in 1:ncol(factors)) {
                       comp <- factors[,i] > 0
                       ## if any of the components is a factor: include in x1 and continue
                       if (any(factor.idx[comp])) {
                           factor.asgn <- c(factor.asgn, i)
                       } else {
                           ## if there is an interaction of this term with a categorical var.
                           tmp <- colSums(factors.cat[comp,,drop=FALSE]) >= sum(comp)
                           if (any(tmp)) {
                               ## if no other continuous variables are involved
                               ## include in x1 and continue
                               ## if (identical(factors[!comp, tmp], factors.cat[!comp, tmp]))
                               if (!all(colSums(factors[!factor.idx & !comp, tmp, drop=FALSE]) > 0))
                                   factor.asgn <- c(factor.asgn, i)
                           }
                       }
                   } },
           ## --- do not include interactions cat * cont in x1:
           f = { factor.asgn <- which(factor.idx %*% factors & !(!factor.idx) %*% factors) },
           stop("unknown split type"))
    x1.idx <- attr(x, "assign") %in% c(0, factor.asgn) ## also include intercept
    names(x1.idx) <- colnames(x)

    ## x1: factors and (depending on type) interactions of / with factors
    ## x2: continuous variables
    list(x1 = x[,  x1.idx, drop=FALSE],
         x1.idx = x1.idx,
         x2 = x[, !x1.idx, drop=FALSE])
}
