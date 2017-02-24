\name{lmrob..M..fit}
\alias{lmrob..M..fit}
\title{Compute M-estimators of regression}
\description{
  This function performs RWLS iterations to find an
  M-estimator of regression.  When started from an S-estimated
  \code{beta.initial}, this results in an MM-estimator.
}
\usage{
lmrob..M..fit(x, y, beta.initial, scale, control, obj,
              mf = obj$model, method = obj$control$method)
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n x p}) typically including a
    column of \code{1}s for the intercept.}
  \item{y}{numeric response vector (of length \eqn{n}).}
  \item{beta.initial}{numeric vector (of length \eqn{p}) of initial
    estimate.  Usually the result of an S-regression estimator.}
  \item{scale}{robust residual scale estimate. Usually
    an S-scale estimator.}
  \item{control}{list of control parameters, as returned
    by \code{\link{lmrob.control}}.  Currently, the components
    \code{c("max.it", "rel.tol","trace.lev", "psi", "tuning.psi", "mts", "subsampling")}
    are accessed.}
  \item{obj}{an optional \code{lmrob}-object.  If specified, this is
    typically used to set values for the other arguments.}
  \item{mf}{(optional) a model frame as returned by \code{\link{model.frame}},
    used only to compute outlier statistics, see \code{\link{outlierStats}}.}
  \item{method}{optional; the \code{method} used for \emph{obj} computation.}
}
\details{
  This function is used by \code{\link{lmrob.fit}} (and
  \code{anova(<lmrob>, type = "Deviance")}) and typically not to be used
  on its own.
}
\value{A list with the following elements:
  \item{coef}{the M-estimator (or MM-estim.) of regression}
  \item{control}{the \code{control} list input used}
  \item{scale}{ The residual scale estimate}
  \item{seed }{ The random number generator seed}
  \item{converged}{ \code{TRUE} if the RWLS iterations converged,
    \code{FALSE} otherwise}
}
\references{
  Yohai, 1987
}
\seealso{
  \code{\link{lmrob.fit}}, \code{\link{lmrob}};
  \code{\link[MASS]{rlm}} from package \pkg{MASS}.
}
\author{Matias Salibian-Barrera and Martin Maechler}
\examples{
data(stackloss)
X <- model.matrix(stack.loss ~ . , data = stackloss)
y <- stack.loss
## Compute manual MM-estimate:
## 1) initial LTS:
m0 <- ltsReg(X[,-1], y)
## 2) M-estimate started from LTS:
m1 <- lmrob..M..fit(X, y, beta.initial = coef(m0), scale = m0$scale, method = "SM",
                    control = lmrob.control(tuning.psi = 1.6, psi = 'bisquare'))
## no 'method' (nor 'obj'):
m1. <- lmrob..M..fit(X, y, beta.initial = coef(m0), scale = m0$scale,
                     control = m1$control)
stopifnot(all.equal(m1, m1., tol = 1e-15)) # identical {call *not* stored!}

cbind(m0$coef, m1$coef)
## the scale is kept fixed:
stopifnot(identical(unname(m0$scale), m1$scale))

##  robustness weights: are
r.s <- with(m1, residuals/scale) # scaled residuals
m1.wts <- Mpsi(r.s, cc = 1.6, psi="tukey") / r.s
summarizeRobWeights(m1.wts)
##--> outliers 1,3,4,13,21
which(m0$lts.wt == 0) # 1,3,4,21 but not 13
\dontshow{stopifnot(which(m0$lts.wt == 0) == c(1,3,4,21))
}
## Manually add M-step to SMD-estimate (=> equivalent to "SMDM"):
m2 <- lmrob(stack.loss ~ ., data = stackloss, method = 'SMD')
m3 <- lmrob..M..fit(obj = m2)

## Simple function that allows custom initial estimates
## (Deprecated; use init argument to lmrob() instead.) %% MM: why deprecated?
lmrob.custom <- function(x, y, beta.initial, scale, terms) {
  ## initialize object
  obj <- list(control = lmrob.control("KS2011"),
              terms = terms) ## terms is needed for summary()
  ## M-step
  obj <- lmrob..M..fit(x, y, beta.initial, scale, obj = obj)
  ## D-step
  obj <- lmrob..D..fit(obj, x)
  ## Add some missing elements
  obj$cov <- TRUE ## enables calculation of cov matrix
  obj$p <- obj$qr$rank
  obj$degree.freedom <- length(y) - obj$p
  ## M-step
  obj <- lmrob..M..fit(x, y, obj=obj)
  obj$control$method <- ".MDM"
  obj
}

m4 <- lmrob.custom(X, y, m2$init$init.S$coef,
                   m2$init$scale, m2$terms)
stopifnot(all.equal(m4$coef, m3$coef))

## Start from ltsReg:
m5 <- ltsReg(stack.loss ~ ., data = stackloss)
m6 <- lmrob.custom(m5$X, m5$Y, coef(m5), m5$scale, m5$terms)
}
\keyword{robust}
\keyword{regression}
