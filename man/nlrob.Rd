\name{nlrob}
\title{Robust Fitting of Nonlinear Regression Models}
\alias{nlrob}
\alias{fitted.nlrob}
\alias{residuals.nlrob}
\alias{predict.nlrob}
\alias{vcov.nlrob}
\description{
  \code{nlrob} fits a nonlinear regression model by robust methods.
  Per default, by an M-estimator, using iterated reweighted least
  squares (called \dQuote{IRLS} or also \dQuote{IWLS}).
}
\usage{
nlrob(formula, data, start, lower, upper,
      weights = NULL, na.action = na.fail,
      method = c("M", "MM", "tau", "CM", "mtl"),
      psi = .Mwgt.psi1("huber", cc=1.345), scale = NULL,
      test.vec = c("resid", "coef", "w"), maxit = 20,
      tol = 1e-06, acc, algorithm = "default", doCov = FALSE,
      control = if(method == "M") nls.control() else
		nlrob.control(method, optArgs = list(trace=trace), ...),
      trace = FALSE, ...)

\method{fitted}{nlrob}(object, ...)
\method{residuals}{nlrob}(object, type = , ...)% FIXME: more 'type's + DOCU
\method{predict}{nlrob}(object, newdata, ...)
}
\arguments{
  \item{formula}{a nonlinear \code{\link{formula}} including variables
    and parameters of the model, such as \code{y ~ f(x, theta)} (cf. \code{\link{nls}}).
    (For some checks: if \eqn{f(.)} is linear, then we need
    parentheses, e.g., \code{y ~ (a + b * x)};
    (note that \code{._nlrob.w} is not allowed as variable or parameter name))
    %% FIXME in code -- long overdue, as nls() is more flexible *SINCE* R 2.2.1
    %% Do not use \code{w} as variable or parameter name!
    %% FIXME: this should really no longer be needed ==> add a check
  }
  \item{data}{an optional data frame containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from
    which \code{nlrob} is called.}
  \item{start}{a named numeric vector of starting parameters estimates,
    only for \code{method = "M"}.}
  \item{lower, upper}{numeric vectors of lower and upper bounds; if
    needed, will be replicated to be as long as the longest of \code{start},
    \code{lower} or \code{upper}.  For (the default) \code{method = "M"},
    if the bounds are unspecified all parameters are assumed to be
    unconstrained; also, for method \code{"M"}, bounds can only be used
    with the \code{"port"} algorithm.  They are ignored, with a warning,
    in cases they have no effect.

    For all other methods, currently these bounds \emph{must} be
    specified as finite values, and one of them must have
    \code{\link{names}} matching the parameter names in \code{formula}.

    For methods \code{"CM"} and \code{"mtl"}, the bounds must
    \emph{additionally} have an entry named \code{"sigma"} as that is
    determined simultaneously in the same optimization, and hence its
    \code{lower} bound must not be negative.
  }
  \item{weights}{an optional vector of weights to be used in the fitting
    process (for intrinsic weights, not the weights \code{w} used in the
    iterative (robust) fit). I.e.,
    \code{sum(w * e^2)} is minimized with \code{e} = residuals,
    \eqn{e_i = y_i - f(xreg_i, \theta)}{e[i] = y[i] - f(xreg[i], theta)},
    where \eqn{f(x,\theta)}{f(x, theta)} is the nonlinear function,
    and \code{w} are the robust weights from \code{resid * weights}.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.  The default action is for the procedure to
    fail.  If NAs are present, use \code{na.exclude} to have residuals with
    \code{length == nrow(data) == length(w)}, where \code{w} are the
    weights used in the iterative robust loop.  This is better if the
    explanatory variables in
    \code{formula} are time series (and so the NA location is important).
    For this reason, \code{na.omit}, which leads to omission of cases
    with missing values on any required variable, is not suitable
    here since the residuals length is different from
    \code{nrow(data) == length(w)}.
  }
  \item{method}{a character string specifying which method to use.  The
    default is \code{"M"}, for historical and back-compatibility
    reasons.  For the other methods, primarily see
    \code{\link{nlrob.algorithms}}. % nlrob-algos.Rd

    \describe{
      \item{"M"}{Computes an M-estimator, using \code{\link{nls}(*,
	  weights=*)} iteratively (hence, IRLS) with weights equal to
	\eqn{\psi(r_i) / r_i}, where \eqn{r_i} is the i-the residual
	from the previous fit.}
      \item{"MM"}{Computes an MM-estimator, starting from \code{init},
	either "S" or "lts".}% more: FIXME
      \item{"tau"}{Computes a Tau-estimator.}
      \item{"CM"}{Computes a \dQuote{Constrained M} (=: CM) estimator.}
      \item{"mtl"}{Compute as \dQuote{Maximum Trimmed Likelihood} (=: MTL)
	estimator.}
    }
    Note that all methods but \code{"M"} are \dQuote{random}, hence
    typically to be preceded by \code{\link{set.seed}()} in usage, see
    also \code{\link{nlrob.algorithms}}. % nlrob-algos.Rd
  }
  \item{psi}{a function (possibly by name) of the form \code{g(x, 'tuning
      constant(s)', deriv)} that for \code{deriv=0} returns
    \eqn{\psi(x)/x}{psi(x)/x} and for \code{deriv=1} returns
    \eqn{\psi'(x)}{psi'(x)}.  Note that tuning constants can \emph{not}
    be passed separately, but directly via the specification of \code{psi},
    typically via a simple \code{\link{.Mwgt.psi1}()} call as per
    default.

    Note that this has been a deliberately non-backcompatible change
    for robustbase version 0.90-0 (summer 2013 -- early 2014).
  }
  \item{scale}{when not \code{NULL} (default), a positive number
    specifying a scale kept \emph{fixed} during the iterations (and
    returned as \code{Scale} component).}
  \item{test.vec}{character string specifying the convergence
    criterion. The relative change is tested for residuals with a value
    of \code{"resid"} (the default), for coefficients with
    \code{"coef"}, and for weights with \code{"w"}.}
  \item{maxit}{maximum number of iterations in the robust loop.}
  \item{tol}{non-negative convergence tolerance for the robust fit.}
  \item{acc}{previous name for \code{tol}, now deprecated.}
  \item{algorithm}{character string specifying the algorithm to use for
    \code{\link{nls}}, see there, only when \code{method = "M"}.  The
    default algorithm is a Gauss-Newton algorithm.}
  \item{doCov}{a logical specifying if \code{nlrob()} should compute the
    asymptotic variance-covariance matrix (see \code{\link{vcov}})
    already.  This used to be hard-wired to \code{TRUE}; however, the
    default has been set to \code{FALSE}, as \code{\link{vcov}(obj)} and
    \code{\link{summary}(obj)} can easily compute it when needed.}
  \item{control}{an optional list of control settings.
    \describe{
      \item{for \code{method = "M"}:}{settings for \code{\link{nls}()}.
	See \code{\link{nls.control}} for the names of the settable control
	values and their effect.}
      \item{for all \code{method}s but \code{"M"}:}{a list, typically
	resulting from \code{\link{nlrob.control}(method, *)}.}
    }
  }
  \item{trace}{logical value indicating if a \dQuote{trace} of
    the \code{nls} iteration progress should be printed.  Default is
    \code{FALSE}. \cr
    If \code{TRUE}, in each robust iteration, the residual
    sum-of-squares and the parameter values are printed at the
    conclusion of each \code{nls} iteration.
    When the \code{"plinear"} algorithm is used, the conditional
    estimates of the linear parameters are printed after the nonlinear
    parameters.}

  \item{object}{an \R object of class \code{"nlrob"}, typically
    resulting from \code{nlrob(..)}.}
  \item{\dots}{for \code{nlrob}: only when \code{method} is \emph{not} \code{"M"},
    optional arguments for \code{\link{nlrob.control}};
    \cr
    for other functions: potentially optional arguments passed to the
    extractor methods.}
  \item{type}{a string specifying the \emph{type} of residuals desired.
    Currently, \code{"response"} and \code{"working"} are supported.
    %% FIXME: 1. document these (here)  2. write and support more types
  }
  \item{newdata}{a data frame (or list) with the same names as the
    original \code{data}, see e.g., \code{\link{predict.nls}}.}
}
\details{
  For \code{method = "M"}, iterated reweighted least squares
  (\dQuote{IRLS} or \dQuote{IWLS}) is used, calling \code{\link{nls}(*,
    weights= .)} where \code{weights} \eqn{w_i} are proportional to
  \eqn{\psi(r_i/ \hat{\sigma})}{psi(r_i/ sig.)}.

  All other methods minimize differently, and work \bold{without}
  \code{\link{nls}}.  See \link{nlrob.algorithms} % -> nlrob-algos.Rd
  for details.
}
\value{
  \code{nlrob()} returns an object of S3 class \code{"nlrob"},
  for \code{method = "M"} also inheriting from class \code{"nls"}, (see
  \code{\link{nls}}).

  It is a list with several components; they are not documented yet,
  as some of them will probably change.
  Instead, rather use \dQuote{accessor} methods, where possible:
  There are methods (at least) for the generic accessor functions
  \code{\link{summary}()}, \code{\link{coefficients}()} (aka \code{coef()})
  \code{fitted.values()}, \code{residuals()}, \code{\link{sigma}()} and
  \code{\link{vcov}()}, the latter for the variance-covariance matrix of
  the estimated parameters, as returned by \code{coef()}, i.e., not
  including the variance of the errors.
  For \code{nlrob()} results, \code{\link{estimethod}()} returns the
  \dQuote{estimation method}, which coincides with the \code{method}
  argument used.

  \code{residuals(.)}, by default \code{type = "response"}, returns
  the residuals \eqn{e_i}, defined above as
  \eqn{e_i = Y_i - f_(x_i, \hat\theta)}{e[i] = Y[i] - f(x[i], theta^)}.
  These differ from the standardized or weighted residuals which, e.g.,
  are assumed to be normally distributed, and a version of which is
  returned in \code{working.residuals} component.
  %% and another is working.residuals/Scale
}

\author{
  \describe{
    \item{\code{method = "M"}:}{
      Andreas Ruckstuhl (inspired by \code{\link[MASS]{rlm}}() and
      \code{\link{nls}}()), in July 1994 for S-plus.\cr
      Christian Sangiorgio did the update to \R and corrected some errors,
      from June 2002 to January 2005, and Andreas contributed slight changes
      and the first methods in August 2005.}
    \item{\code{method = "MM"}, etc:}{Originally all by Eduardo
      L. T. Conceicao, see \code{\link{nlrob.algorithms}}:} % nlrob-algos.Rd
  }

  Since then, the help page, testing, more cleanup, new methods: Martin
  Maechler.
}
\note{
  This function (with the only method \code{"M"}) used to be named
  \code{rnls} and has been in package \pkg{sfsmisc} in the past, but
  been dropped there.
}
\seealso{ \code{\link{nls}}, \code{\link[MASS]{rlm}}.
}
\examples{
DNase1 <- DNase[ DNase$Run == 1, ]

## note that selfstarting models don't work yet % <<< FIXME !!!

##--- without conditional linearity ---

## classical
fmNase1 <- nls( density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                data = DNase1,
                start = list( Asym = 3, xmid = 0, scal = 1 ),
                trace = TRUE )
summary( fmNase1 )

## robust
RmN1  <- nlrob( density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                data = DNase1, trace = TRUE,
                start = list( Asym = 3, xmid = 0, scal = 1 ))
summary( RmN1 )

##--- using conditional linearity ---

## classical
fm2DNase1 <- nls( density ~ 1/(1 + exp(( xmid - log(conc) )/scal ) ),
                  data = DNase1,
                  start = c( xmid = 0, scal = 1 ),
                  alg = "plinear", trace = TRUE )
summary( fm2DNase1 )

## robust
frm2DNase1 <- nlrob(density ~ 1/(1 + exp(( xmid - log(conc) )/scal ) ),
                  data = DNase1, start = c( xmid = 0, scal = 1 ),
                  alg = "plinear", trace = TRUE )
summary( frm2DNase1 )
## Confidence for linear parameter is quite smaller than "Asym" above
c1 <- coef(summary(RmN1))
c2 <- coef(summary(frm2DNase1))
rownames(c2)[rownames(c2) == ".lin"] <- "Asym"
stopifnot(all.equal(c1[,1:2], c2[rownames(c1), 1:2], tol = 0.09)) # 0.07315

### -- new examples -- "moderate outlier":
DN2 <- DNase1
DN2[10,"density"] <- 2*DN2[10,"density"]

fm3DN2 <- nls(density ~  Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
              data = DN2, trace = TRUE,
              start = list( Asym = 3, xmid = 0, scal = 1 ))

## robust
Rm3DN2 <- nlrob(density ~  Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                data = DN2, trace = TRUE,
                start = list( Asym = 3, xmid = 0, scal = 1 ))
Rm3DN2
summary(Rm3DN2) # -> robustness weight of obs. 10 ~= 0.037
confint(Rm3DN2, method = "Wald")

## utility function sfsmisc::lseq() :
lseq <- function (from, to, length)
  2^seq(log2(from), log2(to), length.out = length)
## predict() {and plot}:
h.x <- lseq(min(DN2$conc), max(DN2$conc), length = 100)
nDat <- data.frame(conc = h.x)

h.p  <- predict(fm3DN2, newdata = nDat)# classical
h.rp <- predict(Rm3DN2, newdata = nDat)# robust

plot(density ~ conc, data=DN2, log="x",
     main = format(formula(Rm3DN2)))
lines(h.x, h.p,  col="blue")
lines(h.x, h.rp, col="magenta")
legend("topleft", c("classical nls()", "robust nlrob()"),
       lwd = 1, col= c("blue", "magenta"), inset = 0.05)

## See  ?nlrob.algorithms for examples
\donttest{
DNase1 <- DNase[DNase$Run == 1,]
form <- density ~ Asym/(1 + exp(( xmid -log(conc) )/scal ))
gMM  <- nlrob(form, data = DNase1, method = "MM",
              lower = c(Asym = 0, xmid = 0, scal = 0),
              upper = 3, trace = TRUE)

## "CM" (and "mtl") additionally need bounds for "sigma" :
gCM  <- nlrob(form, data = DNase1, method = "CM",
              lower = c(Asym = 0, xmid = 0, scal = 0, sigma = 0),
              upper = c(3,3,3, sigma = 0.8))
summary(gCM)# (did fail)
}%not (always) tested
}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
