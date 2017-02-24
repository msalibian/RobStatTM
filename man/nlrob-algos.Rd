\name{nlrob-algorithms}
\alias{nlrob.algorithms}
\alias{nlrob.MM}
\alias{nlrob.tau}
\alias{nlrob.CM}
\alias{nlrob.mtl}
\title{MM-, Tau-, CM-, and MTL- Estimators for Nonlinear Robust Regression}
\description{
  \describe{
    \item{"MM":}{Compute an MM-estimator for nonlinear robust (constrained)
      regression.}
    \item{"tau":}{Compute a Tau-estimator for nonlinear robust (constrained)
      regression.}
    \item{"CM":}{Compute a \dQuote{Constrained M} (=: CM) estimator for
      nonlinear robust (constrained) regression.}
    \item{"MTL":}{Compute a \dQuote{Maximum Trimmed Likelihood} (=: MTL)
      estimator for nonlinear robust (constrained) regression.}
  }
}
\usage{
## You can *not* call the  nlrob(*, method = <M>) like this ==> see  help(nlrob)
## ------- ===== ------------------------------------------

nlrob.MM(formula, data, pnames, lower, upper,
	 tol = 1e-06,
	 psi = c("bisquare", "lqq", "optimal", "hampel"),
         init = c("S", "lts"),
	 ctrl = nlrob.control("MM", psi = psi, init = init, fnscale = NULL,
		       tuning.chi.scale = .psi.conv.cc(psi, .Mchi.tuning.defaults[[psi]]),
		       tuning.psi.M     = .psi.conv.cc(psi, .Mpsi.tuning.defaults[[psi]]),
		       optim.control = list(), optArgs = list(...)),
	 ...)

nlrob.tau(formula, data, pnames, lower, upper,
	  tol = 1e-06, psi = c("bisquare", "optimal"),
	  ctrl = nlrob.control("tau", psi = psi, fnscale = NULL,
			tuning.chi.scale = NULL, tuning.chi.tau = NULL,
			optArgs = list(...)),
	  ...)

nlrob.CM(formula, data, pnames, lower, upper,
	 tol = 1e-06,
	 psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
	 ctrl = nlrob.control("CM", psi = psi, fnscale = NULL,
                        tuning.chi = NULL, optArgs = list(...)),
	 ...)

nlrob.mtl(formula, data, pnames, lower, upper,
	  tol = 1e-06,
	  ctrl = nlrob.control("mtl", cutoff = 2.5, optArgs = list(...)),
	  ...)
}
\arguments{
  \item{formula}{nonlinear regression \code{\link{formula}}, using both
    variable names from \code{data} and parameter names from \code{pnames}.}
  \item{data}{data to be used, a \code{\link{data.frame}}}
  \item{pnames}{a \code{\link{character}} vector of parameter names
    (used in \code{formula} above).  \bold{Deprecated}}
  \item{lower, upper}{bounds aka \dQuote{box constraints} for all the
    parameters, in the case "CM" and "mtl" these must include the error
    standard deviation as \code{"sigma"}, see \code{\link{nlrob}()}
    about its \code{\link{names}}, etc.}
  \item{tol}{numerical convergence tolerance.}
  \item{psi, init}{see \code{\link{nlrob.control}}.}
  \item{ctrl}{a \code{\link{list}}, typically the result of a call to
    \code{\link{nlrob.control}}.}

  \item{tuning.psi.M}{..}% FIXME
  \item{optim.control}{..}% FIXME
  \item{optArgs}{a \code{\link{list}} of optional arguments for
    optimization, e.g., \code{trace = TRUE}, passed to to the optimizer,
    which currently must be \code{\link{JDEoptim}(.)}.}
  \item{...}{alternative way to pass the \code{optArgs} above.}
}
\value{
  an \R object of \code{\link{class}} \code{"nlrob.<meth>"}, basically a
  list with components
  %% FIXME
}
\details{
  Copyright 2013, Eduardo L. T. Conceicao.  Available under
  the GPL (>= 2)
}
\author{
  Eduardo L. T. Conceicao; compatibility (to \code{\link{nlrob}})
  tweaks and generalizations, inference, by Martin Maechler.
}
\source{
  For \code{"MTL"}:
  Maronna, Ricardo A., Martin, R. Douglas, and Yohai, Victor J. (2006).
  \emph{Robust Statistics: Theory and Methods} Wiley, Chichester, p. 133.
}
\references{
  \describe{
    \item{"MM":}{
      Yohai, V.J. (1987)
      High breakdown-point and high efficiency robust estimates for
      regression.
      \emph{The Annals of Statistics} \bold{15}, 642--656.
    }
    \item{"tau":}{
      Yohai, V.J., and Zamar, R.H. (1988).
      High breakdown-point estimates of regression by means of the
      minimization of an efficient scale.
      \emph{Journal of the American Statistical Association} \bold{83},
      406--413.
    }
    \item{"CM":}{
      Mendes, B.V.M., and Tyler, D.E. (1996)
      Constrained M-estimation for regression.
      In: \emph{Robust Statistics, Data Analysis and Computer Intensive
	Methods}, Lecture Notes in Statistics 109, Springer, New York, 299--320.
      %% not yet -- e.g. tuning constants for Welsh:
      %% Edlund, O. and Ekblom, H. (2005)
      %% Computing the constrained M-estimates for regression.
      %% Computational Statistics Data Analysis \bold{49}(1): 19--32.
    }
    \item{"MTL":}{
      Hadi, Ali S., and Luceno, Alberto (1997).
      Maximum trimmed likelihood estimators: a unified approach,
      examples, and algorithms.
      Computational Statistics & Data Analysis \bold{25}, 251--272.

      Gervini, Daniel, and Yohai, Victor J. (2002).
      A class of robust and fully efficient regression estimators.
      The Annals of Statistics \bold{30}, 583--616.
    }

  }%describe
}
\examples{%% for more, --> ../tests/nlregrob-tst.R
DNase1 <- DNase[DNase$Run == 1,]
form <- density ~ Asym/(1 + exp(( xmid -log(conc) )/scal ))
pnms <- c("Asym", "xmid", "scal")
set.seed(47) # as these by default use randomized optimization:

fMM <- robustbase:::nlrob.MM(form, data = DNase1,
           lower = setNames(c(0,0,0), pnms), upper = 3,
           ## call to nlrob.control to pass 'optim.control':
           ctrl = nlrob.control("MM", optim.control = list(trace = 1),
                                optArgs = list(trace = TRUE)))

## The same via nlrob() {recommended; same random seed to necessarily give the same}:
set.seed(47)
gMM  <- nlrob(form, data = DNase1, method = "MM",
              lower = setNames(c(0,0,0), pnms), upper = 3, trace = TRUE)
gMM
summary(gMM)
## and they are the same {apart from 'call' and 'ctrl' and new stuff in gMM}:
ni <- names(fMM); ni <- ni[is.na(match(ni, c("call","ctrl")))]
stopifnot(all.equal(fMM[ni], gMM[ni]))
\dontshow{
if(doExtras <- robustbase:::doExtras()) {
 gtau <- nlrob(form, data = DNase1, method = "tau",
	       lower = setNames(c(0,0,0), pnms), upper = 3, trace = TRUE)

 ## these two have "sigma" also as parameter :
 psNms <- c(pnms, "sigma")
 gCM  <- nlrob(form, data = DNase1, method = "CM",
	       lower = setNames(c(0,0,0,0), psNms), upper = 3, trace = TRUE)
 gmtl <- nlrob(form, data = DNase1, method = "mtl",
	       lower = setNames(c(0,0,0,0), psNms), upper = 3, trace = TRUE)
 stopifnot(identical(sapply(list(gMM, gCM, gmtl), estimethod),
                     c("MM", "CM", "mtl")))
}% doExtras
}% dontshow
}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
