\name{psi.findc}
\title{Find Tuning Constant(s) for "lqq" and "ggw" Psi Functions}
\alias{.psi.ggw.findc}% --> ../R/lmrob.MM.R
\alias{.psi.lqq.findc}
\description{
  Find psi function tuning constant sets for \code{"LQQ"} and \code{"GGW"}
  psi (\eqn{\psi}) functions by specifying largest descent (minimal
  slope), efficiency and or breakdown point.
}
\usage{
.psi.ggw.findc(ms, b, eff = NA, bp = NA,
               subdivisions = 100L,
               rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
               tol = .Machine$double.eps^0.25, ms.tol = tol/64, maxiter = 1000)

.psi.lqq.findc(ms, b.c, eff = NA, bp = NA,
               interval = c(0.1, 4), subdivisions = 100L,
               rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
               tol = .Machine$double.eps^0.25, maxiter = 1000)
}
\arguments{
  \item{ms}{number, the minimal slope, typically negative.}
  \item{b, b.c}{number, specifying \eqn{b} or \eqn{b/c} for \code{"ggw"}
    or \code{"lqq"} respectively.}
  \item{eff}{a number (or \code{\link{NA}}), the desired
    \emph{efficiency}, in \eqn{[0,1]} of the estimator.  If \code{NA},
    \code{bp} must be specified as valid number.}
  \item{bp}{a number (or \code{\link{NA}}), the desired
    \emph{\bold{b}reakdown \bold{p}oint} of the estimator, in \eqn{[0,1]}.}

  \item{interval}{for finding \eqn{c} via \code{\link{uniroot}()}.}
  \item{subdivisions}{passed to \code{\link{integrate}()}.}
  \item{rel.tol, abs.tol}{relative and absolute tolerance for \code{\link{integrate}()}.}
  \item{tol}{relative tolerance for \code{\link{uniroot}()}.}
  \item{ms.tol}{relative tolerance for the internal
    \code{.psi.ggw.finda()}, eventually passed to \code{\link{optimize}}
    inside (internal) \code{.psi.ggw.mxs()}.}
  \item{maxiter}{maximal number of iterations for \code{\link{uniroot}()}.}
}
\value{
  constants for \code{"lqq"} psi functions: \eqn{(b, c, s) = (b/c * c, c, s = 1 - min_slope)}
}
\seealso{\code{\link{Mpsi}()} etc for the psi function definitions;
  \code{\link{.Mpsi.tuning.defaults}}, etc, for tuning constants' efaults for\code{\link{lmrob}()}.
}
\author{Manuel Koller (original) and Martin Maechler (arguments, export,
  docs).
}
\examples{
(c.ge95 <- .psi.ggw.findc(ms = -.5, b = 1.5, eff = 0.95))
(c.ge90 <- .psi.ggw.findc(ms = -.5, b = 1.5, eff = 0.90))
(c.gb50 <- .psi.ggw.findc(ms = -.5, b = 1.5,  bp = 0.50))
stopifnot(all.equal(c.ge95, c(0, 1.386362,   1.5, 1.0628199,  4.7773893), tol = 1e-5),
          all.equal(c.ge90, c(0, 1.0282811,  1.5, 0.87086259, 3.2075233), tol = 1e-5),
          all.equal(c.gb50, c(0, 0.20367394, 1.5, 0.29591308, 0.37033962),tol = 1e-5))

(cl.e.95 <- .psi.lqq.findc(ms = -.5, b.c = 1.5, eff = .95))
(cl.b.50 <- .psi.lqq.findc(ms = -.5, b.c = 1.5,  bp = .50))
stopifnot(all.equal(cl.e.95, c(1.4734061,  0.98227073, 1.5), tol = 1e-5),
          all.equal(cl.b.50, c(0.40154568, 0.26769712, 1.5), tol = 1e-5))
}
\keyword{utilities}
