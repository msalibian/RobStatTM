\name{Mpsi}
\title{Psi / Chi / Wgt / Rho Functions for *M-Estimation}
\alias{Mchi}
\alias{Mpsi}
\alias{Mwgt}
\alias{MrhoInf}
\alias{.Mchi}
\alias{.Mpsi}
\alias{.Mwgt}
\alias{.Mwgt.psi1}
\alias{.MrhoInf}
\alias{.psi2ipsi}
\description{
  Compute Psi / Chi / Wgt / Rho functions for M-estimation,
  i.e., including MM, etc.
  %% TODO: More, notably definitions ... but they are all nicely in the
  %% vignette.. How can we link from here to there ???

  \code{MrhoInf(x)} computes \eqn{\rho(\infty)}{rho(Inf)}, i.e., the
  normalizing or scaling constant for the transformation
  from \eqn{\rho(\cdot)}{rho(.)} to
  \eqn{\tilde\rho(\cdot)}{rho~(.)}, where the latter, aka as
  \eqn{\chi()}{chi()} fulfills \eqn{\tilde\rho(\infty) = 1}{rho~(Inf) = 1}
  which makes only sense for \dQuote{redescending} psi functions, i.e.,
  not for \code{"huber"}.

  \code{Mwgt(x, *)} computes \eqn{\psi(x)/x}  (fast and numerically accurately).
}
\usage{
Mpsi(x, cc, psi, deriv = 0)
Mchi(x, cc, psi, deriv = 0)
Mwgt(x, cc, psi)
MrhoInf(cc, psi)

.Mwgt.psi1(psi, cc = .Mpsi.tuning.default(psi))
}
\arguments{
  \item{x}{numeric (\dQuote{abscissa} values) vector, possibly with
    \code{\link{attributes}} such as \code{\link{dim}} or
    \code{\link{names}}, etc.  These are preserved for the
    \code{M*()} functions (but not the \code{.M()} ones).}
  \item{cc}{numeric tuning constant, for some \code{psi} of length
    \eqn{> 1}.}
  \item{psi}{a string specifying the psi / chi / rho / wgt function;
    either \code{"huber"}, or one of the same possible specifiers as for
    \code{psi} in \code{\link{lmrob.control}}, i.e. currently,
    \code{"bisquare"}, \code{"lqq"}, \code{"welsh"}, \code{"optimal"},
    \code{"hampel"}, or \code{"ggw"}.}
  \item{deriv}{an integer, specifying the \emph{order} of derivative to
    consider; particularly, \code{Mpsi(x, *, deriv = -1)} is the
    principal function of \eqn{\psi()}{psi()}, typically denoted
    \eqn{\rho()}{rho()} in the literature.}
  %% FIXME: mention that  deriv = 2 is *partially* implemented
}
\details{
  Theoretically, \code{Mchi()} would not be needed explicitly as it can be computed
  from \code{Mpsi()} and \code{MrhoInf()}, namely, by
  \preformatted{Mchi(x, *, deriv = d)  ==  Mpsi(x, *, deriv = d-1) / MrhoInf(*)}
  for \eqn{d = 0, 1, 2}  (and \sQuote{*} containing \code{par, psi}, and
  equality is in the sense of \code{\link{all.equal}(x,y, tol)} with a
  small \code{tol}.

  Similarly, \code{Mwgt} would not be needed strictly, as it could be
  defined via \code{Mpsi}), but the explicit definition takes care of
  0/0 and typically is of a more simple form.

  For experts, there are slightly even faster versions,
  \code{.Mpsi()}, \code{.Mwgt()}, etc.

  \code{.Mwgt.psi1()} mainly a utility for \code{\link{nlrob}()},
  returns a \emph{\code{\link{function}}} with similar semantics as
  \code{\link[MASS]{psi.hampel}}, \code{\link[MASS]{psi.huber}}, or
  \code{\link[MASS]{psi.bisquare}} from package \pkg{MASS}.  Namely,
  a function with arguments \code{(x, deriv=0)}, which for
  \code{deriv=0} computes \code{Mwgt(x, cc, psi)} and otherwise computes
  \code{Mpsi(x, cc, psi, deriv=deriv)}.

  \code{.Mpsi()}, \code{.Mchi()}, \code{.Mwgt()}, and \code{.MrhoInf()} are
  low-level versions of
  \code{Mpsi()},  \code{Mchi()},  \code{Mwgt()}, and  \code{MrhoInf()}, respectively,
  and \code{.psi2ipsi()} provides the psi-function integer codes needed
  for \code{ipsi} argument of the \code{.M*()} functions.
}
\value{
  a numeric vector of the same length as \code{x}, with corresponding
  function (or derivative) values.
}
\references{
  See the vignette about %% \link{ .. } ????
  \dQuote{\eqn{\psi}{psi}-Functions Available in Robustbase}.
  %% ../vignettes/psi_functions.Rnw
}
\author{
  Manuel Koller, notably for the original C implementation;
  tweaks and speedup via \code{\link{.Call}} and \code{.M*()} etc by
  Martin Maechler.
}
\seealso{
 \code{\link{psiFunc}} and the \code{\linkS4class{psi_func}} class, both
 of which provide considerably more on the \R side, but are less
 optimized for speed.

 \code{\link{.Mpsi.tuning.defaults}}, etc, for tuning constants'
 defaults for\code{lmrob()}, and \code{\link{.psi.ggw.findc}()}
 utilities to construct such constants' vectors.
}
\examples{
x <- seq(-5,7, by=1/8)
matplot(x, cbind(Mpsi(x, 4, "biweight"),
                 Mchi(x, 4, "biweight"),
                 Mwgt(x, 4, "biweight")), type = "l")
abline(h=0, v=0, lty=2, col=adjustcolor("gray", 0.6))

hampelPsi
(ccHa <- hampelPsi @ xtras $ tuningP $ k)
psHa <- hampelPsi@psi(x)
% FIXME: interesting as long as hampelPsi does not use Mpsi(... "hampel") !
## using Mpsi():
Mp.Ha <- Mpsi(x, cc = ccHa, psi = "hampel")
stopifnot(all.equal(Mp.Ha, psHa, tolerance = 1e-15))

psi.huber <- .Mwgt.psi1("huber")
if(getRversion() >= "3.0.0")
stopifnot(identical(psi.huber, .Mwgt.psi1("huber", 1.345),
                    ignore.env=TRUE))
curve(psi.huber(x), -3, 5, col=2, ylim = 0:1)
curve(psi.huber(x, deriv=1), add=TRUE, col=3)

## and show that this is indeed the same as  MASS::psi.huber() :
x <- runif(256, -2,3)
stopifnot(all.equal(psi.huber(x), MASS::psi.huber(x)),
          all.equal(                 psi.huber(x, deriv=1),
                    as.numeric(MASS::psi.huber(x, deriv=1))))

## and how to get  MASS::psi.hampel():
psi.hampel <- .Mwgt.psi1("Hampel", c(2,4,8))
x <- runif(256, -4, 10)
stopifnot(all.equal(psi.hampel(x), MASS::psi.hampel(x)),
          all.equal(                 psi.hampel(x, deriv=1),
                    as.numeric(MASS::psi.hampel(x, deriv=1))))

## M*() preserving attributes :
x <- matrix(x, 32, 8, dimnames=list(paste0("r",1:32), col=letters[1:8]))
comment(x) <- "a vector which is a matrix"
px <- Mpsi(x, cc = ccHa, psi = "hampel")
stopifnot(identical(attributes(x), attributes(px)))

## The "optimal" psi exists in two versions "in the litterature": ---
## Maronna et al. 2006, 5.9.1, p.144f:
psi.M2006 <- function(x, c = 0.013)
  sign(x) * pmax(0, abs(x) - c/dnorm(abs(x)))
## and the other is the one in robustbase from 'robust': via Mpsi(.., "optimal")
## Here are both for 95\% efficiency:
(c106 <- .Mpsi.tuning.default("optimal"))
c1 <- curve(Mpsi(x, cc = c106, psi="optimal"), -5, 7, n=1001)
c2 <- curve(psi.M2006(x), add=TRUE, n=1001, col=adjustcolor(2,0.4), lwd=2)
abline(0,1, v=0, h=0, lty=3)
## the two psi's are similar, but really quite different

## a zoom into Maronna et al's:
c3 <- curve(psi.M2006(x), -.5, 1, n=1001); abline(h=0,v=0, lty=3);abline(0,1, lty=2)
}
\keyword{robust}
