\name{smoothWgt}
\alias{smoothWgt}
\title{Smooth Weighting Function - Generalized Biweight}% The Biweight on a Stick
\description{
  \dQuote{The Biweight on a Stick} ---
  Compute a smooth (when \eqn{h > 0}) weight function typically for
  computing weights from large (robust) \dQuote{distances} using a
  piecewise polynomial function which in fact is a
  2-parameter generalization of Tukey's 1-parameter \dQuote{biweight}.
}
\usage{
smoothWgt(x, c, h)
}
\arguments{
  \item{x}{numeric vector of abscissa values}
  \item{c}{\dQuote{cutoff}, a typically positive number.}
  \item{h}{\dQuote{bandwidth}, a positive number.}
}
\details{
  Let \eqn{w(x;c,h) := }\code{smoothWgt(x, c, h)}.  Then,
  \deqn{% "FIXME": rather use amsmath package \cases{.}
    w(x; c,h) := 0 \ \ \ \ \ \mathrm{if}\ |x| \ge c + h/2,}{%
    w(x; c,h) := 0   if |x| >= c + h/2,}
  \deqn{
    w(x; c,h) := 1 \ \ \ \ \ \mathrm{if}\ |x| \le c - h/2,}{%
    w(x; c,h) := 1   if |x| <= c - h/2,}
  \deqn{
    w(x; c,h) := \bigl((1 - |x| - (c-h/2))^2\bigr)^2 \ \mathrm{if}\ c-h/2 < |x| < c+h/2,}{%
    w(x; c,h) := (1 - (|x| - (c-h/2))^2)^2   if c-h/2 < |x| < c+h/2.}

  \code{smoothWgt()} is \emph{scale invariant} in the sense that
  \deqn{w(\sigma x; \sigma c, \sigma h) = w(x; c, h),}{%
    w(S x; S c, S h) = w(x; c, h),} when \eqn{\sigma > 0}{S > 0}.
}
\value{
  a numeric vector of the same length as \code{x} with weights between
  zero and one.  Currently all \code{\link{attributes}} including
  \code{\link{dim}} and \code{\link{names}} are dropped.
}
%% \references{ TODO: Write a small vignette !
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{Mwgt}(.., psi = "bisquare")} of which \code{smoothWgt()}
  is a generalization, and
  \code{\link{Mwgt}(.., psi = "optimal")} which looks similar for larger
  \code{c} with its constant one part around zero, but also has only
  one parameter.
}
\examples{
## a somewhat typical picture:
curve(smoothWgt(x, c=3, h=1), -5,7, n = 1000)

csW <- curve(smoothWgt(x, c=1/2, h=1), -2,2) # cutoff 1/2, bandwidth 1
## Show that the above is the same as
## Tukey's "biweight" or "bi-square" weight function:
bw <- function(x) pmax(0, (1 - x^2))^2
cbw <- curve(bw,                     col=adjustcolor(2, 1/2), lwd=2, add=TRUE)
cMw <- curve(Mwgt(x,c=1,"biweight"), col=adjustcolor(3, 1/2), lwd=2, add=TRUE)
stopifnot(## proving they are all the same:
   all.equal(csW, cbw, tol=1e-15),
   all.equal(csW, cMw, tol=1e-15))
}
\keyword{arith}
\keyword{robust}
