\name{psiFunc}
\alias{psiFunc}
\alias{huberPsi}
\alias{hampelPsi}
\title{Constructor for Objects "Psi Function" Class}
%% Rd Problem:  'psi_func' cannot be used in title!
\description{
  \code{psiFunc(..)} is a convenience interface to
  \code{new("psi_func",..)}, i.e. for constructing objects of class
  \code{"psi_func"}.
}
\usage{
psiFunc(rho, psi, wgt, Dpsi,Dwgt, Erho = NULL, Epsi2 = NULL, EDpsi = NULL, name, ...)

huberPsi
hampelPsi
}
\arguments{
  \item{rho, psi, wgt, Dpsi, Dwgt}{each a \code{\link{function}} of \code{x} and
    tuning parameters typically. Specification of Dwgt is optional.}
  \item{Erho, Epsi2, EDpsi}{see \code{\linkS4class{psi_func}}, and note
    that these may change in the future.}
  \item{name}{Name of \eqn{\psi}{psi}-function used for printing.}
  \item{\dots}{potential further arguments for specifying tuning
    parameter names and defaults.} %% FIXME; give more details
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
%%\references{ ~put references to the literature/web site here ~ }
\author{Martin Maechler}
\seealso{The description of class \code{\linkS4class{psi_func}}.
}
\examples{
## classical {trivial, not interesting}:
F1 <- function(x, .) rep.int(1, length(x))
FF <- function(.) rep.int(1, length(.))
cPsi <- psiFunc(rho = function(x,.) x^2 / 2, psi = function(x, .) x,
                wgt = F1, Dpsi = F1,
                Erho = function(.) rep.int(1/2, length(.)),
                Epsi2 = FF, EDpsi = FF, name = "classic", . = Inf)
show(cPsi)

## Hampel's psi and rho:
H.38 <- chgDefaults(hampelPsi, k = c(1.5, 3.5, 8))
k. <- H.38@xtras$tuningP$k ; k.. <- as.vector(outer(c(-1,1), k.))
c.t <- adjustcolor("skyblue3", .8)
.ax.k <- function(side) { abline(h=0, v=0, lty=2)
  axis(side, at = k.., labels=formatC(k..), pos=0, col=c.t, col.axis=c.t) }
op <- par(mfrow=c(2,1), mgp = c(1.5, .6, 0), mar = .6+c(2,2,1,.5))
curve(H.38@psi(x), -10, 10, col=2, lwd=2, n=512)
lines(k.., H.38@psi(k..), type = "h", lty=3, col=c.t); .ax.k(1)
curve(H.38@rho(x), -10, 10, col=2, lwd=2, n=512); abline(h=0, v=0, lty=2)
lines(k.., H.38@rho(k..), type = "h", lty=3, col=c.t); .ax.k(1)
title(expression("Hampel's " ~~~ psi(x) ~~ "and" ~~ rho(x) ~~~ " functions"))
par(op)

## Not the same, but similar, directly using the plot() method:
plot(H.38)
}
\keyword{classes}
\keyword{robust}
