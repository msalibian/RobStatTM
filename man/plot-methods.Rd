\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,psi_func-method}
\title{Plot an Object of the "Psi Function" Class}
\description{
  The \code{\link{plot}} method objects of class
  \code{\linkS4class{psi_func}} simply visualizes the
  \eqn{\rho()}, \eqn{\psi()}, and weight functions and their
  derivatives.
}
\usage{
\S4method{plot}{psi_func}(x, y,
     which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
     main = "full", % shortMain = FALSE,
     col = c("black", "red3", "blue3", "dark green", "light green"),
     leg.loc = "right", ...)
}
\arguments{
  \item{x}{object of class \code{\linkS4class{psi_func}} to be plotted}
  \item{y}{(optional) vector of abscissa values (to plot object at).}
  \item{which}{\code{\link{character}} vector of slots to be included in
    plot; by default, all of the slots are included}
  \item{main}{string or logical indicating the kind of plot title;
  either \code{"full"}, \code{"short"} or \code{FALSE} which chooses a
  full, a short or no main title at all.}
  %% \item{shortMain}{use short or long plot title}%-------- deprecated!
  \item{col}{colors to be used for the different slots}
  \item{leg.loc}{legend placement, see also \code{x} argument
    of \code{\link{legend}}}
  \item{...}{passed to \code{\link{matplot}}}
}
\note{
  An earlier version had argument \code{shortMain} which is deprecated
  now. Use \code{main = "short"} instead of \code{shortMain = TRUE}.

  If you want to specify your own title, use \code{main=FALSE}, and a
  subsequent \code{\link{title}(...)} call.
}
\seealso{
  \code{\link{psiFunc}()} and the \code{\link{class}}
  \code{\linkS4class{psi_func}}.
}
\examples{
plot(huberPsi)
plot(huberPsi, which=c("psi", "Dpsi", "wgt"),
     main="short", leg = "topleft")

plot(hampelPsi)
## Plotting aspect ratio = 1:1 :
plot(hampelPsi, asp=1, main="short",
     which = c("psi", "Dpsi", "wgt", "Dwgt"))
}
\keyword{methods}
