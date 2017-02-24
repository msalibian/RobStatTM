\name{r6pack}
\alias{r6pack}
\title{Robust Distance based observation orderings based on robust "Six pack"}
\description{
  Compute six initial robust estimators of multivariate location and
  \dQuote{scatter} (scale); then, for each, compute the distances
  \eqn{d_{ij}}{d_ij} and take the \code{h} (\eqn{h > n/2}) observations
  with smallest distances.  Then compute the statistical distances based
  on these h observations.

  Return the indices of the observations sorted in increasing order.
}
\usage{
r6pack(x, h, full.h, scaled = TRUE, scalefn = rrcov.control()$scalefn)
}
\arguments{
  \item{x}{n x p data matrix}
  \item{h}{integer, typically around (and slightly larger than) \eqn{n/2}.}
  \item{full.h}{logical specifying if the full (length n) observation
    ordering should be returned; otherwise only the first \code{h} are.
    For \code{.detmcd()}, \code{full.h=FALSE} is typical.}
  \item{scaled}{logical indicating if the data \code{x} is
    already scaled; if false, we apply \code{x <- doScale(x, median,
      scalefn)}.}
  \item{scalefn}{a \code{\link{function}(u)} to compute a robust
    univariate scale of u.}
}% args
\details{%% --> ../R/detmcd.R
  The six initial estimators are
  \enumerate{
    \item{Hyperbolic tangent of standardized data}
    \item{Spearmann correlation matrix}
    \item{Tukey normal scores}
    \item{Spatial sign covariance matrix}
    \item{BACON}
    \item{Raw OGK estimate for scatter}
  }
}
\references{
  Hubert, M., Rousseeuw, P. J. and Verdonck, T. (2012)
  A deterministic algorithm for robust location and scatter.
  Journal of Computational and Graphical Statistics \bold{21}, 618--637.
}
\value{
  a \eqn{h' \times 6}{h' x 6} \code{\link{matrix}} of observation
  indices, i.e., with values from \eqn{1,\dots,n}{1..n}.  If
  \code{full.h} is true, \eqn{h' = n}, otherwise \eqn{h' = h}.
}
\author{Valentin Todorov, based on the original Matlab code by
  Tim Verdonck and Mia Hubert.  Martin Maechler for tweaks
  (performance etc), and \code{full.h}.
}
\seealso{
  \code{\link{covMcd}(*, nsamp = "deterministic")};
  \code{\link[rrcov]{CovSest}(*, nsamp = "sdet")} from package \pkg{rrcov}.
}
\examples{
data(pulpfiber)
dim(m.pulp <- data.matrix(pulpfiber)) #  62 x 8
dim(fr6  <- r6pack(m.pulp, h = 40, full.h= FALSE)) #  h x 6  = 40 x 6
dim(fr6F <- r6pack(m.pulp, h = 40, full.h= TRUE )) #  n x 6  = 62 x 6
stopifnot(identical(fr6, fr6F[1:40,]))
\dontshow{
stopifnot(apply(fr6[1:10,], 2L,
   function(col) c(1,3,6,35,36,38) \%in\% col))
}
}
\keyword{robust}
\keyword{multivariate}

