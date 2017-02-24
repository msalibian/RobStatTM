\name{plot.lmrob}
\alias{plot.lmrob}
\title{Plot Method for "lmrob" Objects}
\description{
  Diagnostic plots for elements of class lmrob
}
\usage{
\method{plot}{lmrob}(x, which = 1:5,
     caption = c("Standardized residuals vs. Robust Distances",
   "Normal Q-Q vs. Residuals", "Response vs. Fitted Values",
   "Residuals vs. Fitted Values" , "Sqrt of abs(Residuals) vs. Fitted Values"),
   panel = if(add.smooth) panel.smooth else points,
   sub.caption = deparse(x$call), main = "",
   compute.MD = TRUE,
   ask = prod(par("mfcol")) < length(which) && dev.interactive(),
   id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
   label.pos = c(4,2), qqline = TRUE, add.smooth = getOption("add.smooth"),
   \dots, p=0.025)
}
\arguments{
  \item{x}{ an object as created by \code{lmrob} }
  \item{which}{ integer number between 1 and 5 to specify which
    plot is desired }
  \item{caption}{Caption for the different plots}
  \item{panel}{panel function.  The useful alternative to
    \code{\link{points}}, \code{\link{panel.smooth}} can be chosen
    by \code{add.smooth = TRUE}.}
  \item{main}{main title}
  \item{sub.caption}{sub titles}
  \item{compute.MD}{logical indicating if the robust Mahalanobis
    distances should be recomputed, using \code{\link{covMcd}()} when
    needed, i.e., if \code{which} contains \code{1}.}
  \item{ask}{waits for user input before displaying each plot }
  \item{id.n}{number of points to be labelled in each plot, starting
    with the most extreme.}
  \item{labels.id}{vector of labels, from which the labels for extreme
    points will be chosen.  \code{NULL} uses observation numbers.}
  \item{cex.id}{magnification of point labels.}
  \item{label.pos}{positioning of labels, for the left half and right
    half of the graph respectively.}% for plots ...
  \item{qqline}{logical indicating if a \code{\link{qqline}()} should be
    added to the normal Q-Q plot.}
  \item{add.smooth}{logical indicating if a smoother should be added to
    most plots; see also \code{panel} above.}
  \item{\dots}{optional arguments for \code{\link{par}},
    \code{\link{title}}, etc.}
  \item{p}{threshold for distance-distance plot}
}
\references{
  Robust diagnostic plots as in Rousseeuw and van Zomeren (1990), see
  \sQuote{References} in \code{\link{ltsPlot}}.
}
\details{
  if \code{compute.MD = TRUE} and the robust Mahalanobis distances need
  to be computed, they are stored (\dQuote{cached}) with the object
  \code{x} when this function has been called from top-level.
}
\seealso{ \code{\link{lmrob}}, also for examples, \code{\link{plot.lm}}.
}
\examples{
data(starsCYG)
## Plot simple data and fitted lines
plot(starsCYG)
 lmST <-    lm(log.light ~ log.Te, data = starsCYG)
RlmST <- lmrob(log.light ~ log.Te, data = starsCYG)
RlmST
abline(lmST, col = "red")
abline(RlmST, col = "blue")

op <- par(mfrow = c(2,2), mgp = c(1.5, 0.6, 0), mar= .1+c(3,3,3,1))
plot(RlmST, which = c(1:2, 4:5))
par(op)
}
\keyword{robust}
\keyword{regression}

