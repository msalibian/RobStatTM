\name{wagnerGrowth}
\encoding{utf8}
\Rdversion{1.1}
\alias{wagnerGrowth}
\docType{data}
\title{
  Wagner's Hannover Employment Growth Data
}
\description{
  Wagner (1994) investigates the rate of employment growth (\code{y}) as
  function of percentage of people engaged in \bold{p}roducation
  \bold{a}ctivities (\code{PA}) and \bold{h}igher \bold{s}ervices
  (\code{HS}) and of the \bold{g}rowth of these percentages (\code{GPA},
  \code{GHS}) during three time periods in 21 geographical regions of
  the greater Hannover area.
}
\usage{data(wagnerGrowth)}
\format{
  A data frame with \eqn{21 \times 3 = 63}{21 * 3 = 63} observations
  (one per \code{Region x Period}) on the following 7 variables.
  \describe{
    \item{\code{Region}}{a \code{\link{factor}} with 21 levels, denoting
      the corresponding region in Hannover (conceptually a \dQuote{block
	factor}).}
    \item{\code{PA}}{numeric: percent of people involved in production activities.}
    \item{\code{GPA}}{\bold{g}rowth of \code{PA}.}
    \item{\code{HS}}{a numeric vector}
    \item{\code{GHS}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{Period}}{a \code{\link{factor}} with levels \code{1:3},
      denoting the time period, 1 = 1979-1982, 2 = 1983-1988, 3 = 1989-1992.}
  }
}
% \details{

% }
\source{
  Hubert, M. and Rousseeuw, P. J.  (1997).
  Robust regression with both continuous and binary regressors,
  \emph{Journal of Statistical Planning and Inference} \bold{57}, 153--163.
}
\references{
  Wagner J. (1994).
  Regionale Beschäftigungsdynamik und höherwertige
  Produktionsdienste: Ergebnisse für den Grossraum Hannover (1979-1992).
  \emph{Raumforschung und Raumordnung} \bold{52}, 146--150.
}
\examples{
data(wagnerGrowth)
## maybe
str(wagnerGrowth)


require(lattice)
(xyplot(y ~ Period | Region, data = wagnerGrowth,
         main = "wagnerGrowth: 21 regions @ Hannover"))

(dotplot(y ~ reorder(Region,y,median), data = wagnerGrowth,
         main = "wagnerGrowth",
         xlab = "Region [ordered by  median(y | Region) ]"))
}
\keyword{datasets}
