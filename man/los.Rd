\name{los}
\alias{los}
\docType{data}
\title{ Length of Stay Data }
\description{
  Length of stay for 201 patients that stayed at the University Hospital
  of Lausanne during the year 2000.
}
\usage{data(los)}
\format{
  Vector of integer values giving the length of stay (days):

  int [1:201] 16 13 17 4 15 24 59 18 33 8 ...
}
\details{
  These data may be used to estimate and predict the total resource
  consumption of this group of patients.

  Cf. Ruffieux, Paccaud and Marazzi (2000).
}
\source{
  The data were kindly provided by A. Marazzi.

  Cf. Hubert, M. and Vandervieren, E. (2006), p. 13--15.
}
\references{
  Ruffieux, C., Paccaud, F. and A. Marazzi (2000)
  Comparing rules for truncating hospital length of stay;
  \emph{Casemix Quarterly} \bold{2}, n. 1.

  See also those for \code{\link{adjbox}}.
}
\examples{
 summary(los) # quite skewed, with median(.) = 8
 plot(table(los))
 boxplot(los, horizontal=TRUE, add=TRUE, col = "red", axes=FALSE)
 ##-> "outliers" instead of "just skewed"

 hist(log(los))
 boxplot(log(los), add=TRUE, col=2, border=2, horizontal = TRUE, at = -1)

 ## Hubert and Vandervieren (2006), p. 15, Fig. 11.
 adjbox(los, col = "gray", staplecol="red", outcol = "red",
        main = "(Skewness-)Adjusted and original boxplot for 'los' data")
 boxplot(los, add = TRUE, staplewex= 0.2, outcex= 0.5, outpch= 4,
         staplecol = "blue", outcol = "blue", staplelwd=2)
 legend("topright", c("adjbox(los)", "boxplot(los)"),
        col=c("red","blue"), lwd = 1:2, bty="n")
}
\keyword{datasets}
