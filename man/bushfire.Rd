\name{bushfire}
\alias{bushfire}
\docType{data}
\title{ Campbell Bushfire Data }
\description{
  This data set was used by Campbell (1984) to locate bushfire scars.
  The dataset contains satelite measurements on five frequency bands,
  corresponding to each of 38 pixels.
}
\usage{data(bushfire)}
\format{
  A data frame with 38 observations on 5 variables.
%  \describe{}
}
%\Note{}
\source{
  Maronna, R.A. and Yohai, V.J. (1995)
  The Behavoiur of the Stahel-Donoho Robust Multivariate Estimator.
  \emph{Journal of the American Statistical Association} \bold{90}, 330--341.
}
%\seealso{}
\examples{
data(bushfire)
plot(bushfire)
covMcd(bushfire)
}
\keyword{datasets}
