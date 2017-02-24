\name{cloud}
\alias{cloud}
\docType{data}
\title{Cloud point of a Liquid}
\description{
  This data set contains the measurements concerning the
  cloud point of a Liquid, from Draper and Smith (1969).
  The cloud point is a measure of the degree of crystallization in a
  stock.
}
\usage{data(cloud)}
\format{
  A data frame with 19 observations on the following 2 variables.
  \describe{
    \item{\code{Percentage}}{Percentage of I-8}
    \item{\code{CloudPoint}}{Cloud point}
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, p.96, table 10.
}
\examples{
data(cloud)
summary(lm.cloud <- lm(CloudPoint ~., data=cloud))
%% summary(lts.cloud <- ltsReg(CloudPoint ~., data=cloud))
}
\keyword{datasets}
