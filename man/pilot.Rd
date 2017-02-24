\name{pilot}
\alias{pilot}
\docType{data}
\title{Pilot-Plant Data}
\description{
  Pilot-Plant data from Daniel and Wood (1971). The response variable
  corresponds to the acid content determined by titration and the
  explanatory variable is the organic acid content determined by
  extraction and weighing. This data set was analyzed also by Yale and
  Forsythe (1976).
}
\usage{data(pilot)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{X}}{Organic acid content - extraction}
    \item{\code{Y}}{Acid content - titration }
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, page 21, table 1.
}
\examples{
data(pilot)
summary(lm.pilot <- lm(Y ~.,data=pilot))
%%summary(lts.pilot <- ltsReg(Y~.,data=pilot))
}
\keyword{datasets}
