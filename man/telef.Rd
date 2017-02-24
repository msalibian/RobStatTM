\name{telef}
\alias{telef}
\docType{data}
\title{Number of International Calls from Belgium}
\description{
  Number of international calls from Belgium, taken from the Belgian
  Statistical Survey, published by the Ministry of Economy.
}
\usage{data(telef)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{Calls}}{Number of Calls (in tens of millions) }
    \item{\code{Year}}{Year (1950 - 1973)}
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, page 26, table 2.
}
\examples{
data(telef)
summary(lm.telef <- lm(Year~., data=telef))
%%summary(lts.telef <- ltsReg(Year~., data=telef))
}
\keyword{datasets}
