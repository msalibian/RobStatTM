\name{lactic}
\alias{lactic}
%% FIXME: shorter
\title{Lactic Acid Concentration Measurement Data}
\description{
  Data on the Calibration of an Instrument that Measures
  Lactic Acid Concentration in Blood, from Afifi and Azen (1979) -
  comparing the true concentration X with the measured value Y.
}
\usage{data(lactic)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{X}}{True Concentration}
    \item{\code{Y}}{Instrument}
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, p.62, table 10.
}
\examples{
data(lactic)
summary(lm.lactic <- lm(Y ~., data=lactic))
%% summary(lts.lactic <- ltsReg(Y ~., data=lactic))
}
\keyword{datasets}
