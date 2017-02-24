\name{heart}
\alias{heart}
\docType{data}
\title{Heart Catherization Data}
\description{
  This data set was analyzed by Weisberg (1980) and Chambers et
  al. (1983).  A catheter is passed into a major vein or artery at the
  femoral region and moved into the heart.  The proper length of the
  introduced catheter has to be guessed by the physician. The aim of the
  data set is to describe the relation between the catheter length and
  the patient's height (X1) and weight (X2).

  This data sets is used to demonstrate the effects caused by collinearity.
  The correlation between height and weight is so high that either
  variable almost completely determines the other.
}
\usage{
data(heart)
%> QA bug: would want:
%> data(heart, package="robustbase")
%> but that gives two warnings
}
\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{height}}{Patient's height in inches}
    \item{\code{weight}}{Patient's weights in pounds}
    \item{\code{clength}}{Y: Catheter Length (in centimeters)}
  }
}
\note{There are other \code{heart} datasets in other \R packages,
  notably \pkg{survival}, hence considering using
  \code{package = "robustbase"}, see examples.
}
\source{
  Weisberg (1980)

  Chambers et al. (1983)

  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, p.103, table 13.
}
\examples{
data(heart, package="robustbase")
heart.x <- data.matrix(heart[, 1:2]) # the X-variables
plot(heart.x)
covMcd(heart.x)
summary( lm.heart <-     lm(clength ~ . , data = heart))
summary(lts.heart <- ltsReg(clength ~ . , data = heart))
}
\keyword{datasets}
