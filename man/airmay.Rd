\name{airmay}
\alias{airmay}
\docType{data}
\title{Air Quality Data}
\description{
   Air Quality Data Set for May 1973, from Chambers et al. (1983).
   The whole data set consists of daily readings of air quality
   values from May 1, 1973 to September 30, 1973,
   but here are included only the values for May. This data set
   is an example of the special treatment of the missing values.
}
\usage{data(airmay)}
\format{
  A data frame with 31 observations on the following 4 variables.
  \describe{
    \item{\code{X1}}{Solar Radiation in Longleys in the frequency
      band 4000-7700 from 0800 to 1200 hours at Central Park}
    \item{\code{X2}}{Average windspeed (in miles per hour) between 7000
      and 1000 hours at La Guardia Airport}
    \item{\code{X3}}{Maximum daily temperature (in degrees Fahrenheit)
      at La Guardia Airport}
    \item{\code{Y}}{Mean ozone concentration (in parts per billion)
      from 1300 to 1500 hours at Roosevelt Island}
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, p.86, table 6.
}
\examples{
data(airmay)
summary(lm.airmay <- lm(Y ~ ., data=airmay))
%%summary(lts.airmay <- ltsReg(Y ~ ., data=airmay))

airmay.x <- data.matrix(airmay[,1:3])
%%covMcd(airmay.x)
}
\keyword{datasets}
