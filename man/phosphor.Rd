\name{phosphor}
\alias{phosphor}
\docType{data}
\title{Phosphorus Content Data}
\description{
  This dataset investigates the effect from inorganic and organic
  Phosphorus in the soil upon the phosphorus content of the corn grown
  in this soil, from Prescott (1975).
}
\usage{data(phosphor)}
\format{
  A data frame with 18 observations on the following 3 variables.
  \describe{
    \item{\code{inorg}}{Inorganic soil Phosphorus}
    \item{\code{organic}}{Organic soil Phosphorus}
    \item{\code{plant}}{Plant Phosphorus content}
  }
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987)
 \emph{Robust Regression and Outlier Detection.} Wiley, p.156, table 24.
}
\examples{
data(phosphor)
plot(phosphor)
summary(lm.phosphor <- lm(plant ~ ., data = phosphor))
summary(lts.phosphor <- ltsReg(plant ~ ., data = phosphor))

phosphor.x <- data.matrix(phosphor[, 1:2])
cPh <- covMcd(phosphor.x)
plot(cPh, "dd")
}
\keyword{datasets}
