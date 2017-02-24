\name{pension}
\alias{pension}
\title{Pension Funds Data}
\description{
 The total 1981 premium income of pension funds of Dutch firms,
 for 18 Professional Branches, from de Wit (1982).
}
\usage{data(pension)}
\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{\code{Income}}{Premium Income (in millions of guilders)}
    \item{\code{Reserves}}{Premium Reserves (in millions of guilders)}
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, p.76, table 13.
}
\examples{
data(pension)
plot(pension)

summary(lm.p  <-    lm(Reserves ~., data=pension))
summary(lmR.p <- lmrob(Reserves ~., data=pension))
summary(lts.p <- ltsReg(Reserves ~., data=pension))
abline( lm.p)
abline(lmR.p, col=2)
abline(lts.p, col=2, lty=2)

## MM: "the" solution is much simpler:
plot(pension, log = "xy")
lm.lp  <-    lm(log(Reserves) ~ log(Income), data=pension)
lmR.lp <- lmrob(log(Reserves) ~ log(Income), data=pension)
plot(log(Reserves) ~ log(Income), data=pension)
## no difference between LS and robust:
abline( lm.lp)
abline(lmR.lp, col=2)
}
\keyword{datasets}
