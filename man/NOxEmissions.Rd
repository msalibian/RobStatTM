\name{NOxEmissions}
\alias{NOxEmissions}
\docType{data}
\encoding{utf8}
\title{NOx Air Pollution Data}
\description{
  A typical medium sized environmental data set with hourly measurements
  of \eqn{NOx} pollution content in the ambient air.
}
\usage{data(NOxEmissions)}
\format{
  A data frame with 8088 observations on the following 4 variables.
  \describe{
    \item{\code{julday}}{day number, a factor with levels \code{373}
      \dots \code{730}, typically with 24 hourly measurements.}
    \item{\code{LNOx}}{\eqn{\log} of hourly mean of NOx concentration in
      ambient air [ppb] next to a highly frequented motorway.}
    \item{\code{LNOxEm}}{\eqn{\log} of hourly sum of NOx emission of
      cars on this motorway in arbitrary units.}
    \item{\code{sqrtWS}}{Square root of wind speed [m/s].}
  }
}
\details{
  The original data set had more observations, but with missing values.
  Here, all cases with missing values were omitted
  (\code{\link{na.omit}(.)}), and then only those were retained that
  belonged to days with at least 20 (fully) observed hourly
  measurements.
}
\source{
  René Locher (at ZHAW, Switzerland).
%% E-mail to R-SIG-robust mailing list, on 2006-04-20.
}
\seealso{another NOx dataset, \code{\link{ambientNOxCH}}.
}
% \references{
%   ~~ possibly secondary sources and usages ~~
% }
\examples{
data(NOxEmissions)
plot(LNOx ~ LNOxEm, data = NOxEmissions, cex = 0.25, col = "gray30")

\dontrun{## these take too much time --
 ## p = 340  ==> already Least Squares is not fast
 (lmNOx <- lm(LNOx ~ . ,data = NOxEmissions))
 plot(lmNOx) #->  indication of 1 outlier

 M.NOx <- MASS::rlm(LNOx ~ . , data = NOxEmissions)
 ## M-estimation works
 ## whereas  MM-estimation fails:
 try(MM.NOx <- MASS::rlm(LNOx ~ . , data = NOxEmissions, method = "MM"))
 ## namely because S-estimation fails:
 try(lts.NOx <- ltsReg(LNOx ~ . , data = NOxEmissions))
 try(lmR.NOx <- lmrob (LNOx ~ . , data = NOxEmissions))
}% don't run
}
\keyword{datasets}
