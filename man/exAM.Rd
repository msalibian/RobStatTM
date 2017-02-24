\name{exAM}
\alias{exAM}
\docType{data}
\title{Example Data of Antille and May - for Simple Regression}
\description{
  This is an artificial data set, cleverly construced and used by
  Antille and May to demonstrate \sQuote{problems} with LMS and LTS.
}
\usage{data(exAM)}
\format{
  A data frame with 12 observations on 2 variables, \code{x} and \code{y}.
}
\details{
  Because the points are not in general position, both LMS and LTS
  typically \emph{fail}; however, e.g., \code{\link[MASS]{rlm}(*,
    method="MM")} \dQuote{works}.
}
\source{
  Antille, G. and El May, H. (1992)
  The use of slices in the LMS and the method of density slices:
  Foundation and comparison.\cr
  In Yadolah Dodge and Joe Whittaker, editors, \emph{COMPSTAT:
  Proc. 10th Symp. Computat. Statist., Neuchatel}, \bold{1}, 441--445;
  Physica-Verlag.
}
% \references{
%   ~~ possibly secondary sources and usages ~~
% }
\examples{
data(exAM)
plot(exAM)
summary(ls <- lm(y ~ x, data=exAM))
abline(ls)
}
\keyword{robust}
\keyword{datasets}
