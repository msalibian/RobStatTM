\name{alcohol}
\alias{alcohol}
\docType{data}
\title{Alcohol Solubility in Water Data}
\description{
  The solubility of alcohols in water is important in understanding
  alcohol transport in living organisms.  This dataset from (Romanelli
  et al., 2001) contains physicochemical characteristics of 44 aliphatic
  alcohols.  The aim of the experiment was the prediction of the
  solubility on the basis of molecular descriptors.
}
\usage{data(alcohol)}
\format{
  A data frame with 44 observations on the following 7 numeric variables.
  \describe{
    \item{\code{SAG}}{solvent accessible surface-bounded molecular volume.}
    \item{\code{V}}{volume}
    \item{\code{logPC}}{Log(PC); PC = octanol-water partitions coefficient}
    \item{\code{P}}{polarizability}
    \item{\code{RM}}{molar refractivity}
    \item{\code{Mass}}{the mass}
    \item{\code{logSolubility}}{ln(Solubility), the response.}
  }
}
% \details{
% }
\source{
  The website accompanying the MMY-book:
  \url{http://www.wiley.com/legacy/wileychi/robust_statistics}
}
\references{
  %% FIXME: I don't think this is the correct reference
  %   Romanelli, J.R., Kelly, J.J. and Litwein, D.E.M (2001)
  %   Hand-assisted laparoscopic surgery in the United States: An overview
  %   \emph{Seminars in Laparoscopic Surgery} \bold{8} 96--103.

  Maronna, R.A., Martin, R.D. and Yohai, V.J. (2006)
  \emph{Robust Statistics, Theory and Methods}, Wiley.
}
\examples{
data(alcohol)
## version of data set with trivial names, as
s.alcohol <- alcohol
names(s.alcohol) <- paste("Col", 1:7, sep='')
}
\keyword{datasets}
